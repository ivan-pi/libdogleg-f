!> lmder example from MINPACK
!>
!> This example is taken from the documentation of the function `lmder` in 
!> the User Guide for MINPACK-1 (pgs. 79-85):
!>
!>    The problem is to determine the values of x(1), x(2), and x(3) which 
!>    provide the best fit (in the least squares sense) of
!>       
!>       x(1) + u(i)/(v(i)*x(2) + w(i)*x(3)),   i = 1, 15
!>
!>    to the data
!>
!>       y = [0.14,0.18,0.22,0.25,0.29,0.32,0.35,0.39,
!>            0.37,0.58,0.73,0.96,1.34,2.10,4.39]
!>
!>    where u(i) = i, v(i) = 16 - i, and w(i) = min(u(i),v(i)). The
!>    j-th component of FVEC is thus defined by
!>
!>       y(i) - (x(1) + u(i)/(v(i)*x(2) + w(i)*x(3)))
!>
!> Reference:
!>   J. J. Moré, B. S. Garbow, and K. E. Hillstrom, User Guide for MINPACK-1, 
!>   Argonne National Laboratory Report ANL-80-74, Argonne, Ill., 1980.
!>   https://cds.cern.ch/record/126569/files/CM-P00068642.pdf
!>   
module lmder_fcn

   implicit none
   private

   public :: dp, nstate, nres, dl_fcn

   integer, parameter :: dp = kind(1.0d0)

   integer, parameter :: nstate = 3
      ! Number of state parameters

   integer, parameter :: nres = 15
      ! Number of residuals, equal to the number of data points

   real(dp), parameter :: y(nres) = &
      [1.4d-1,1.8d-1,2.2d-1,2.5d-1,2.9d-1,3.2d-1,3.5d-1,3.9d-1, &
      3.7d-1, 5.8d-1,7.3d-1,9.6d-1,1.34d0,2.1d0,4.39d0]

contains

   subroutine fcn(m,n,x,fvec,fjac,iflag)
      integer, intent(in) :: m, n, iflag
      real(dp), intent(in) :: x(n)
      real(dp), intent(inout) :: fvec(m)
      real(dp), intent(inout) :: fjac(n, m)
         ! note: the fjac array has been transposed compared to the 
         !       original example from MINPACK; the storage pattern used by
         !       MINPACK is opposite of what libdogleg expects

      integer :: i
      real(dp) :: tmp1, tmp2, tmp3, tmp4

      if (iflag == 1) then
         do i = 1, 15
            tmp1 = i
            tmp2 = 16 - i
            tmp3 = tmp1
            if (i > 8) tmp3 = tmp2
            fvec(i) = y(i) - (x(1) + tmp1/(x(2)*tmp2 + x(3)*tmp3))
         end do
      else if (iflag == 2) then
         do i = 1, 15
            tmp1 = i
            tmp2 = 16 - i
            tmp3 = tmp1
            if (i > 8) tmp3 = tmp2
            tmp4 = (x(2)*tmp2 + x(3)*tmp3)**2
            fjac(1,i) = -1.0_dp
            fjac(2,i) = tmp1*tmp2/tmp4
            fjac(3,i) = tmp1*tmp3/tmp4
         end do
      end if
   end subroutine

   subroutine dl_fcn(p, x, J, params)
      real(dp), intent(in) :: p(:)     ! shape(p) := [nstate]
      real(dp), intent(out) :: x(:)    ! shape(x) := [nres]
      real(dp), intent(out) :: J(:,:)  ! shape(J) := [nstate, nres]
      class(*), intent(in), optional :: params

      integer :: m, n

      m = size(x)   ! nmeas or nresiduals
      n = size(p)   ! nstate

      !
      ! Evaluate objective function
      !
      call fcn(m, n, p, x, J, iflag = 1)

      !
      ! Evaluate Jacobian
      !
      call fcn(m, n, p, x, J, iflag = 2)

   end subroutine

end module


program lmder_fcn_main

   use dogleg, only: dl_parameters2_t, dl_getDefaultParameters, dl_optimize
   use lmder_fcn, only: dp, dl_fcn, nstate, nres

   implicit none

   real(dp) :: p(nstate), optf
   type(dl_parameters2_t) :: settings

   p = [1.0_dp, 1.0_dp, 1.0_dp]

   call dl_getDefaultParameters(settings)
   settings%Jt_x_threshold        = sqrt(epsilon(optf))
   settings%update_threshold      = sqrt(epsilon(optf))
   settings%trustregion_threshold = sqrt(epsilon(optf))

   write(*,'(A)') "Solving..."
   call print_state("initial:", p)

   call dl_optimize(optf, p, nres, dl_fcn, &
      parameters = settings)

   write(*, '(A, G0)') "Done. Optimum = ", optf
   call print_state("final:", p)

   ! Values given in User Guide for MINPACK
   p = [0.8241058d-01, 0.1133037d+01, 0.2343695d+01]
   call print_state("expected:", p)

contains

   !> Helper function to print results
   subroutine print_state(key, p)
      character(len=*), intent(in) :: key
      real(dp), intent(in) :: p(:)
      integer :: i
      write(*, '(A)') key
      do i = 1, size(p)
         write(*, '("  p(",I0,") = ",G12.7)')  i, p(i)
      end do
   end subroutine

end program