!> MathWorks Example
!>
!> This example is borrowed from the documentation of the MATLAB function 
!> `lsqnonlin`. The example appears under the title "Nonlinear Least Squares 
!> Solution and Residual Norm"
!>
!> The goal is to find the pair of values x = (x1, x2) that minimize
!>
!>    \sum_{k=1}^10 (2 + 2*k - exp(k*x1) - exp(k*x2))^2
!>
!> Note the residual is composed of 10 functions of the form
!>
!>    F_k(x) = 2 + 2*k - exp(k*x1) - exp(k*x2)
!>
!> These are to be returned in the output vector of the callback function.
!>
module mathworks_example

   use, intrinsic :: iso_c_binding, only: wp => c_double

   implicit none
   private

   public :: wp, mw_fcn, nstate, nres

   integer, parameter :: nstate = 2
   integer, parameter :: nres   = 10

contains

   subroutine mw_fcn(x, F, J, params)
      real(wp), intent(in) :: x(:)     ! shape(x) := [nstate]
      real(wp), intent(out) :: F(:)    ! shape(F) := [nres]
      real(wp), intent(out) :: J(:,:)  ! shape(J) := [nstate, nres]
      class(*), intent(in), optional :: params
      integer :: k
      associate(x1 => x(1), x2 => x(2))
         do k = 1, nres
            F(k) = 2 + 2*k - exp(k*x1) - exp(k*x2)
            J(1,k) = -exp(k*x1)*k
            J(2,k) = -exp(k*x2)*k
         end do
      end associate
   end subroutine

end module

program mathworks_example_main

   use dogleg, only: dl_parameters2_t, dl_getDefaultParameters, dl_optimize
   use mathworks_example, only: wp, mw_fcn, nstate, nres

   implicit none

   real(wp) :: x(nstate), optf
   type(dl_parameters2_t) :: settings

   ! initial guess
   x = [0.3_wp, 0.4_wp]

   call dl_getDefaultParameters(settings)

   write(*,'(A)') "Solving..."
   call print_state("initial:", x)

   call dl_optimize(optf, x, nres, mw_fcn, &
      parameters = settings)

   write(*, '(A, F9.4)') "Done. Optimum = ", optf
   call print_state("final:", x)

   ! Values given in MathWorks documentation (not to full precision)
   x = [0.2578_wp, 0.2578_wp]
   call print_state("expected:", x)

contains

   !> Helper function to print results
   subroutine print_state(key, x)
      character(len=*), intent(in) :: key
      real(wp), intent(in) :: x(:)
      integer :: i
      write(*, '(A)') key
      do i = 1, size(x)
         write(*, '("  x(",I0,") = ",G9.4)')  i, x(i)
      end do
   end subroutine

end program