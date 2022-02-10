!> Powell's function example
!>
!> A full description can be found in the Ceres library:
!>   http://ceres-solver.org/nnls_tutorial.html#powell-s-function
!>
!> The full code for the original example can be found at:
!>   https://ceres-solver.googlesource.com/ceres-solver/+/master/examples/powell.cc
!>
module powell

   use, intrinsic :: iso_c_binding, only: wp => c_double
   implicit none

contains

   !> This example considers the optimization of Powell's function
   !> where p = [x1, x2, x3, x4] (nstate = 4), and
   !>
   !>   f1(p) = x1 + 10*x2
   !>   f2(p) = sqrt(5) * (x3 - x4)
   !>   f3(p) = (x2 - 2*x3)**2
   !>   f4(p) = sqrt(10) * (x1 - x4)**2
   !>
   !> The objective is to minimize 0.5*|F(p)| (nmeas = 4) where
   !>    
   !>    F(p) = [f1(p), f2(p), f3(p), f4(p)]
   !>
   !> The exact solution to this problem is p = [0, 0, 0, 0].
   !>
   subroutine powells_function(p, x, J, params)
      real(wp), intent(in) :: p(:)     ! shape(p) := [nstate]
      real(wp), intent(out) :: x(:)    ! shape(x) := [nmeas]
      real(wp), intent(out) :: J(:,:)  ! shape(J) := [nstate, nmeas]
      class(*), intent(in), optional :: params

      real(wp) :: f1, f2, f3, f4

      associate( &
         ! Unpack the optimization variables
         x1 => p(1), &
         x2 => p(2), &
         x3 => p(3), &
         x4 => p(4))

         !
         ! Evaluate objective function
         !
         f1 = x1 + 10*x2
         f2 = sqrt(5.0_wp)*(x3 - x4)
         f3 = (x2 - 2*x3)**2
         f4 = sqrt(10.0_wp)*(x1 - x4)**2
         x = [f1, f2, f3, f4]

         !
         ! Evaluate the dense Jacobian
         !

         ! Derivatives of f1 with respect to [x1, x2, x3, x4]
         J(1,1) = 1.0_wp
         J(2,1) = 10.0_wp
         J(3,1) = 0.0_wp 
         J(4,1) = 0.0_wp
         
         ! Derivatives of f2 with respect to [x1, x2, x3, x4]
         J(1,2) = 0.0_wp
         J(2,2) = 0.0_wp
         J(3,2) = sqrt(5.0_wp) 
         J(4,2) = -J(3,2)
         
         ! Derivatives of f3 with respect to [x1, x2, x3, x4]
         J(1,3) = 0.0_wp
         J(2,3) = 2.0_wp*(x2 - 2*x3)
         J(3,3) = -2.0_wp*J(2,3)
         J(4,3) = 0.0_wp
         
         ! Derivatives of f4 with respect to [x1, x2, x3, x4]
         J(1,4) = 2.0_wp*sqrt(10.0_wp)*(x1 - x4)
         J(2,4) = 0.0_wp
         J(3,4) = 0.0_wp 
         J(4,4) = -J(1,4)

      end associate

   end subroutine


end module

program powell_main

   use dogleg, only: dl_parameters2_t, dl_getDefaultParameters, dl_optimize
   use powell, only: wp, powells_function

   implicit none

   integer, parameter :: nstate = 4, nmeas = 4
   real(wp) :: p(nstate), optf
   type(dl_parameters2_t) :: settings

   p = [3.0_wp, -1.0_wp, 0.0_wp, 1.0_wp]

   call dl_getDefaultParameters(settings)
   settings%Jt_x_threshold        = 1.e-10_wp
   settings%update_threshold      = 1.e-10_wp
   settings%trustregion_threshold = 1.e-10_wp

   write(*,'(A)') "Solving..."
   call print_state("initial:", p)

   call dl_optimize(optf, p, nmeas, powells_function, &
      parameters = settings)

   write(*, '(A, G0)') "Done. Optimum = ", optf
   call print_state("final:", p)

contains

   !> Helper function to print results
   subroutine print_state(key, p)
      character(len=*), intent(in) :: key
      real(wp), intent(in) :: p(:)
      integer :: i
      write(*, '(A)') key
      do i = 1, size(p)
         write(*, '("  p(",I0,") = ",G0)')  i, p(i)
      end do
   end subroutine

end program