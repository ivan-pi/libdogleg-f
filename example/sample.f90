module sample
   
   use, intrinsic :: iso_c_binding, only: c_double
   implicit none

   integer, parameter :: wp = c_double

!  This is a trivial sample application to demonstrate libdogleg in action.
!  Let's say that I have a simple non-linear model
! 
!  a*b * x**2 + b*c * y**2 + c * x*y + d * x + e * y + f = measurements
! 
!  here I'm trying to estimate the vector (a,b,c,d,e,f) to most closely fit the
!  data vector measurements. This problem is clearly non-sparse, but both sparse
!  and dense versions of libdogleg are demonstrated here.
! 
!  First I generate some noise-corrupted data, and then use libdogleg to solve
!  the problem.

!  My state vector (a,b,c,d,e,f) has 6 elements
   integer, parameter :: nstate = 6

!  I simulate measurements using these as the true values
   real(wp), parameter :: ref_a = 1.0_wp
   real(wp), parameter :: ref_b = 2.0_wp
   real(wp), parameter :: ref_c = 3.0_wp
   real(wp), parameter :: ref_d = 4.0_wp
   real(wp), parameter :: ref_e = 5.0_wp
   real(wp), parameter :: ref_f = 6.0_wp

!  I simulate by sampling x-y space in a grid. This grid is defined here
   integer, parameter :: grid_width = 10
   real(wp), parameter :: grid_min = -10.0_wp
   real(wp), parameter :: grid_delta = 2.0_wp
   integer, parameter :: nmeas = grid_width**2

   real(wp) :: allx(nmeas)
   real(wp) :: ally(nmeas)
   real(wp) :: allm_simulated_noisy(nmeas)

contains

   subroutine generate_simulation_grid()

      integer :: j, ix, iy
      real(wp) :: x, y

      j = 1
      do ix = 1, grid_width
         x = grid_min + (ix - 1)*grid_delta
         do iy = 1, grid_width
            y = grid_min + (iy - 1)*grid_delta
            allx(j) = x
            ally(j) = y
            j = j + 1
         end do
      end do

   end subroutine

   subroutine simulate()

      associate( &
         a => ref_a, &
         b => ref_b, &
         c => ref_c, &
         d => ref_d, &
         e => ref_e, &
         f => ref_f, &
         x => allx, &
         y => ally)
         
         allm_simulated_noisy = a*b*x**2 + b*c*y**2 + c*x*y + d*x + e*y + f

      end associate
      
      ! add +- 0.5 units of uniformly-random noise
      add_noise: block
         real(wp) :: noise(nmeas)
         call random_number(noise)
         allm_simulated_noisy = allm_simulated_noisy + (noise - 0.5_wp)
      end block add_noise

   end subroutine


   subroutine optimizer_callback(p, x, J, params)
      real(wp), intent(in) :: p(:)
      real(wp), intent(out) :: x(:)
      real(wp), intent(out) :: J(:,:)
      class(*), intent(in), optional :: params

      integer :: i

      associate( &
         a => p(1), &
         b => p(2), &
         c => p(3), &
         d => p(4), &
         e => p(5), &
         f => p(6))

         do i = 1, size(x)
            x(i) = &
               a*b*allx(i)**2    + &
               b*c*ally(i)**2    + &
               c*allx(i)*ally(i) + &
               d*allx(i)         + &
               e*ally(i)         + &
               f                 - &
               allm_simulated_noisy(i)

            J(1,i) = b*allx(i)**2
            J(2,i) = a*allx(i)**2 + c*ally(i)**2
            J(3,i) = b*ally(i)**2 + allx(i)*ally(i)
            J(4,i) = allx(i)
            J(5,i) = ally(i)
            J(6,i) = 1.0_wp 

         end do

      end associate

   end subroutine

end module

program sample_main

   use, intrinsic :: iso_fortran_env, only: stderr => error_unit
   
   use dogleg

   use sample, only: nstate, nmeas, wp, &
      generate_simulation_grid, simulate, optimizer_callback

   implicit none

   real(wp) :: p(nstate), optimum
   type(dl_parameters2_t) :: dl_params
   
   integer :: i, nseed
   integer, allocatable :: seed(:)
   

   ! We want a deterministic seed
   call random_seed(size = nseed)
   allocate(seed(nseed))
   seed = 0
   call random_seed(put = seed)


   call generate_simulation_grid()
   call simulate()

!  I start solving with all my state variables set to some random noise
   call random_number(p)
   p = p - 0.1

   call dl_getDefaultParameters(dl_params)
!   dl_params%dogleg_debug = 1
   call display_dl_params(dl_params)

   write(stderr,'(A)') "SOLVING:"
   call dl_optimize(optimum, p, nmeas, optimizer_callback, &
      parameters = dl_params)

   write(stderr, *) "Done. Optimum = ", optimum
   write(stderr, *) "optimal state:"
   do i = 1, size(p)
      write(stderr, '("  p(",I0,") = ",G0)')  i, p(i)
   end do

contains

   subroutine display_dl_params(p)
      type(dl_parameters2_t), intent(in) :: p

      print *, "p%max_iterations                 = ", p%max_iterations
      print *, "p%dogleg_debug                   = ", p%dogleg_debug
      print *, "p%trustregion0                   = ", p%trustregion0
      print *, "p%trustregion_decrease_factor    = ", p%trustregion_decrease_factor
      print *, "p%trustregion_decrease_threshold = ", p%trustregion_decrease_threshold
      print *, "p%trustregion_increase_factor    = ", p%trustregion_increase_factor
      print *, "p%trustregion_increase_threshold = ", p%trustregion_increase_threshold
      print *, "p%Jt_x_threshold                 = ", p%Jt_x_threshold
      print *, "p%update_threshold               = ", p%update_threshold
      print *, "p%trustregion_threshold          = ", p%trustregion_threshold

   end subroutine

end program