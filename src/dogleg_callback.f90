!
! This file is part of libdogleg-f.
!
! libdogleg-f is free software: you can redistribute it and/or modify it under 
! the terms of the GNU Lesser General Public License as published by the Free 
! Software Foundation, either version 3 of the License, or (at your option) 
! any later version.
!
! libdogleg-f is distributed in the hope that it will be useful, but WITHOUT 
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
! FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License 
! for more details.
!
! You should have received a copy of the GNU Lesser General Public License 
! along with libdogleg-f. If not, see <https://www.gnu.org/licenses/>. 
!
module dogleg_callback

   use, intrinsic :: iso_c_binding, only: &
      c_double, c_int, c_f_pointer, c_ptr
   
   implicit none

   abstract interface
      !> Fortran callback interface for dense problems
      !>
      !> Evaluates the vector function x = f(p) and jacobian J = df/dp
      !>
      !> Arguments
      !>
      !>   p (in): an input array of size nstate
      !>   x (out): an output array of size nmeas, containing the function
      !>            values x = f(p)
      !>   J (out): a rank-2 array of shape [nstate, nmeas] containing the Jacobian
      !>            matrix of function f with respect to the parameters p
      !>   params (in, optional): a container for any additional callback parameters 
      !>                          which are not part of the optimization problem
      !>
      subroutine dl_callback_dense(p, x, J, params)
         import c_double
         implicit none
         real(c_double), intent(in) :: p(:)     ! shape(p) := [nstate]
         real(c_double), intent(out) :: x(:)    ! shape(x) := [nmeas]
         real(c_double), intent(out) :: J(:,:)  ! shape(J) := [nstate, nmeas]
         class(*), intent(in), optional :: params
      end subroutine
   end interface

   !> Dense callback container.
   !>
   !> Used to pass the callback function and parameters over to the C library.
   !>
   type :: dl_dense_t
      private
      integer(c_int) :: nstate = 0, nmeas = 0
      procedure(dl_callback_dense), nopass, pointer :: f => null()
      class(*), pointer :: fparams => null()
   end type

   interface dl_dense_t
      module procedure :: construct_dogleg_dense_t
   end interface

contains

!> Overloaded structure constructor for the dense callback container type
   function construct_dogleg_dense_t(nstate, nmeas, f, fparams) result(this)
      integer(c_int), intent(in) :: nstate, nmeas
      procedure(dl_callback_dense) :: f
      class(*), intent(in), target, optional :: fparams
      type(dl_dense_t) :: this

      this%nstate = nstate
      this%nmeas = nmeas

      this%f => f
      nullify(this%fparams)
      if (present(fparams)) &
         this%fparams => fparams
   end function

!> We use an adaptor routine to conform to the interface
!> expected by libdogleg. The actual Fortran callback function and any
!> additional parameters are passed through the `void*` dummy argument `cookie`.
!>
!> The callback function prototype expected by `libdogleg` is:
!>
!>```C
!>   void (dogleg_callback_dense_t)(const double* p, 
!>                                  double*       x, 
!>                                  double*       J, 
!>                                  void*         cookie)
!>```
   subroutine dl_callback_adaptor(p, x, J, cookie) bind(c)
      real(c_double), intent(in), target :: p(*)
         !> p is the current state vector
      real(c_double), intent(inout) :: x(*)
      real(c_double), intent(inout), target :: J(*)
      type(c_ptr), value :: cookie

      type(dl_dense_t), pointer :: cb 
      real(c_double), pointer :: Jf(:,:) => null()

      ! Retrieve Fortran callback object from c_ptr
      call c_f_pointer(cookie, cb)

      associate(nstate => cb%nstate, nmeas => cb%nmeas)

         ! Pointer remapping of the Jacobian to a rank-2 array (matrix)
         !
         ! Luckily, libdogleg uses a column-first ordering of the Jacobian
         ! which is well-suited to Fortran
         !
         Jf(1:nstate,1:nmeas) => J(1:nstate*nmeas)

         !
         ! Call the actual Fortran callback function
         !
         call cb%f(p(1:nstate), x(1:nmeas), Jf, cb%fparams)

      end associate

   end subroutine 

end module