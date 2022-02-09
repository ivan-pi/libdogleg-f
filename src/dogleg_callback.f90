module dogleg_callback

   use, intrinsic :: iso_c_binding, only: &
      c_double, c_int, c_f_pointer, c_ptr
   
   implicit none

   abstract interface
      subroutine dl_callback_dense(p, x, J, params)
         import c_double
         implicit none
         real(c_double), intent(in) :: p(:)
         real(c_double), intent(out) :: x(:)
         real(c_double), intent(out) :: J(:,:)
         class(*), intent(in), optional :: params
      end subroutine
   end interface

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

!> Overloads the structure constructor for the callback wrapper type
!>
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
!> expected by libdogleg. Both the callback function and the
!> data are passed through the `void* cookie`
!>
!> The callback prototype expected by libdogleg is:
!>
!>   void (dogleg_callback_dense_t)(const double* p, 
!>                                  double*       x, 
!>                                  double*       J, 
!>                                  void*         cookie)
!>
   subroutine dl_callback_adaptor(p, x, J, cookie) bind(c)
      real(c_double), intent(in), target :: p(*)
         !> p is the current state vector
      real(c_double), intent(inout) :: x(*)
      real(c_double), intent(inout), target :: J(*)
      type(c_ptr), value :: cookie

      type(dl_dense_t), pointer :: cb 
         ! Callback container

      real(c_double), pointer :: J_f(:,:) => null()
         ! Jacobian matrix as Fortran rank-2 array

      ! Retrieve Fortran callback object from c_ptr
      call c_f_pointer(cookie, cb)

      associate(nstate => cb%nstate, nmeas => cb%nmeas)

         ! Pointer remapping of Jacobian to rank-2 matrix
         J_f(1:nstate,1:nmeas) => J(1:nstate*nmeas)

         ! Execute Fortran callback function
         call cb%f(p(1:nstate), x(1:nmeas), J_f, cb%fparams)

      end associate

   end subroutine 

end module