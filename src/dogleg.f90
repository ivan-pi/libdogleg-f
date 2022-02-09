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
module dogleg

    use iso_fortran_env, only: int64

    use, intrinsic :: iso_c_binding, only: &
        wp => c_double, c_int, c_ptr

    use dogleg_callback, only: &
        dl_dense_t, dl_callback_dense, dl_callback_adaptor
    use dogleg_c, dl_parameters2_t => dogleg_parameters2_t

    implicit none
    private

    !
    ! interfaces
    !
    public :: dl_callback_dense

    !
    ! types
    !
    public :: dl_parameters2_t
    public :: dl_solver_context_t

    !
    ! procedures
    !
    public :: dl_getDefaultParameters
    public :: dl_setMaxIterations
    public :: dl_setTrustregionUpdateParameters
    public :: dl_setDebug
    public :: dl_setInitialTrustregion
    public :: dl_setThresholds
    public :: dl_optimize
    public :: dl_testGradient
    public :: dl_freeContext

! -------------------------------------------------------


    interface dl_optimize
        module procedure :: dl_optimize_dense
        module procedure :: dl_optimize_dense2
    end interface


    interface dl_testGradient
        module procedure :: dl_testGradient_dense
    end interface


    type :: dl_solver_context_t
        private
        type(c_ptr) :: handle = c_null_ptr
    contains
        final :: dl_freeContext
    end type

contains

    subroutine dl_getDefaultParameters(parameters)
        type(dl_parameters2_t), intent(out), target :: parameters
        call dogleg_getDefaultParameters(c_loc(parameters))
    end subroutine

    subroutine dl_setMaxIterations(n)
        integer, intent(in) :: n
        call dogleg_setMaxIterations(int(n, c_int))
    end subroutine  

    subroutine dl_setTrustregionUpdateParameters( &
        downFactor, downThreshold, upFactor, upThreshold)
        real(wp), intent(in) :: downFactor, downThreshold
        real(wp), intent(in) :: upFactor, upThreshold
        call dogleg_setTrustregionUpdateParameters( &
            downFactor, downThreshold, upFactor, upThreshold)
    end subroutine

    subroutine dl_setDebug(debug)
        integer, intent(in) :: debug
        call dogleg_setDebug(int(debug,c_int))
    end subroutine

    subroutine dl_setInitialTrustregion(t)
        real(wp), intent(in) :: t
        call dogleg_setInitialTrustregion(t)
    end subroutine

    subroutine dl_setThresholds(Jt_x, update, trustregion)
        real(wp), intent(in) :: Jt_x, update, trustregion
        call dogleg_setThresholds(Jt_x, update, trustregion)
    end subroutine

    subroutine dl_optimize_dense(optimum, p, nmeas, f, fparams, ctx)
        real(wp), intent(out) :: optimum
        real(wp), intent(inout), contiguous :: p(:)
        integer, intent(in) :: nmeas
        procedure(dl_callback_dense) :: f
        class(*), intent(in), target, optional :: fparams
        type(dl_solver_context_t), intent(out), optional :: ctx

        type(dl_dense_t), target :: cb
        integer(c_int) :: nstate_, nmeas_
        type(c_ptr) :: ctx_handle

        nstate_ = size(p, kind=c_int)
        nmeas_ = int(nmeas, kind=c_int)

        cb = dl_dense_t(nstate_, nmeas_, f, fparams)

        if (present(ctx)) then
            ctx_handle = ctx%handle
        else
            ctx_handle = c_null_ptr
        end if


        optimum = dogleg_optimize_dense(p, nstate_, nmeas_, &
            f = c_funloc(dl_callback_adaptor), &
            cookie = c_loc(cb), &
            returnContext = ctx_handle) 

    end subroutine 

    subroutine dl_optimize_dense2(optimum, p, nmeas, f, fparams, &
            parameters, ctx)
        
        real(wp), intent(out) :: optimum
        real(wp), intent(inout), contiguous :: p(:)
        integer, intent(in) :: nmeas
        procedure(dl_callback_dense) :: f
        class(*), intent(in), target, optional :: fparams
        type(dl_parameters2_t), intent(in), target :: parameters
        type(dl_solver_context_t), intent(out), optional :: ctx

        type(dl_dense_t), target :: cb
        integer(c_int) :: nstate_, nmeas_
        type(c_ptr) :: ctx_handle
        integer(int64) :: i64

        nstate_ = size(p, kind=c_int)
        nmeas_ = int(nmeas, kind=c_int)

        cb = dl_dense_t(nstate_, nmeas_, f, fparams)

        if (present(ctx)) then
            ctx_handle = ctx%handle
        else
            ctx_handle = c_null_ptr
        end if

        optimum = dogleg_optimize_dense2(p, nstate_, nmeas_, &
            f = c_funloc(dl_callback_adaptor), &
            cookie = c_loc(cb), &
            parameters = c_loc(parameters), &
            returnContext = ctx_handle) 

    end subroutine

    subroutine dl_testGradient_dense(var, p, nmeas, f, fparams)
        integer, intent(in) :: var
        integer, intent(in) :: nmeas
        procedure(dl_callback_dense) :: f
        real(wp), intent(in), contiguous :: p(:)
        class(*), intent(in), target, optional :: fparams

        type(dl_dense_t), target :: cb
        integer(c_int) :: nstate_, nmeas_

        nstate_ = size(p, kind=c_int)
        nmeas_ = int(nmeas, kind=c_int)

        cb = dl_dense_t(nstate_, nmeas_, f, fparams)

        call dogleg_testGradient_dense(var, p, nstate_, nmeas_, &
            c_funloc(dl_callback_adaptor), c_loc(cb))

    end subroutine

    subroutine dl_freeContext(ctx)
        type(dl_solver_context_t), intent(inout) :: ctx
        call dogleg_freeContext(ctx%handle)
    end subroutine

end module