!
! The size of C unsigned int is supposed to be equal to int
! see https://docs.microsoft.com/en-us/cpp/cpp/data-type-ranges?view=msvc-170
!
!

module dogleg_c

    use, intrinsic :: iso_c_binding
    implicit none

    public

    type, bind(c) :: dogleg_parameters2_t
        integer(c_int) :: max_iterations
        integer(c_int) :: dogleg_debug

        real(c_double) :: trustregion0

        ! These are probably OK to leave alone. Tweaking them can maybe result in
        ! slighlty faster convergence
        real(c_double) :: trustregion_decrease_factor
        real(c_double) :: trustregion_decrease_threshold
        real(c_double) :: trustregion_increase_factor
        real(c_double) :: trustregion_increase_threshold

        ! The termination thresholds. Documented in the header
        real(c_double) :: Jt_x_threshold
        real(c_double) :: update_threshold
        real(c_double) :: trustregion_threshold
    end type

    interface

        subroutine dogleg_getDefaultParameters(parameters) &
            bind(c, name = "dogleg_getDefaultParameters")
            import c_ptr
            type(c_ptr), value :: parameters
        end subroutine


        subroutine dogleg_setMaxIterations(n) &
            bind(c, name="dogleg_setMaxIterations")
            import c_int
            integer(c_int), value :: n
        end subroutine

        subroutine dogleg_setTrustregionUpdateParameters( &
            downFactor, downThreshold, upFactor, upThreshold) &
            bind(c, name="dogleg_setTrustregionUpdateParameters")
            import c_double
            real(c_double), value :: downFactor, downThreshold
            real(c_double), value :: upFactor, upThreshold
        end subroutine

        subroutine dogleg_setDebug(debug) bind(c,name="dogleg_setDebug")
            import c_int
            integer(c_int), value :: debug
        end subroutine

        subroutine dogleg_setInitialTrustregion(t) &
            bind(c, name="dogleg_setInitialTrustregion")
            import c_double
            real(c_double), value :: t
        end subroutine

        subroutine dogleg_setThresholds(Jt_x, update, trustregion) &
            bind(c, name='dogleg_setThresholds')
            import c_double
            real(c_double), value :: Jt_x, update, trustregion
        end subroutine

        function dogleg_optimize_dense(p, Nstate, Nmeas, f, cookie, returnContext) &
            bind(c, name='dogleg_optimize_dense')
            import c_int, c_double, c_funptr, c_ptr
            real(c_double) :: dogleg_optimize_dense
                !> Return value

            real(c_double), intent(inout) :: p(*)
                !> Initial estimate of the state vector (and holds the final result)
            integer(c_int), value, intent(in) :: Nstate
                !> Number of optimization variables (elements of p)
            integer(c_int), value, intent(in) :: Nmeas
                !> Number of measurements (elements of x). `Nmeas >= Nstate` is a requirement
            type(c_funptr), value, intent(in) :: f
                !> The callback function that the optimization routine calls
            type(c_ptr), value, intent(in) :: cookie
                !> Arbitrary data pointer passed untouched to `f`
            type(c_ptr) :: returnContext
                !> Used to retrieve the full context structure from the solver.
        end function

        function dogleg_optimize_dense2( &
            p, Nstate, Nmeas, f, cookie, parameters, returnContext) &
            bind(c, name='dogleg_optimize_dense2')
            import c_int, c_double, c_funptr, c_ptr
            real(c_double) :: dogleg_optimize_dense2
                !> Return value

            real(c_double), intent(inout) :: p(*)
                !> Initial estimate of the state vector (and holds the final result)
            integer(c_int), value :: Nstate
                !> Number of optimization variables (elements of p)
            integer(c_int), value :: Nmeas
                !> Number of measurements (elements of x). `Nmeas >= Nstate` is a requirement
            type(c_funptr), value, intent(in) :: f
                !> The callback function that the optimization routine calls
            type(c_ptr), value :: cookie
                !> Arbitrary data pointer passed untouched to `f`
            type(c_ptr), value, intent(in) :: parameters
                !> Pointer to the set of parameters to use. 
            type(c_ptr) :: returnContext
                !> Used to retrieve the full context structure from the solver.
        end function

!       subroutine dogleg_computeJtJfactorization(point, ctx) &
!           bind(c, name="dogleg_computeJtJfactorization")
!           import c_ptr
!           type(c_ptr) :: point
!           type(c_ptr) :: ctx
!       end subroutine

        subroutine dogleg_testGradient_dense(var, p0, Nstate, Nmeas, f, cookie) &
            bind(c, name="dogleg_testGradient_dense")
            import c_int, c_double, c_funptr, c_ptr
            integer(c_int), value :: var
            real(c_double), intent(in) :: p0(*)
            integer(c_int), value :: Nstate
            integer(c_int), value :: Nmeas
            type(c_funptr) :: f
            type(c_ptr), value :: cookie
        end subroutine

        subroutine dogleg_freeContext(ctx) bind(c,name="dogleg_freeContext")
            import c_ptr
            type(c_ptr) :: ctx
        end subroutine

    end interface

end module
