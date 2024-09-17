module interfaces
    use env, only: wp

    implicit none

    interface
        function multivariate_fun(x) result(f)
            ! Interface for the multivariate objective function
            import wp
            real(wp)    :: x(:)
            real(wp)    :: f
        end function multivariate_fun

        function constraints(x) result(c)
            !! Interface for the constrains functions. Result is vector of the values of the left hand sides for each
            !! constraint g_i(x) <= 0
            import wp
            real(wp)                :: x(:)
            real(wp), allocatable   :: c(:)
        end function constraints
    end interface


    type optimization_result
        real(wp), allocatable           :: x(:)         !! The resulting vector of optimized parameters
        real(wp)                        :: fun          !! The value of objective function at `x`
        logical                         :: success      !! Optimization success
        integer                         :: iterations   !! Number of iterations to reach a solution
        character(len=:), allocatable   :: message      !! Optimization process info
    end type optimization_result

contains

    function no_constraints(x) result(c)
        !--------------------------------------------------------------------------------------------------------------
        !! Function with `constraints` interface for empty constraints
        !--------------------------------------------------------------------------------------------------------------
        real(wp)                :: x(:)
        real(wp), allocatable   :: c(:)
        !--------------------------------------------------------------------------------------------------------------
        allocate(c(1))
        c = 0.0_wp
    end function no_constraints

end module interfaces