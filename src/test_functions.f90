module test_functions
    !------------------------------------------------------------------------------------------------------------------
    !! Set of objective functions and constraints for solvers testing
    !------------------------------------------------------------------------------------------------------------------
    use env, only: wp, pi

contains

    function rosenbrock(x) result(f)
        !--------------------------------------------------------------------------------------------------------------
        !! Multivariate Rosenbrock function
        !--------------------------------------------------------------------------------------------------------------
        real(wp) :: x(:)
        real(wp) :: f
        !--------------------------------------------------------------------------------------------------------------
        integer :: n, i
        n = size(x, 1)
        f = 0.0_wp
        do i = 1, n - 1
            f = f + (1.0_wp - x(i)) ** 2 + 100.0_wp * (x(i+1) - x(i) ** 2) ** 2
        end do
    end function rosenbrock


    function square(x) result(f)
        !--------------------------------------------------------------------------------------------------------------
        !! Multivariate square function
        !--------------------------------------------------------------------------------------------------------------
        real(wp) :: x(:)
        real(wp) :: f
        !--------------------------------------------------------------------------------------------------------------
        f = sum(x ** 2)
    end function square


    function constraints_example(x) result(c)
        !--------------------------------------------------------------------------------------------------------------
        !! Multivariate square function
        !! Example of constraints function for constraints:
        !! x(1) >= 1; x(2) <= 3; x(2) - x(1) >= 0
        !--------------------------------------------------------------------------------------------------------------
        real(wp)                :: x(:)
        real(wp), allocatable   :: c(:)
        !--------------------------------------------------------------------------------------------------------------
        allocate(c(3))
        c(1) = 1.0_wp - x(1)
        c(2) = x(2) - 3.0_wp
        c(3) = x(1) - x(2)
    end function constraints_example


    function simonescu_fun(x) result(f)
        !--------------------------------------------------------------------------------------------------------------
        !! Simonescu function of 2 variables
        !--------------------------------------------------------------------------------------------------------------
        real(wp) :: x(:)
        real(wp) :: f
        !--------------------------------------------------------------------------------------------------------------
        f = 0.1_wp * x(1) * x(2)
    end function simonescu_fun


    function simonescu_constr(x) result(c)
        !--------------------------------------------------------------------------------------------------------------
        !! Constraints for Simonescu function
        !--------------------------------------------------------------------------------------------------------------
        real(wp)                :: x(:)
        real(wp), allocatable   :: c(:)
        !--------------------------------------------------------------------------------------------------------------
        allocate(c(1))
        c(1) = x(1) ** 2 + x(2) ** 2 - (1.0_wp + 0.2_wp * cos(8 * atan(x(1) / x(2)))) ** 2
    end function simonescu_constr


    function holder_table_fun(x) result(f)
        !--------------------------------------------------------------------------------------------------------------
        !! Hölder table function. f(+-8.00502, +-9.66459) = -19.285
        !--------------------------------------------------------------------------------------------------------------
        real(wp)    :: x(:)
        real(wp)    :: f
        !--------------------------------------------------------------------------------------------------------------
        f = -abs(sin(x(1)) * cos(x(2)) * exp(abs(1 - sqrt(x(1) ** 2 + x(2) ** 2) / pi)))
    end function holder_table_fun


    function holder_table_constr(x) result(c)
        !--------------------------------------------------------------------------------------------------------------
        !! Constraints for Hölder table function
        !--------------------------------------------------------------------------------------------------------------
        real(wp)                :: x(:)
        real(wp), allocatable   :: c(:)
        !--------------------------------------------------------------------------------------------------------------
        allocate(c(4))
        c(1) = x(1) - 10.0_wp
        c(2) = -10.0_wp - x(1)
        c(3) = x(2) - 10.0_wp
        c(4) = -10.0_wp - x(2)
    end function holder_table_constr


    function schaffer_function_n2(x) result(f)
        !--------------------------------------------------------------------------------------------------------------
        !! Hölder table function. f(+-8.00502, +-9.66459) = -19.285
        !--------------------------------------------------------------------------------------------------------------
        real(wp)    :: x(:)
        real(wp)    :: f
        !--------------------------------------------------------------------------------------------------------------
        f = 0.5_wp + (sin(x(1) ** 2 - x(2) ** 2) ** 2 - 0.5_wp) / &
                (1.0_wp + 0.001_wp * (x(1) ** 2 + x(2) ** 2)) ** 2
    end function schaffer_function_n2

end module test_functions