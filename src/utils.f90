module utils
    use env, only: wp
    use interfaces, only: multivariate_fun, constraints

contains

    function randint(m, n) result(r)
        !--------------------------------------------------------------------------------------------------------------
        !! Get random integer from `m` to `n` inclusive
        !--------------------------------------------------------------------------------------------------------------
        integer :: m, n
        integer :: r
        !--------------------------------------------------------------------------------------------------------------
        real(wp) :: u
        call random_number(u)

        r = floor(m + (n - m + 1) * u)
    end function randint

    subroutine rand_mat(lower, upper, mat)
        !--------------------------------------------------------------------------------------------------------------
        !! Generate random matrix of values with lower and upper bounds lower(j) <= mat(i, j) <= upper(j).
        !--------------------------------------------------------------------------------------------------------------
        real(kind=8), intent(in)    :: lower(:), upper(:)   !! Arrays of lower and upper bounds for matrix elements
        real(kind=8), intent(out)   :: mat(:, :)            !! Random matrix
        !--------------------------------------------------------------------------------------------------------------
        integer :: n
        call random_number(mat)
        n = size(mat, 1)
        mat = spread(lower, 1, n) + spread(upper - lower, 1, n) * mat
    end subroutine rand_mat


    function penalty_function(fun, constr, x, penalty) result(f)
        !--------------------------------------------------------------------------------------------------------------
        !! Calculate objective function with penalty for violating constraints
        !--------------------------------------------------------------------------------------------------------------
        procedure(multivariate_fun) :: fun
        procedure(constraints)      :: constr
        real(wp)                    :: x(:)
        real(wp)                    :: penalty
        !--------------------------------------------------------------------------------------------------------------
        f = fun(x) + penalty * sum(max(0.0_wp, constr(x)) ** 2)
    end function penalty_function


    subroutine grad(f, x, h, g)
        !--------------------------------------------------------------------------------------------------------------
        !! Numerical estimation of the gradient of a multivariate function using finite diference scheme
        !--------------------------------------------------------------------------------------------------------------
        procedure(multivariate_fun) :: f     !! Multivariate function with 'fun' interface
        real(kind=8), intent(in)    :: x(:)  !! Point for calculation of the gradient
        real(kind=8), intent(in)    :: h     !! Argument increment << 1
        real(kind=8), intent(out)   :: g(:)  !! Regulting gradient vector
        !--------------------------------------------------------------------------------------------------------------
        integer :: i, n
        real(kind=8), allocatable :: x1(:), x2(:)

        n = size(x)
        allocate(x1(n), x2(n))
        x1 = x
        x2 = x

        do i = 1, n
            if (i > 1) then
                x1(i - 1) = x(i - 1)
                x2(i - 1) = x(i - 1)
            end if
            x1(i) = x(i) - h
            x2(i) = x(i) + h
            g(i) = (f(x2) - f(x1)) / h / 2.0d0
        end do
    end subroutine grad


end module utils