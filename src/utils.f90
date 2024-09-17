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


end module utils