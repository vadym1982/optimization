module utils
    use env, only: wp
    use interfaces, only: multivariate_fun

contains

    subroutine rand_mat(lower, upper, mat)
        !--------------------------------------------------------------------------------------------------------------
        !! Generate random matrix of values with lower and upper bounds lower(j) <= mat(i, j) <= upper(j)
        !--------------------------------------------------------------------------------------------------------------
        real(kind=8), intent(in)    :: lower(:), upper(:)   !! Arrays of lower and upper bounds for matrix elements
        real(kind=8), intent(out)   :: mat(:, :)            !! Random matrix
        !--------------------------------------------------------------------------------------------------------------
        integer :: n
        call random_number(mat)
        n = size(mat, 1)
        mat = spread(lower, 1, n) + spread(upper - lower, 1, n) * mat
    end subroutine rand_mat


!    subroutine calc_fun(f, c, x, y, mask)
!        !--------------------------------------------------------------------------------------------------------------
!        !! Calculate multivariate function with `fun` interface for each row in matrix `x`
!        !--------------------------------------------------------------------------------------------------------------
!        procedure(multivariate_fun) :: f        !! Multivariate function with `multivariate_fun` interface
!        procedure(constraints)      :: c        !! Constraints function with `constraints` interface
!        real(kind=8)                :: x(:, :)  !! Matrix of arguments. Each row represents a vector of function
!                                                !! arguments values
!        real(kind=8)                :: y(:)     !! Array of function values y(i) = f(x(i, :))
!        logical                     :: mask(:)  !! An array of results for satisfying constraints for each row of `x`
!        !--------------------------------------------------------------------------------------------------------------
!        integer :: i
!
!        do i = 1, size(x, 1)
!            y(i) = f(x(i, :))
!            mask(i) = all(c(x(i, :)))
!        end do
!    end subroutine calc_fun

end module utils