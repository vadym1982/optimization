module conjugate_gradient
    !------------------------------------------------------------------------------------------------------------------
    !! Implementation of the Conjugate Gradient (Fletcher–Reeves) optimization algorithm with multistart
    !! for finding global minumum of the multivariate function with boud constraints lower <= x <= upper
    !------------------------------------------------------------------------------------------------------------------
    use env, only: wp
    use interfaces, only: multivariate_fun, constraints, optimization_result
    use utils, only: rand_mat, grad

    implicit none

    real(kind=8), parameter :: inc_factor = 1.3_wp, dec_factor = 0.5_wp


    type cg_solver
        !! Conjugate Gradient solver class
        integer     :: attempts     !! Number of attempts to find minimum from different seeds
        real(wp)    :: delta        !! Argument delta for gradient calculation
        integer     :: iter_limit   !! Max number of iterations
    contains
        procedure :: optimize => optimize_cg
    end type cg_solver


    interface cg_solver
        procedure :: cg_solver_init
    end interface cg_solver

contains

    function cg_solver_init(attempts, delta, iter_limit) result(self)
        !--------------------------------------------------------------------------------------------------------------
        !! Constructor for PSO solver
        !--------------------------------------------------------------------------------------------------------------
        type(cg_solver) :: self
        integer, intent(in), optional   :: attempts     !! Number of attempts
        real(wp), intent(in), optional  :: delta        !! Argument delta for gradient calculation
        integer, intent(in), optional   :: iter_limit   !! Max number of iterations
        !--------------------------------------------------------------------------------------------------------------
        self%attempts = 80; if(present(attempts)) self%attempts = attempts
        self%delta = 1.0e-7_wp; if(present(delta)) self%delta = delta
        self%iter_limit = 500; if(present(iter_limit)) self%iter_limit = iter_limit
    end function cg_solver_init


    subroutine opt_cg(f, x0, n, tol, max_iter, x_opt, y_opt, iter, success)
        !--------------------------------------------------------------------------------------------------------------
        !! Nonlinear unconstraibed optimization with Fletcher–Reeves conjugate gradient method
        !--------------------------------------------------------------------------------------------------------------
        procedure(multivariate_fun) :: f            !! Objective function with 'fun' interface
        real(kind=8), intent(in)    :: x0(:)        !! Seed vector
        integer, intent(in)         :: n            !! Problen dimension
        real(kind=8), intent(in)    :: tol          !! Absolute convergence tolerance
        integer, intent(in)         :: max_iter     !! Maximum number of iteration for sesrching
                                                    !! When this number of iterations is reached,
                                                    !! the algorithm is stopped and the success
                                                    !! argument is set to false.
        real(kind=8), intent(out)   :: x_opt(:)     !! Vector if the optimal values
        real(kind=8), intent(out)   :: y_opt        !! Value of the objective function at the optimum
        integer, intent(out)        :: iter         !! Number of iterations to achive algorithm convergence
        logical, intent(out)        :: success   !  ! Optimuzation success
        !--------------------------------------------------------------------------------------------------------------
        integer :: i
        real(kind=8), allocatable :: x(:), g1(:), g2(:), s1(:), s2(:)
        real(kind=8) :: w, a = 1.0d0, y, g1dot, g2dot

        allocate(x(n), g1(n), g2(n), s1(n), s2(n))
        iter = 0
        x = x0
        y_opt = f(x0)
        success = .true.

        ! Main algorithm loop
        do
            call grad(f, x, 1.0d-7, g2)
            g2dot = dot_product(g2, g2)

            if (norm2(g2) <= epsilon(1.0d0)) then
                x_opt = x
                y_opt = f(x)
                return
            end if

            if (mod(iter, n) == 0) then
                ! Refreshing search direction
                s2 = -g2
            else
                ! Updating search direction
                w = (g2dot / g1dot)
                s2 = -g2 + w * s1
                if (dot_product(s2, g2) >= 0d0) then
                    ! Refreshing if the function does not decrease in the direction
                    s2 = -g2
                end if
            end if

            ! Search along s2
            outer: do
                y = f(x + a * s2)
                if (y < y_opt) then
                    ! Increasing step
                    i = 0
                    do
                        y_opt = y
                        a = a * inc_factor
                        y = f(x + a * s2)
                        if (y >= y_opt) then
                            a = a / inc_factor
                            exit outer
                        end if
                        i = i + 1
                        if (i > max_iter) then
                            y_opt = y
                            exit outer
                        end if
                    end do
                else
                    i = 0
                    do
                        ! Decreasing step
                        a = a * dec_factor
                        y = f(x + a * s2)
                        if (y < y_opt) then
                            y_opt = y
                            exit outer
                        end if
                        i = i + 1
                        if (i > max_iter * 10) then
                            exit outer
                        end if
                    end do
                end if
            end do outer

            ! Accepting improved solution
            x_opt = x + a * s2

            ! Convergence checking
            if (norm2(x - x_opt) <= tol) then
                exit
            end if

            x = x_opt
            g1 = g2
            s1 = s2
            g1dot = g2dot
            iter = iter + 1

            ! Iterations controll
            if (iter >= max_iter) then
                success = .false.
                return
            end if
        end do

    end subroutine opt_cg


    subroutine optimize_cg(self, f, lower, upper, tol, solution)
        !--------------------------------------------------------------------------------------------------------------
        !! Solve optimization problem with Differential Evolution solver.
        !! Minimize f(x), with bounds: lower <= x <= upper
        !--------------------------------------------------------------------------------------------------------------
        class(cg_solver)                        :: self
        procedure(multivariate_fun)             :: f            !! Objective function with `multivariate_fun` interface
        real(wp), intent(in)                    :: tol          !! Absolute convergence tolerance
        real(wp), intent(in)                    :: lower(:)     !! Vector of lower bounds for variable vector
        real(wp), intent(in)                    :: upper(:)     !! Vector of upper bounds for variable vector
        type(optimization_result), intent(out)  :: solution     !! Solution result
        !--------------------------------------------------------------------------------------------------------------
        real(wp), allocatable :: x_opt(:), best_x(:), x0(:, :)
        real(wp) :: y_opt, best_y
        integer :: i, n, iter, best_iter
        logical :: success, is_found

        n = size(lower, dim=1)
        is_found = .false.
        best_y = huge(best_y)
        allocate(x0(self%attempts, n), x_opt(n), best_x(n))

        do i = 1, self%attempts
            call rand_mat(lower, upper, x0)
            call opt_cg(f, x0(i, :), n, tol, self%iter_limit, x_opt, y_opt, iter, success)

            if (.not. success) cycle
            if (any(x_opt < lower) .or. any(x_opt > upper)) cycle

            if (y_opt < best_y) then
                best_x = x_opt
                best_y = y_opt
                is_found = .true.
                best_iter = iter
            end if
        end do

        if (is_found) then
            solution%x = best_x
            solution%fun = best_y
            solution%success = .true.
            solution%iterations = best_iter
            solution%message = ""
        else
            solution%success = .false.
            solution%message = "Unable to find solution"
        end if
    end subroutine optimize_cg


end module conjugate_gradient