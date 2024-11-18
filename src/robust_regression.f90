module robust_regeression
    use env, only: wp
    use differential_evolution

    implicit none

    interface
        function regression_fun(x, c) result(y)
            import wp
            real(wp), intent(in)    :: x(:, :)
            real(wp), intent(in)    :: c(:)
            real(wp), allocatable   :: y(:, :)
        end function regression_fun

        function loss_fun(res, alpha) result(loss)
            import wp
            real(wp), intent(in)    :: res(:, :)
            real(wp), intent(in)    :: alpha
            real(wp), allocatable   :: loss(:)
        end function loss_fun
    end interface

contains

    function squered_loss(res, alpha) result(loss)
        !--------------------------------------------------------------------------------------------------------------
        !! Squared loss
        !--------------------------------------------------------------------------------------------------------------
        real(wp), intent(in)    :: res(:, :)
        real(wp), intent(in)    :: alpha
        real(wp), allocatable   :: loss(:)
        !--------------------------------------------------------------------------------------------------------------
        loss = sum(res ** 2, dim=2)
    end function squered_loss


    function euclidean_loss(res, alpha) result(loss)
        !--------------------------------------------------------------------------------------------------------------
        !! Absolute loss
        !--------------------------------------------------------------------------------------------------------------
        real(wp), intent(in)    :: res(:, :)
        real(wp), intent(in)    :: alpha
        real(wp), allocatable   :: loss(:)
        !--------------------------------------------------------------------------------------------------------------
        loss = sqrt(sum(res ** 2, dim=2))
    end function euclidean_loss


    function atan_loss(res, alpha) result(loss)
        !--------------------------------------------------------------------------------------------------------------
        !! Atan loss
        !--------------------------------------------------------------------------------------------------------------
        real(wp), intent(in)    :: res(:, :)
        real(wp), intent(in)    :: alpha
        real(wp), allocatable   :: loss(:)
        !--------------------------------------------------------------------------------------------------------------
        loss = atan(sqrt(sum(res ** 2, dim=2)) / alpha)
    end function atan_loss


    subroutine regression(fun, loss, x, y, alpha, lower, upper, population, iterations, c)
        !--------------------------------------------------------------------------------------------------------------
        !! Solve regression problem
        !--------------------------------------------------------------------------------------------------------------
        procedure (regression_fun)      :: fun
        procedure (loss_fun)            :: loss
        real(wp), intent(in)            :: x(:, :)
        real(wp), intent(in)            :: y(:, :)
        real(wp), intent(in)            :: alpha
        real(wp), intent(in)            :: lower(:)     !! Vector of lower bounds for variable vector
        real(wp), intent(in)            :: upper(:)     !! Vector of upper bounds for variable vector
        integer, optional               :: population
        integer, intent(in), optional   :: iterations   !! Number of iterations
        real(wp), intent(out)           :: c(:)
        !--------------------------------------------------------------------------------------------------------------
        type (de_solver)            :: solver
        type (optimization_result)  :: solution
        integer :: iter, p

        iter = 1000; if (present(iterations)) iter = iterations
        p = 80; if (present(population)) p = population

        solver = de_solver(population=p)
        call solver%optimize(obj, lower, upper, iter, solution)
        c = solution%x

    contains

        function obj(c) result(f)
            real(wp)    :: c(:)
            real(wp)    :: f
            f = sum(loss(fun(x, c) - y, alpha))
        end function obj

    end subroutine regression

end module robust_regeression