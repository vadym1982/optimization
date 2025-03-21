module particle_swarm
    !------------------------------------------------------------------------------------------------------------------
    !! Implementation of the Particle Swarm Optimization algorithm
    !------------------------------------------------------------------------------------------------------------------
    use env, only: wp
    use interfaces, only: multivariate_fun, constraints, optimization_result
    use utils, only: rand_mat

    implicit none


    type pso_solver
        !! PSO solver class
        integer     :: p_count     !! Number of partickles
        real(wp)    :: weight      !! Old speed veights factor
        real(wp)    :: c1          !! Local trend weight factor
        real(wp)    :: c2          !! Global trend weight factor
    contains
        procedure :: optimize => optimize_pso
    end type pso_solver


    interface pso_solver
        procedure :: pso_solver_init
    end interface pso_solver

contains

    function pso_solver_init(p_count, weight, c1, c2) result(self)
        !--------------------------------------------------------------------------------------------------------------
        !! Constructor for PSO solver
        !--------------------------------------------------------------------------------------------------------------
        type(pso_solver) :: self
        integer, intent(in), optional   :: p_count  !! Number of particles
        real(wp), intent(in), optional  :: weight   !! Momentum factor
        real(wp), intent(in), optional  :: c1       !! Local trend factor
        real(wp), intent(in), optional  :: c2       !! Global trend factor
        !--------------------------------------------------------------------------------------------------------------
        self%p_count = 80; if(present(p_count)) self%p_count = p_count
        self%weight = 0.9_wp; if(present(weight)) self%weight = weight
        self%c1 = 0.5_wp; if(present(c1)) self%c1 = c1
        self%c2 = 0.3_wp; if(present(c2)) self%c1 = c2
    end function pso_solver_init


    subroutine optimize_pso(self, f, constr, lower, upper, iterations, solution)
        !--------------------------------------------------------------------------------------------------------------
        !! Solve optimization problem with PSO solver.
        !! Minimize f(x), subjected to constraints:
        !!    constr(x) <= 0
        !!    lower <= x <= upper
        !--------------------------------------------------------------------------------------------------------------
        class(pso_solver)                       :: self
        procedure(multivariate_fun)             :: f            !! Objective function with `multivariate_fun` interface
        procedure(constraints)                  :: constr       !! Constraints function with `constraints` interface
        integer, intent(in)                     :: iterations   !! Number of iterations
        real(wp), intent(in)                    :: lower(:)     !! Vector of lower bounds for variable vector
        real(wp), intent(in)                    :: upper(:)     !! Vector of upper bounds for variable vector
        type(optimization_result), intent(out)  :: solution     !! Solution result
        !--------------------------------------------------------------------------------------------------------------
        integer :: m, n, i, i_min, j, tmp
        real(wp), allocatable :: x(:, :), v(:, :), p(:, :), g(:)
        real(wp), allocatable :: yp(:), delta(:), r(:)
        real(wp) :: y, y_min
        logical, allocatable :: mask_p(:)
        logical :: valid_solution, valid

        m = self%p_count
        n = size(lower, dim=1)
        allocate(x(m, n), v(m, n), p(m, n), g(n))
        allocate(yp(m), r(n))
        allocate(mask_p(m))

        call rand_mat(lower, upper, x)
        delta = (upper - lower)
        call rand_mat(-delta, delta, v)
        valid_solution = .false.
        y_min = huge(y_min)
        i_min = 0

        ! Initial population
        do i = 1, m
            mask_p(i) = all(constr(x(i, :)) <= 0.0_wp) .and. all(x(i, :) >= lower) .and. all(x(i, :) <= upper)

            if (mask_p(i)) then
                valid_solution = .true.
                yp(i) = f(x(i, :))
                p(i, :) = x(i, :)

                if (yp(i) < y_min) then
                    y_min = yp(i)
                    i_min = i
                end if
            end if
        end do

        if (valid_solution) g(:) = x(i_min, :)

        ! PSO iterations
        do i = 1, iterations
            ! Momentum speed
            v = self%weight * v

            do j = 1, m
                ! Local trend
                if (mask_p(j)) then
                    call random_number(r)
                    v(j, :) = v(j, :) + self%c1 * r * (p(j, :) - x(j, :))
                end if

                ! Global trend
                if (valid_solution) then
                    call random_number(r)
                    v(j, :) = v(j, :) + self%c2 * r * (g - x(j, :))
                end if
            end do

            x = x + v
            i_min = 0

            ! Update trends
            do j = 1, m
                valid = all(constr(x(j, :)) <= 0.0_wp) .and. all(x(j, :) >= lower) .and. all(x(j, :) <= upper)

                if (valid) then
                    valid_solution = .true.
                    y = f(x(j, :))

                    if ((y < yp(j)) .or. (.not. mask_p(j))) then
                        yp(j) = y
                        p(j, :) = x(j, :)
                        mask_p(j) = .true.
                    end if

                    if (y < y_min) then
                        i_min = j
                        y_min = y
                    end if
                end if
            end do
            if (i_min > 0) g(:) = x(i_min, :)
        end do

        ! Optimization result
        if (valid_solution) then
            solution%x = g
            solution%fun = f(solution%x)
            solution%success = .true.
            solution%iterations = iterations
            solution%message = ""
        else
            solution%success = .false.
            solution%iterations = iterations
            solution%message = "Could not find a solution that satisfies the constraints"
        end if

    end subroutine optimize_pso

end module particle_swarm