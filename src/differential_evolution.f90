module differential_evolution
    !------------------------------------------------------------------------------------------------------------------
    !! Implementation of the Differential Evolution optimization algorithm for finding global minumum of the
    !! multivariate function with boud constraints lower <= x <= upper
    !------------------------------------------------------------------------------------------------------------------
    use env, only: wp
    use interfaces, only: multivariate_fun, constraints, optimization_result
    use utils, only: rand_mat, penalty_function, randint


    type de_solver
        !! Differential Evolution solver class
        integer     :: population   !! Size of population
        real(wp)    :: cp           !! Crossover probability
        real(wp)    :: dw           !! Differential weight
    contains
        procedure :: optimize => optimize_de
    end type de_solver


    interface de_solver
        procedure :: de_solver_init
    end interface de_solver

contains

    function de_solver_init(population, cp, dw) result(self)
        !--------------------------------------------------------------------------------------------------------------
        !! Constructor for PSO solver
        !--------------------------------------------------------------------------------------------------------------
        type(de_solver) :: self
        integer, intent(in), optional   :: population   !! Number of particles
        real(wp), intent(in), optional  :: cp           !! Crossover probability
        real(wp), intent(in), optional  :: dw           !! Differential weight
        !--------------------------------------------------------------------------------------------------------------
        self%population = 80; if(present(population)) self%population = population
        self%cp = 0.9_wp; if(present(cp)) self%cp = cp
        self%dw = 0.8_wp; if(present(dw)) self%dw = dw
    end function de_solver_init


    subroutine optimize_de(self, f, lower, upper, iterations, solution)
        !--------------------------------------------------------------------------------------------------------------
        !! Solve optimization problem with Differential Evolution solver.
        !! Minimize f(x), with bounds: lower <= x <= upper
        !--------------------------------------------------------------------------------------------------------------
        class(de_solver)                        :: self
        procedure(multivariate_fun)             :: f            !! Objective function with `multivariate_fun` interface
        integer, intent(in)                     :: iterations   !! Number of iterations
        real(wp), intent(in)                    :: lower(:)     !! Vector of lower bounds for variable vector
        real(wp), intent(in)                    :: upper(:)     !! Vector of upper bounds for variable vector
        type(optimization_result), intent(out)  :: solution     !! Solution result
        !--------------------------------------------------------------------------------------------------------------
        real(wp), allocatable :: x(:, :), xp(:, :), v1(:), v2(:), v3(:), y(:), yp(:)
        logical, allocatable :: mask(:)
        integer :: ind(4)
        integer :: i, j, num, k, d
        real(kind=8) :: r

        ! Initial population
        m = self%population
        n = size(lower, dim=1)
        allocate(x(m, n), xp(m, n), v1(n), v2(n), v3(n), y(m), yp(m), mask(m))
        call rand_mat(lower, upper, x)

        ! Calculate objective function values
        do i = 1, m
            y(i) = f(x(i, :))
        end do

        ! DE iterations
        do i = 1, iterations
            mask(:) = .false.

            ! Loop over each vector
            do j = 1, m
                ind(1) = j

                ! Selection of the 3 different random vectors indices
                do num = 1, 3
                    do while (.true.)
                        k = randint(1, m)
                        if (any(ind(1: num) == k)) cycle
                        ind(num + 1) = k
                        exit
                    end do
                end do

                v1 = x(ind(2), :)
                v2 = x(ind(3), :)
                v3 = x(ind(4), :)
                d = randint(1, n)

                ! Mutant vector
                do k = 1, n
                    call random_number(r)
                    if (r < self%cp) then
                        xp(j, k) = v1(k) + self%dw * (v2(k) - v3(k))
                    else
                        xp(j, k) = x(j, k)
                    end if
                end do

                xp(j, d) = v1(d) + self%dw * (v2(d) - v3(d))
                where (xp(j, :) < lower) xp(j, :) = lower
                where (xp(j, :) > upper) xp(j, :) = upper
                yp(j) = f(xp(j, :))
                if (yp(j) <= y(j)) mask(j) = .true.
            end do

            ! New generation
            do j = 1, m
                if (mask(j)) then
                    x(j, :) = xp(j, :)
                    y(j) = yp(j)
                end if
            end do
        end do

        k = minloc(y, dim=1)
        solution%x = x(k, :)
        solution%fun = f(solution%x)
        solution%success = .true.
        solution%iterations = iterations
        solution%message = ""
    end subroutine optimize_de

end module differential_evolution