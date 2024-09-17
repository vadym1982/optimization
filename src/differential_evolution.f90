module differential_evolution
    !------------------------------------------------------------------------------------------------------------------
    !! Implementation of the Differential Evolution optimization algorithm
    !------------------------------------------------------------------------------------------------------------------
    use env, only: wp
    use interfaces, only: multivariate_fun, constraints, optimization_result
    use utils, only: rand_mat, penalty_function


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
        type(pso_solver) :: self
        integer, intent(in), optional   :: population   !! Number of particles
        real(wp), intent(in), optional  :: cp           !! Crossover probability
        real(wp), intent(in), optional  :: dw           !! Differential weight
        !--------------------------------------------------------------------------------------------------------------
        self%population = 80; if(present(population)) self%population = population
        self%cp = 0.9_wp; if(present(cp)) self%cp = cp
        self%dw = 0.8_wp; if(present(dw)) self%dw = dw
    end function de_solver_init

end module differential_evolution