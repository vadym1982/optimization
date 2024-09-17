program test
    use env, only: wp
    use interfaces, only: optimization_result, no_constraints
    use particle_swarm, only: pso_solver
    use differential_evolution, only: de_solver
    use test_functions, only: rosenbrock, square, constraints_example, simonescu_fun, simonescu_constr, &
            holder_table_fun, holder_table_constr,  schaffer_function_n2

    implicit none

    type(pso_solver) :: solver1
    type(de_solver) :: solver2
    type(optimization_result) :: solution1, solution2
    integer :: i
    real(wp) :: t1, t2

    call cpu_time(t1)
    solver1 = pso_solver(p_count=80, weight=0.9_wp, c1=0.4_wp, c2=0.6_wp)

    call solver1%optimize(rosenbrock, no_constraints, [(-5.0_wp, i=1,7)], [(5.0_wp, i=1,7)], 1000, &
            solution1)
    call cpu_time(t2)

    print *, solution1%x, solution1%fun
    print *, solution1%message
    print "(A,F12.9)", "PSO. Optimization time: ", t2 - t1

    call cpu_time(t1)
    solver2 = de_solver(population=80)

    call solver2%optimize(rosenbrock, [(-5.0_wp, i=1,7)], [(5.0_wp, i=1,7)], 1000, &
            solution2)
    call cpu_time(t2)

    PRINT *, "ok"

    print *, solution2%x, solution2%fun
    print *, solution2%message
    print "(A,F12.9)", "DE. Optimization time: ", t2 - t1

end program test