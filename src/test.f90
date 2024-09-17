program test
    use env, only: wp
    use interfaces, only: optimization_result, no_constraints
    use particle_swarm, only: pso_solver
    use test_functions, only: rosenbrock, square, constraints_example, simonescu_fun, simonescu_constr, &
            holder_table_fun, holder_table_constr,  schaffer_function_n2

    implicit none

    type(pso_solver) :: solver1
    type(optimization_result) :: solution1
    integer :: i
    real(wp) :: t1, t2

    call cpu_time(t1)
    solver1 = pso_solver(p_count=80, weight=0.9_wp, c1=0.4_wp, c2=0.6_wp)
    call solver1%optimize(schaffer_function_n2, no_constraints, [(-50.0_wp, i=1,2)], [(100.0_wp, i=1,2)], 1000, &
            solution1)
    call cpu_time(t2)

    print *, solution1%x, solution1%fun
    print *, solution1%message
    print "(A,F12.9)", "Optimization time: ", t2 - t1

end program test