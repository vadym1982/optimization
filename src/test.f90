program test
    use env, only: wp
    use interfaces, only: optimization_result, no_constraints
    use particle_swarm, only: pso_solver
    use test_functions, only: rosenbrock, square, constraints_example, simonescu_fun, simonescu_constr, &
            holder_table_fun, holder_table_constr

    implicit none

    type(pso_solver) :: solver1
    type(optimization_result) :: solution1
    integer :: i

    solver1 = pso_solver(p_count=80, weight=0.9_wp, c1=0.4_wp, c2=0.6_wp)
    call solver1%optimize(holder_table_fun, holder_table_constr, [(-10d0, i=1,2)], [(10d0, i=1,2)], 1000, solution1)

    print *, solution1%x, solution1%fun
    print *, solution1%message


end program test