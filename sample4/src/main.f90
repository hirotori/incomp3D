program main
    use fractional_step_m
    use simulator_m
    use case_common_m
    implicit none
    type(case_common_t) :: this_case
    type(simulator_t) :: sim
    type(solver_fs) :: solver

    call sim%run(this_case, solver)
    
end program main