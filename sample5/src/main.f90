program main
    use case_channel_m
    ! use les_solver_m
    use fractional_step_implicit_m
    use simulator_m
    implicit none
    type(case_channel_t) case_
    type(solver_fs_imp_t) solver_
    type(simulator_t) simulator_

    call solver_%set_parameter(1000, 1.0d-6, 1.3d0)
    call simulator_%run(case_, solver_)
    
end program main