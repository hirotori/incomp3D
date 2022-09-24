program fs3duns
    !/////////////////////////////////////////////////////////////////
    !! version: 1.0.0
    !! author: T.Ikeda
    !! summary:
    !! - メイン関数. 
    !/////////////////////////////////////////////////////////////////
    use case_common_m
    use fractional_step_m
    use simulator_m

    type(case_common_t) :: current_case
    type(simulator_t) :: simulator
    type(solver_fs) :: solver
    call simulator%run(current_case, solver)
   
end program fs3duns