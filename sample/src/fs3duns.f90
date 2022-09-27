program fs3duns
    !/////////////////////////////////////////////////////////////////
    !! version: 1.0.0
    !! author: T.Ikeda
    !! summary:
    !! - メイン関数. 
    !/////////////////////////////////////////////////////////////////
    use case_tgv_m
    use fractional_step_implicit_m
    use simulator_m

    type(case_tgv_t) :: current_case
    type(simulator_t) :: simulator
    type(solver_fs_imp_t) :: solver
    call simulator%run(current_case, solver)
   
end program fs3duns