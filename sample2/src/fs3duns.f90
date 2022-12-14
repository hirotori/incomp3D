program fs3duns
    !/////////////////////////////////////////////////////////////////
    !! version: 1.0.0
    !! author: T.Ikeda
    !! summary:
    !! - メイン関数. 
    !! - あらかじめ用意されたケースクラス, シミュレータ, アルゴリズムにより
    !!   解析を行うサンプルプログラム. 
    !/////////////////////////////////////////////////////////////////
    use case_cavity_m
    use fractional_step_implicit_m
    use simulator_m

    type(case_cavity_t) :: current_case
    type(simulator_t) :: simulator
    type(solver_fs_imp_t) :: solver

    call solver%set_parameter(1000, 1.0d-4, 1.0d0)
    call simulator%run(current_case, solver)
   
end program fs3duns