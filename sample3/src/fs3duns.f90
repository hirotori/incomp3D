program fs3duns_pois
    !/////////////////////////////////////////////////////////////////
    !! version: 1.0.0
    !! author: T.Ikeda
    !! summary:
    !! - メイン関数. 
    !! - ポアズイユ流のサンプル.
    !/////////////////////////////////////////////////////////////////
    use command_line_args_m
    use case_poiseuille_m
    use fractional_step_m
    use fractional_step_implicit_m
    use simulator_m
    type(args_t),allocatable :: args(:)
    type(case_poiseuille_t) :: current_case
    type(simulator_t) :: simulator
    class(solver_fs),allocatable :: solver

    args = argv()
    if ( size(args) /= 2 ) then
        print "('USAGE :: ')"
        print "('./fs3duns_pois.exe [solver type]')"
        print "('[solver type] = ""explicit"" or ""implicit""')"
        error stop
    end if

    if ( args(2)%v == "explicit" ) then
        allocate(solver_fs :: solver)
    else if ( args(2)%v == "implicit" ) then
        allocate(solver_fs_imp_t :: solver)
        select type (solver)
        type is (solver_fs_imp_t)
            call solver%set_parameter(1000, 1.0d-9, 1.0d0)
        end select
    else 
        error stop "invalid argument :: "//args(2)%v
    end if

    call simulator%run(current_case, solver)
   
end program fs3duns_pois