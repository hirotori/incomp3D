module simulator_m
    use floating_point_parameter_m, only : ip, dp
    use case_common_m
    use boundary_condition_m
    use mesh_m
    use fluid_field_m
    use fractional_step_m
    use time_monitor_m

    implicit none
    type, public :: simulator_t
        real(dp),allocatable :: v_old(:,:,:,:)

        contains
        procedure :: run
    
    end type 
contains

subroutine run(this, current_case, solver)
    !!各ケースに合わせた非圧縮流れ計算を実施する.
    use,intrinsic :: ieee_arithmetic, only : ieee_is_nan
    class(simulator_t),intent(inout) :: this
    class(case_common_t),intent(inout) :: current_case
        !!計算ケース. 
    class(solver_fs),intent(inout) :: solver
        !!非圧縮流れ計算のアルゴリズム.

    type(rectilinear_mesh_t) :: grid
        !!計算空間.
    type(fluid_field_t) fluid
        !!流れ場.

    integer(ip) nstep, extents(3)
    real(dp) dt, re
    character(*),parameter :: time_bar = repeat("=", 35)
    type(date_t) :: start_date, end_date
    logical sim_diverged

    start_date = get_current_date_and_time()

    !ループ前の処理
    call current_case%phase_pre_process(grid, fluid)
    extents(:) = grid%get_extents()
    call boundary_condition_velocity(extents, fluid%velocity, current_case%bc_types)
    call boundary_condition_pressure(extents, fluid%pressure, current_case%bc_types)

    !陰解法のため, ループ前の初期値で面流束を求める. 初期値は連続の式を満たすと仮定して, 
    !線形補間だけ行う. 厳密な圧力場が分かっていればその直後に補正をかけても良いが, だいたい一様初期値なので意味が無いとして保留.
    call cal_face_velocity(extents(1), extents(2), extents(3), fluid%velocity, fluid%mflux_i, fluid%mflux_j, fluid%mflux_k)

    allocate(this%v_old, source = fluid%velocity)
    call solver%init(fluid, grid, current_case%settings_solver, current_case%settings_case) 

    call current_case%phase_writeout(grid, fluid, 0)

    dt = current_case%settings_case%dt
    re = current_case%settings_case%reynolds_number
    do nstep = current_case%settings_case%nstart, current_case%settings_case%nend
        print "(A)", time_bar
        print "('nstep = ',i0, ', time = ',g0)", nstep, real(nstep,dp)*current_case%settings_case%dt
        print "(A)", time_bar

        !どのアルゴリズムでも共通して計算させる.
        call calc_gradient_tensor(extents, grid%dv, grid%ds, grid%dx, fluid%velocity, fluid%dudr)

        call solver%predict_pseudo_velocity(extents, grid%ds, grid%dv, grid%dx, dt, re, &
                                            current_case%settings_case%body_force, &
                                            this%v_old, &
                                            fluid%velocity, fluid%mflux_i, fluid%mflux_j, fluid%mflux_k, &
                                            current_case%bc_types)

        call boundary_condition_velocity(extents, fluid%velocity, current_case%bc_types)

        call solver%calc_corrected_velocity(extents, grid%ds, grid%dv, dt, &
                                           fluid%mflux_i, fluid%mflux_j, fluid%mflux_k, fluid%pressure, fluid%velocity, &
                                           current_case%settings_solver, current_case%settings_case%p_ref, current_case%bc_types, &
                                           sim_diverged)

        call boundary_condition_velocity(extents, fluid%velocity, current_case%bc_types)

        this%v_old(:,:,:,:) = fluid%velocity(:,:,:,:)

        if ( sim_diverged ) then
            print "('ERROR :: Calculation Terminated.')"
            call current_case%process_diverged(grid, fluid, nstep)
            return
        end if

        call current_case%phase_post_process(grid, fluid)

        call current_case%check_flow_field(grid, fluid)

        if ( mod(nstep, current_case%settings_case%nwrite) == 0) then
            call current_case%phase_writeout(grid, fluid, nstep)
        end if

    end do

    end_date = get_current_date_and_time()
    call start_date%print_self()
    call end_date%print_self()

end subroutine
    
end module simulator_m