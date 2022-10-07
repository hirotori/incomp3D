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
    !$ use omp_lib
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
    integer(ip) :: nthread = 1
        !!コンソールに表示するための並列スレッド数. 並列化されない場合はシングルスレッドを表す1.
    real(dp) dt, re
    character(*),parameter :: time_bar = repeat("=", 35)
    type(date_t) :: start_date, end_date
    logical sim_diverged

    !並列スレッド数の指定.
    !$ call omp_set_num_threads(omp_get_max_threads())
    !$ nthread = omp_get_max_threads()

    start_date = get_current_date_and_time()

    !!ループ前の処理. ケースクラスで場を初期化かつメッシュの生成. 
    call current_case%phase_pre_process(grid, fluid)
    extents(:) = grid%get_extents()
    !!その後ソルバの共通パラメータ(空間スキームの選択)を設定する.
    call solver%init(fluid, grid, current_case%settings_solver, current_case%settings_case) 
    !!次に境界条件を適用させ, 仮想セルを更新する. この処理はシミュレーターでなくソルバに任せる.
    !!@note ここはシミュレーターに任せても問題は無いが, 面フラックスの修正有り/無しがソルバに委ねられているため.
    call solver%process_before_loop(grid, fluid, current_case%settings_case, current_case%bc_types)
    
    !既知の速度場を保持しておく.
    allocate(this%v_old, source = fluid%velocity)
    
    !初期場を書き出す.
    call current_case%phase_writeout(grid, fluid, 0)

    dt = current_case%settings_case%dt
    re = current_case%settings_case%reynolds_number
    do nstep = current_case%settings_case%nstart, current_case%settings_case%nend
        !!メインループでは,  1.中間速度の計算, 2.圧力補正の順に行う.
        print "(A)", time_bar
        print "('nstep = ',i0, ', time = ',g0)", nstep, real(nstep,dp)*current_case%settings_case%dt
        print "('max thread = ',i0)", nthread
        print "(A)", time_bar

        call current_case%set_current_step(nstep)

        !どのアルゴリズムでも共通して計算させる.
        call calc_gradient_tensor(extents, grid%dv, grid%ds, grid%dx, fluid%velocity, fluid%dudr)

        call solver%predict_pseudo_velocity(extents, grid%ds, grid%dv, grid%dx, dt, re, &
                                            current_case%settings_case%body_force, &
                                            this%v_old, &
                                            fluid%velocity, fluid%mflux_i, fluid%mflux_j, fluid%mflux_k, fluid%dudr, &
                                            current_case%bc_types)

        !中間速度に対する境界条件の適用. ソルバの外部で行う理由は, この処理が共通なため.
        call boundary_condition_velocity(extents, fluid%velocity, grid%dx, current_case%bc_types)

        call solver%calc_corrected_velocity(extents, grid%ds, grid%dv, grid%dx, dt, &
                                           fluid%mflux_i, fluid%mflux_j, fluid%mflux_k, fluid%pressure, fluid%velocity, &
                                           current_case%settings_solver, current_case%settings_case%p_ref, current_case%bc_types, &
                                           sim_diverged)

        !補正された新しい時間段階の速度に対する境界条件の適用.
        call boundary_condition_velocity(extents, fluid%velocity, grid%dx, current_case%bc_types)

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