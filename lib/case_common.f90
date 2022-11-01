module case_common_m
    !/////////////////////////////////////////////////////////////////
    !! version: 1.0.0
    !! author: T.Ikeda
    !! summary:
    !! - 解析のループ前処理, メインループの後処理, データ書きだしに関するモジュール.
    !! @note
    !! 追加の処理が必要な場合は, 本モジュール内のクラスをオーバーライドする.
    !! @endnote
    !! @todo
    !! リファクタリング
    !! @endtodo
    !/////////////////////////////////////////////////////////////////
    use floating_point_parameter_m, only : ip, dp
    use mesh_reader_m
    use mesh_m
    use rectilinear_mesh_m
    use fluid_field_m
    use setting_parameter_m
    use boundary_condition_m, only : bc_t, set_bc_type
    use vtk_field_data_m
    use output_m
    use system_operator_m
    implicit none

    private

    type case_common_t
        !/////////////////////////////////////////////////////////////////
        !! version: 1.0.0
        !! author: T.Ikeda
        !! summary:
        !! - 解析ケースの前処理や後処理, データ保存処理などを担う.
        !/////////////////////////////////////////////////////////////////
        character(:),allocatable :: config_file
            !!計算設定ファイル名.
        type(case_setting) :: settings_case
            !!解析ケースに関連する設定パラメータ.
        type(bc_t) :: bc_types(6)
            !!各境界面の境界条件および境界値（速度, 圧力）を扱う.
        character(:),allocatable,private :: outdir_
            !!出力ディレクトリ.
        integer(ip),private :: current_step_ = 0
        logical,private :: is_equil_mesh = .true.
        contains
        procedure,non_overridable :: get_current_step
        procedure,non_overridable :: get_current_time
        procedure,non_overridable :: set_current_step
        procedure,non_overridable :: set_output_directory
        procedure,non_overridable :: get_current_output_directory
        procedure,non_overridable :: phase_pre_process
        procedure,non_overridable :: uses_an_equil_mesh
        procedure :: add_on_pre_process
        procedure :: phase_post_process
        procedure :: phase_writeout
        procedure :: check_flow_field
        procedure,non_overridable :: process_diverged
    end type

    public case_common_t

contains

pure integer(ip) function get_current_step(this)
    !!現在の計算ステップを取得する.
    class(case_common_t),intent(in) :: this

    get_current_step = this%current_step_

end function

pure real(dp) function get_current_time(this)
    !!現在の計算ステップを取得する.
    class(case_common_t),intent(in) :: this

    get_current_time = this%current_step_*this%settings_case%dt

end function

subroutine set_output_directory(this, path)
    !!出力ディレクトリを指定する.
    class(case_common_t),intent(inout) :: this
    character(*),intent(in) :: path
        !!計算結果ファイルを出力するディレクトリのパス.

    this%outdir_ = path


end subroutine

function get_current_output_directory(this) result(path)
    !!現在の出力先ディレクトリを返す.
    class(case_common_t),intent(in) :: this
    character(:),allocatable :: path

    path = this%outdir_

end function

subroutine set_current_step(this, current_step)
    !!現在の計算ステップを記録する.
    class(case_common_t),intent(inout) :: this
    integer(ip),intent(in) :: current_step

    this%current_step_ = current_step

end subroutine

logical function uses_an_equil_mesh(this)
    class(case_common_t),intent(in) :: this

    uses_an_equil_mesh = this%is_equil_mesh

end function

subroutine phase_pre_process(this, grid, fld)
    !!計算開始前の処理を行う. 
    !!計算の設定ロードや格子の初期化, 流体クラスの初期化を担う.
    !!@note 追加の処理が必要な場合は `add_on_pre_process` を使用する.
    !!@warning このメンバはオーバーライド出来ない.
    class(case_common_t),intent(inout) :: this
    class(equil_mesh_t),allocatable,intent(inout) :: grid
    type(fluid_field_t),intent(inout) :: fld

    integer(ip) imx, jmx, kmx
    real(dp) lengths(3)
    integer(ip) :: bc_ids(2,6)
    real(dp) :: bc_properties(4,6)

    real(dp),allocatable :: x(:), y(:), z(:)
    
    this%config_file = "config.txt"
    this%outdir_ = "output"
    if ( .not. mkdir(this%outdir_) )  print "(A)", "command skipped."
    call read_config(this%config_file, this%settings_case, bc_ids, bc_properties)
    
    call get_mesh_from_file(this%settings_case%grid_file_name, imx, jmx, kmx, lengths, x, y, z)
    !座標配列が割り当てられているかどうかで分岐する.
    if ( allocated(x) .and. allocated(y) .and. allocated(z)) then
        allocate(rectilinear_mesh_t :: grid)
        call grid%init(imx, jmx, kmx, x, y, z)
        this%is_equil_mesh = .false.
    else
        allocate(equil_mesh_t :: grid)
        call grid%init(imx, jmx, kmx, lengths)
    end if

    call set_bc_type(this%bc_types, bc_ids, bc_properties, grid%get_extents())
    call fld%init(imx, jmx, kmx, this%settings_case%u_ic, this%settings_case%p_ref, &
                  this%settings_case%reynolds_number, this%settings_case%body_force)

    call this%add_on_pre_process(grid, fld)
    
end subroutine

subroutine add_on_pre_process(this, grid, fld)
    !!追加で初期状態を管理したい場合に呼び出す.
    class(case_common_t),intent(inout) :: this
    class(equil_mesh_t),intent(inout) :: grid
    type(fluid_field_t),intent(inout) :: fld

    print "('Add-On process is being called from **case_common_t**')"
end subroutine

subroutine phase_post_process(this, grid, fld)
    !!現時間段階の計算が終わった後に呼び出される. 
    class(case_common_t),intent(inout) :: this
    class(equil_mesh_t),intent(in) :: grid
    type(fluid_field_t),intent(in) :: fld

    print "('post process is being called from **case_common_t** ')"


end subroutine

subroutine phase_writeout(this, grid, fld, nstep)
    !!データの書き出し.
    !!@note 現状, 書き出しサブルーチンは本モジュールのプライベートルーチンとなっている.
    !!@todo 書き出しサブルーチンの隔離.
    class(case_common_t),intent(in) :: this
    class(equil_mesh_t),intent(in) :: grid
    type(fluid_field_t),intent(in) :: fld
    integer(ip),intent(in) :: nstep
    ! character(:),allocatable :: fname
    
    call writeout_(this, grid, fld, this%outdir_//"/"//"result_", nstep)

end subroutine

subroutine writeout_(this, grid, fld, basename, current_step)
    !!データを書き出す. フォーマットはvtkファイルとなっている.
    !!@note データ節約のため等間隔メッシュの場合格子座標データを必要としないvtk structured points フォーマットを使用している.
    class(case_common_t),intent(in) :: this
    class(equil_mesh_t),intent(in) :: grid
    type(fluid_field_t),intent(in) :: fld
    character(*),intent(in) :: basename
    integer(ip),intent(in) :: current_step

    type(attrib_data_holder_t),allocatable :: holders(:)

    allocate(holders(2))
    
    holders(1)%name = "Velocity"
    call holders(1)%register_vector(fld%velocity, [2,2,2], grid%get_extents())
    holders(2)%name = "Pressure"
    call holders(2)%register_scalar(fld%pressure, [2,2,2], grid%get_extents())

    if ( this%uses_an_equil_mesh() ) then
        call writeout_single_vtk_str_points(basename, current_step, grid%get_extents(), & 
             [real(dp) :: 0, 0, 0], grid%get_equil_dx(), holders)
    else
        call writeout_single_vtk_recti_grid(basename, current_step, grid%get_extents(), &
                                            grid%xp, grid%yp, grid%zp, holders=holders)
    end if
    
end subroutine


subroutine check_flow_field(this, grid, fld)
    !!現段階の流れ場計算結果の確認.
    !!@todo オーバーライド属性の検討. モニタできることは決まっているので不要? check_resultをパブリックにしておけば十分か.
    class(case_common_t),intent(in) :: this
    class(equil_mesh_t),intent(in) :: grid
    type(fluid_field_t),intent(in) :: fld

    integer(ip) extents(3)

    extents(:) = grid%get_extents()
    !check conservation law
    call check_result(extents(1), extents(2), extents(3), grid%dsx, grid%dsy, grid%dsz, grid%dv, &
                      fld%velocity, fld%pressure)    

end subroutine

subroutine check_result(imx, jmx, kmx, dsi, dsj, dsk, dv, v, p)
    !!速度, 圧力の値域, および保存性のチェック.
    implicit none
    integer(IP) imx, jmx, kmx
    real(DP),intent(in) :: dsi(2:,2:), dsj(2:,2:), dsk(2:,2:)
    real(DP),intent(in) :: dv(2:,2:,2:)
    real(DP),intent(in) :: v(:,:,:,:), p(:,:,:)
    
    integer(IP) i, j, k

    real(DP) surf_q(6), div_u_c, uf(6)
    real(dp) ds(3)
    
    print "('***checking current result***')"

    print "(g0, ' <= p <= ', g0)", minval(p(2:imx,2:jmx,2:kmx)), maxval(p(2:imx,2:jmx,2:kmx))
    print "(g0, ' <= u <= ', g0)", minval(v(1,2:imx,2:jmx,2:kmx)), maxval(v(1,2:imx,2:jmx,2:kmx))
    print "(g0, ' <= v <= ', g0)", minval(v(2,2:imx,2:jmx,2:kmx)), maxval(v(2,2:imx,2:jmx,2:kmx))
    print "(g0, ' <= w <= ', g0)", minval(v(3,2:imx,2:jmx,2:kmx)), maxval(v(3,2:imx,2:jmx,2:kmx))
    
    surf_q(:) = 0.0_dp
    print "(':-- SURFACE FLOW --:')" !境界面の流束積分. 
    do k = 2, kmx
    do j = 2, jmx
        ds(1) = dsi(j,k)
        surf_q(1) = surf_q(1) + 0.5_dp*(v(1,2,j,k) + v(1,1,j,k))*ds(1)
        surf_q(2) = surf_q(2) + 0.5_dp*(v(1,imx,j,k) + v(1,imx+1,j,k))*ds(1)
    end do    
    end do

    do k = 2, kmx
    do i = 2, imx
        ds(2) = dsj(i,k)
        surf_q(3) = surf_q(3) + 0.5_dp*(v(2,i,1,k) + v(2,i,2,k))*ds(2)
        surf_q(4) = surf_q(4) + 0.5_dp*(v(2,i,jmx,k) + v(2,i,jmx+1,k))*ds(2)
    end do    
    end do

    do j = 2, jmx
    do i = 2, imx
        ds(3) = dsk(i,j)
        surf_q(5) = surf_q(5) + 0.5_dp*(v(3,i,j,1) + v(3,i,j,2))*ds(3)
        surf_q(6) = surf_q(6) + 0.5_dp*(v(3,i,j,kmx) + v(3,i,j,kmx+1))*ds(3)
    end do    
    end do

    print "('xmin-xmax ', *(g0,1x,:,':'))", surf_q(1:2)
    print "('ymin-ymax ', *(g0,1x,:,':'))", surf_q(3:4)
    print "('zmin-zmax ', *(g0,1x,:,':'))", surf_q(5:6)

    print "(':-- Divergence --:')"
    ! div_u_f = 0.0_dp
    div_u_c = 0.0_dp
    do k = 2, kmx
    do j = 2, jmx
    do i = 2, imx
        ds(1) = dsi(j,k)
        ds(2) = dsj(i,j)
        ds(3) = dsk(i,k)
        uf(1) = -0.5_dp*(v(1,i-1,j  ,k  ) + v(1,i  ,j  ,k  ))*ds(1)
        uf(2) =  0.5_dp*(v(1,i  ,j  ,k  ) + v(1,i+1,j  ,k  ))*ds(1)
        uf(3) = -0.5_dp*(v(2,i  ,j-1,k  ) + v(2,i  ,j  ,k  ))*ds(2)
        uf(4) =  0.5_dp*(v(2,i  ,j  ,k  ) + v(2,i  ,j+1,k  ))*ds(2)
        uf(5) = -0.5_dp*(v(3,i  ,j  ,k-1) + v(3,i  ,j  ,k  ))*ds(3)
        uf(6) =  0.5_dp*(v(3,i  ,j  ,k  ) + v(3,i  ,j  ,k+1))*ds(3)
        div_u_c = div_u_c + sum(uf)
    end do    
    end do
    end do
    print "('sum of div(u) ', g0)", div_u_c !ゼロに近いはず.

end subroutine


subroutine process_diverged(this, grid, fld, current_step)
    !!計算が発散したときの処理. 
    !!@note オーバーライド不可.
    class(case_common_t),intent(in) :: this
    class(equil_mesh_t),intent(in) :: grid
    type(fluid_field_t),intent(in) :: fld
    integer(ip),intent(in) :: current_step

    call writeout_(this, grid, fld, this%outdir_//"/"//"diverged_", current_step)
end subroutine
    
end module case_common_m