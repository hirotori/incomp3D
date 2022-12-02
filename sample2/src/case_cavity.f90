module case_cavity_m
    use floating_point_parameter_m
    use case_common_m
    use mesh_m
    use fluid_field_m
    use output_m
    implicit none
    type, extends(case_common_t), public :: case_cavity_t
        real(dp),allocatable :: vorticity_(:,:,:,:)
    contains
        procedure add_on_pre_process
        procedure phase_post_process
        procedure phase_writeout
    end type
    
    character(*),parameter :: dir_ = "cavity_result"
        !!計算結果とは別のファイルの出力先ディレクトリ.
    integer(ip),parameter,private :: nwrite_uv = 1000
        !!データ比較用のファイルを出力する間隔. 
        !!@note 定常流れを非定常アルゴリズムで解くので, 逐一確認しながら計算を進めた方が好ましい.

contains
subroutine add_on_pre_process(this, grid, fld)
    !!追加で初期状態を管理したい場合に呼び出す.
    use system_operator_m
    class(case_cavity_t),intent(inout) :: this
    class(equil_mesh_t),intent(inout) :: grid
    type(fluid_field_t),intent(inout) :: fld

    !!セル数が奇数個でない場合ストップさせる.
    integer(ip) extents(3)

    extents(:) = grid%get_extents()

    if ( mod(extents(1)-1, 2) == 0 .or. mod(extents(2)-1, 2) == 0) then
        error stop "case_cavity :: grid points must be even numbers."
    end if

    allocate(this%vorticity_(3,extents(1),extents(2),extents(3)))

    if ( .not. mkdir(dir_) ) error stop "directory creating error." 

end subroutine

subroutine phase_post_process(this, grid, fld, step, time)
    !!追加で初期状態を管理したい場合に呼び出す.
    class(case_cavity_t),intent(inout) :: this
    class(equil_mesh_t),intent(in) :: grid
    type(fluid_field_t),intent(in) :: fld
    integer(ip),intent(in) :: step
    real(dp),intent(in) :: time

    !!キャビティの重心を通る, 鉛直線方向の速度x成分と, 水平線方向の速度y成分をファイルに書き出す.
    integer(ip) x_gc, y_gc
    integer(ip) extents(3)
    character(:),allocatable :: fname
    integer(ip) nstep, unit
    integer(ip) i, j ,k
    real(dp),save :: ke_prev = 0.0d0
    real(dp) :: ke_, dt_
    !簡単のため, セル数は計算領域内に奇数個用意する. 

    extents(:) = grid%get_extents()
    dt_ = this%settings_case%dt
    x_gc = (2 + extents(1))/2
    y_gc = (2 + extents(2))/2

    if ( mod(step,nwrite_uv) == 0 ) then
        fname = get_filename_with_digit_("cavity_u",step,".txt")
        open(newunit=unit, file=dir_//'/'//fname, status="replace")
            write(unit,"(A)") "y u"
            do k = 2, extents(3)
            do j = 2, extents(2)
                write(unit, "(2(g0,1x))") grid%yc(j), fld%velocity(1,x_gc,j,k)
            end do
            end do
        close(unit)

        fname = get_filename_with_digit_("cavity_v",step,".txt")
        open(newunit=unit, file=dir_//'/'//fname, status="replace")
            write(unit,"(A)") "x v"
            do k = 2, extents(3)
            do i = 2, extents(1)
                write(unit, "(2(g0,1x))") grid%xc(i), fld%velocity(2,i,y_gc,k)
            end do
            end do
        close(unit)
    end if

    call calculate_vorticity(this%vorticity_, fld%dudr)

    !定常判定をしてみる.
    !1. 運動エネルギの総和の変化率.
    ke_ = 0.0d0
    do k = 2, extents(3)
    do j = 2, extents(2)
    do i = 2, extents(1)
        ke_ = ke_ + 0.5d0*norm2(fld%velocity(:,i,j,k))**2.0d0*grid%dv(i,j,k)
    end do
    end do
    end do

    print"('>>> kinetic energy')"
    print "('(ke_ - ke_prev)/(ke_prev*dt) = ', g0)", (ke_ - ke_prev)/(ke_*dt_)
    ke_prev = ke_
end subroutine

subroutine phase_writeout(this, grid, fld, nstep)
    !!データの書き出し.
    !!@note 現状, 書き出しサブルーチンは本モジュールのプライベートルーチンとなっている.
    !!@todo 書き出しサブルーチンの隔離.
    class(case_cavity_t),intent(in) :: this
    class(equil_mesh_t),intent(in) :: grid
    type(fluid_field_t),intent(in) :: fld
    integer(ip),intent(in) :: nstep
    character(:),allocatable :: fname

    fname = this%output_directory()
    
    call writeout_(this, grid, fld, fname//"/"//"result_", nstep)

end subroutine

subroutine writeout_(this, grid, fld, basename, current_step)
    !!データを書き出す. フォーマットはvtkファイルとなっている.
    !!@note データ節約のため格子座標データを必要としないvtk structured points フォーマットを使用している.
    class(case_cavity_t),intent(in) :: this
    type(equil_mesh_t),intent(in) :: grid
    type(fluid_field_t),intent(in) :: fld
    character(*),intent(in) :: basename
    integer(ip),intent(in) :: current_step

    type(attrib_data_holder_t),allocatable :: holders(:)

    allocate(holders(3))
    
    holders(1)%name = "Velocity"
    call holders(1)%register_vector(fld%velocity, [2,2,2], grid%get_extents())
    holders(2)%name = "Pressure"
    call holders(2)%register_scalar(fld%pressure, [2,2,2], grid%get_extents())
    holders(3)%name = "Vorticity"
    call holders(3)%register_vector(this%vorticity_, [2,2,2], grid%get_extents())

    call writeout_single_vtk_recti_grid(basename, current_step, grid%get_extents(), &
        grid%xp, grid%yp, grid%zp, holders=holders)
    
end subroutine
    
end module case_cavity_m