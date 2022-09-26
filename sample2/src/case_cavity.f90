module case_cavity_m
    use floating_point_parameter_m
    use case_common_m
    use mesh_m
    use fluid_field_m
    use output_m, only : get_filename_with_digit_
    implicit none
    type, extends(case_common_t), public :: case_cavity_t
        
    contains
        procedure add_on_pre_process
        procedure phase_post_process
    end type
    
contains
subroutine add_on_pre_process(this, grid, fld)
    !!追加で初期状態を管理したい場合に呼び出す.
    class(case_cavity_t),intent(inout) :: this
    type(rectilinear_mesh_t),intent(inout) :: grid
    type(fluid_field_t),intent(inout) :: fld

    !!セル数が奇数個でない場合ストップさせる.
    integer(ip) extents(3)

    extents(:) = grid%get_extents()

    if ( mod(extents(1)-1, 2) == 0 .or. mod(extents(2)-1, 2) == 0) then
        error stop "case_cavity :: grid points must be odd numbers."
    end if

end subroutine

subroutine phase_post_process(this, grid, fld)
    !!追加で初期状態を管理したい場合に呼び出す.
    class(case_cavity_t),intent(in) :: this
    type(rectilinear_mesh_t),intent(in) :: grid
    type(fluid_field_t),intent(in) :: fld

    !!キャビティの重心を通る, 鉛直線方向の速度x成分と, 水平線方向の速度y成分をファイルに書き出す.
    integer(ip) x_gc, y_gc
    integer(ip) extents(3)
    character(:),allocatable :: fname
    integer(ip) nstep, unit
    integer(ip) i, j ,k
    !簡単のため, セル数は計算領域内に奇数個用意する. 

    nstep = this%get_current_step()
    extents(:) = grid%get_extents()

    x_gc = (2 + extents(1))/2
    y_gc = (2 + extents(2))/2

    if ( mod(nstep,100) == 0 ) then
        fname = get_filename_with_digit_("cavity_u",nstep,".txt")
        open(newunit=unit, file=fname, status="replace")
            write(unit,"(A)") "y u"
            do k = 2, extents(3)
            do j = 2, extents(2)
                write(unit, "(2(g0,1x))") grid%rc(2,x_gc,j,k), fld%velocity(1,x_gc,j,k)
            end do
            end do
        close(unit)

        fname = get_filename_with_digit_("cavity_v",nstep,".txt")
        open(newunit=unit, file=fname, status="replace")
            write(unit,"(A)") "x v"
            do k = 2, extents(3)
            do i = 2, extents(1)
                write(unit, "(2(g0,1x))") grid%rc(1,i,y_gc,k), fld%velocity(2,i,y_gc,k)
            end do
            end do
        close(unit)
    end if

end subroutine
    
end module case_cavity_m