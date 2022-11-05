module output_m
    use floating_point_parameter_m
    use vtk_field_data_m
    use vtk_structured_points_writer_m
    use vtk_rectilinear_grid_m
    use mesh_reader_m, only : mesh_form_equil, mesh_form_rectil
    implicit none

    interface writeout_as_msh_format
        module procedure :: writeout_as_msh_equil_format
        module procedure :: writeout_as_msh_rectil_format
    end interface

contains
function get_filename_with_digit_(filename_no_ext, digit, ext) result(filename_)
    !!番号付きファイル名を返す. 数字は0埋めされない.
    implicit none
    character(*),intent(in) :: filename_no_ext
        !!拡張子のないファイル名.
    integer,intent(in) :: digit
        !!ファイル名に付けたい番号.
    character(*),intent(in) :: ext
        !!拡張子

    character(:),allocatable :: filename_

    character(32) digit_with_ext

    write(digit_with_ext, "(i0,A)") digit, ext
    filename_ = trim(filename_no_ext) // trim(digit_with_ext)

end function

subroutine writeout_single_vtk_str_points(basename, nstep, extents, origin, spacings, holders)
    !!データをvtk structured pointsフォーマットで書き出す.
    character(*),intent(in) :: basename
    integer(ip),intent(in) :: nstep
    integer(ip),intent(in) :: extents(3)
    real(dp),intent(in) :: origin(3)
    real(dp),intent(in) :: spacings(3)
    type(attrib_data_holder_t),intent(in),allocatable,optional :: holders(:)

    character(:),allocatable :: fname
    type(vtk_structured_points_t) vtp_writer

    fname = get_filename_with_digit_(basename, nstep, ".vtk")
    call vtp_writer%init(extents, origin, spacings)
    call vtp_writer%writeout(fname, holders)

end subroutine

subroutine writeout_single_vtk_recti_grid(basename, nstep, extents, x, y, z, holders, save_as_bin)
    !!データをvtk rectilinear gridフォーマットで書き出す.
    character(*),intent(in) :: basename
    integer(ip),intent(in) :: nstep
    integer(ip),intent(in) :: extents(3)
    real(dp),intent(in) :: x(:), y(:), z(:)
    type(attrib_data_holder_t),intent(in),allocatable,optional :: holders(:)
    logical,intent(in),optional :: save_as_bin

    character(:),allocatable :: fname
    type(vtk_rectilinear_grid_t) writer
    ! logical is_bin

    fname = get_filename_with_digit_(basename, nstep, ".vtk")
    call writer%set_dimensions(extents)
    call writer%set_coordinates(x, y, z)
    if ( present(save_as_bin) ) then
        if ( save_as_bin ) call writer%save_as_binary()
    end if
    call writer%writeout(fname, holders)

end subroutine

subroutine writeout_as_msh_equil_format(basename, extents, x, y, z)
    !!等間隔格子をオリジナルフォーマット(.msh)で書き出す.
    character(*),intent(in) :: basename
    integer(ip),intent(in) :: extents(3)
    real(dp),intent(in) :: x, y, z
        !!領域の寸法.

    integer unit

    open(newunit=unit, file=trim(basename)//".msh", status="replace")
    write(unit, "(A)") mesh_form_equil
    write(unit, "(3(i0,1x))") extents(1:3) 
    write(unit, "(3(g0,1x))") x, y, z
    close(unit) 

end subroutine

subroutine writeout_as_msh_rectil_format(basename, extents, x, y, z)
    !!不等間隔格子をオリジナルフォーマット(.msh)で書き出す.
    character(*),intent(in) :: basename
    integer(ip),intent(in) :: extents(3)
    real(dp),intent(in) :: x(:), y(:), z(:)
        !!節点座標配列.
    integer i, j, k
    integer unit

    open(newunit=unit, file=trim(basename)//".msh", status="replace")
    write(unit, "(A)") mesh_form_rectil
    write(unit, "(3(i0,1x))") extents(1:3) 
    do k = 1, extents(3)
    do j = 1, extents(2)
    do i = 1, extents(1)
        write(unit, "(3(g0,1x))") x(i), y(j), z(k)
    end do
    end do    
    end do

    close(unit) 

end subroutine

end module output_m