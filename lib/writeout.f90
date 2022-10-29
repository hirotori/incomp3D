module output_m
    use floating_point_parameter_m
    use vtk_field_data_m
    use vtk_structured_points_writer_m
    use vtk_rectilinear_grid_m
    implicit none

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

subroutine writeout_single_vtk_recti_grid(basename, nstep, extents, rp, holders, save_as_bin)
    !!データをvtk rectilinear gridフォーマットで書き出す.
    character(*),intent(in) :: basename
    integer(ip),intent(in) :: nstep
    integer(ip),intent(in) :: extents(3)
    real(dp),intent(in) :: rp(:,:,:,:)
    type(attrib_data_holder_t),intent(in),allocatable,optional :: holders(:)
    logical,intent(in),optional :: save_as_bin

    character(:),allocatable :: fname
    type(vtk_rectilinear_grid_t) writer
    ! logical is_bin

    fname = get_filename_with_digit_(basename, nstep, ".vtk")
    call writer%set_dimensions(extents)
    call writer%set_coordinates(rp)
    if ( present(save_as_bin) ) then
        if ( save_as_bin ) call writer%save_as_binary()
    end if
    call writer%writeout(fname, holders)

end subroutine

end module output_m