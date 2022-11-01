module vtk_rectilinear_grid_m
    use floating_point_parameter_m
    use abst_vtk_legacy_writer_m
    use vtk_field_data_m
    implicit none
    
    character(*),parameter :: DATASET_ =  "DATASET RECTILINEAR_GRID"
    character(*),parameter :: DIMENSIONS_ =  "DIMENSIONS"
    character(*),parameter :: ERR_MSG_HEADER = "vtk_rectilinear_grid_writer_t ::"

    type, extends(base_vtk_legacy_str_grid_t), public :: vtk_rectilinear_grid_t 
        real(dp),pointer :: x(:) => null(), y(:) => null(), z(:) => null()
            !!頂点座標.
        contains
        procedure set_coordinates
        procedure writeout
    
    end type 


contains

subroutine set_coordinates(this, x, y, z)
    class(vtk_rectilinear_grid_t),intent(inout) :: this
    real(dp),contiguous,target,intent(in) :: x(:), y(:), z(:)
        !!頂点座標.

    integer(ip) imx, jmx, kmx
    imx = ubound(x, dim=1)
    jmx = ubound(y, dim=1)
    kmx = ubound(z, dim=1)

    this%x(1:imx) => x
    this%y(1:jmx) => y
    this%z(1:kmx) => z

end subroutine

subroutine writeout(this, path, holder)
    implicit none
    class(vtk_rectilinear_grid_t),intent(in) :: this
    type(attrib_data_holder_t),optional,allocatable,intent(in) :: holder(:)
    character(*),intent(in) :: path

    integer(ip) unit

    call writeout_common(this, unit, path)
    if (.not. this%is_binary ) then
        write(unit, "(A)") DATASET_
        write(unit, "(A,1x,*(i0,1x))") DIMENSIONS_, this%dimensions
        write(unit, "(A,i0,1x,A)") "X_COORDINATES ", this%dimensions(1), " double"
        write(unit, "(*(g0,1x))") this%x(:)
        write(unit, "(A,i0,1x,A)") "Y_COORDINATES ", this%dimensions(2), " double"
        write(unit, "(*(g0,1x))") this%y(:)
        write(unit, "(A,i0,1x,A)") "Z_COORDINATES ", this%dimensions(3), " double"
        write(unit, "(*(g0,1x))") this%z(:)
    end if

    if ( present(holder) ) then
        call writeout_common_field_data(this, unit, holder)
    end if

end subroutine
    
end module 
