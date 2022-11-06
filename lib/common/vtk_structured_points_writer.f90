module vtk_structured_points_writer_m
    use floating_point_parameter_m
    use abst_vtk_legacy_writer_m
    use vtk_field_data_m
    implicit none
    private
    character(*),parameter :: DATASET_ =  "DATASET STRUCTURED_POINTS"
    character(*),parameter :: ORIGIN_ =  "ORIGIN"
    character(*),parameter :: SPACING_ =  "SPACING"

    type,extends(base_vtk_legacy_str_grid_t):: vtk_structured_points_t
        real(DP) :: origin(3)
        real(DP) :: spacings(3)
        contains
        procedure :: init => init_
        procedure :: writeout
    end type

    public vtk_structured_points_t
contains

    subroutine init_(this, dimensions, origin, spacings)
        implicit none
        class(vtk_structured_points_t),intent(inout) :: this
        integer(IP),intent(in) :: dimensions(3)
        real(DP),intent(in) :: origin(3)
        real(DP),intent(in) :: spacings(3)
        
        call this%set_dimensions(dimensions)
        this%origin = origin
        this%spacings = spacings

    end subroutine
    
    subroutine writeout(this, path, holder)
        implicit none
        class(vtk_structured_points_t),intent(in) :: this
        type(attrib_data_holder_t),optional,allocatable,intent(in) :: holder(:)
        character(*),intent(in) :: path

        integer(IP) unit, n_, ic_
        
        call writeout_common(this, unit, path)

        if ( .not. this%is_binary ) then
            write(unit,"(A)") DATASET_
            write(unit,"(A,1x,*(g0.6,1x))") "DIMENSIONS", this%dimensions
            write(unit,"(A,1x,*(g0.6,1x))") ORIGIN_, this%origin
            write(unit,"(A,1x,*(g0.6,1x))") SPACING_, this%spacings
    
            if ( present(holder) ) then
                call writeout_common_field_data(this, unit, holder)
            end if    
        end if
        close(unit)
        
        print"(A)", ">>> output complete. path :: ", path
    end subroutine

end module vtk_structured_points_writer_m