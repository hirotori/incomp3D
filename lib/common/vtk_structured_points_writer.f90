module vtk_structured_points_writer_m
    use floating_point_parameter_m
    use vtk_field_data_m
    implicit none
    private
    character(*),parameter :: HEADER_ = "# vtk DataFile Version 3.0"
    character(*),parameter :: COMMENT_ = "OUTPUT"
    character(*),parameter :: FORMAT_ASCII_ =  "ASCII"
    character(*),parameter :: FORMAT_BINARY_ =  "BINARY"
    character(*),parameter :: DATASET_ =  "DATASET STRUCTURED_POINTS"
    character(*),parameter :: DIMENSIONS_ =  "DIMENSIONS"
    character(*),parameter :: ORIGIN_ =  "ORIGIN"
    character(*),parameter :: SPACING_ =  "SPACING"

    character(*),parameter :: ERR_MSG_HEADER = "field_file_t ::"

    type vtk_structured_points_t
        integer(IP) :: dimensions(3)
        real(DP) :: origin(3)
        real(DP) :: spacings(3)
        integer(IP),private :: cell_count_ 
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
        
        this%dimensions = dimensions
        this%origin = origin
        this%spacings = spacings

        this%cell_count_ = product(dimensions(:)-1)

    end subroutine
    
    subroutine writeout(this, path, holder)
        implicit none
        class(vtk_structured_points_t),intent(in) :: this
        type(attrib_data_holder_t),optional,allocatable,intent(in) :: holder(:)
        character(*),intent(in) :: path

        integer(IP) unit, n_, ic_

        if(present(holder) .and. .not. allocated(holder)) error stop ERR_MSG_HEADER//" unallocated holder recieved."
        open(newunit=unit, file=path, status="replace", action="write", form="formatted")
        write(unit,"(A)") HEADER_
        write(unit,"(A)") COMMENT_
        write(unit,"(A)") FORMAT_ASCII_
        write(unit,"(A)") DATASET_
        write(unit,"(A,1x,*(g0.6,1x))") DIMENSIONS_, this%dimensions
        write(unit,"(A,1x,*(g0.6,1x))") ORIGIN_, this%origin
        write(unit,"(A,1x,*(g0.6,1x))") SPACING_, this%spacings

        if ( allocated(holder) ) then
            write(unit,"(A)")
            write(unit,"(A,1x,i0)") "CELL_DATA",this%cell_count_
            do n_ = 1, size(holder)
                call holder(n_)%write_data(unit, this%cell_count_)
                write(unit,"(A)")
            end do
        end if

        close(unit)
        
        print*, "output complete. path :: ", path
    end subroutine

end module vtk_structured_points_writer_m