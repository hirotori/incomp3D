module abst_vtk_legacy_writer_m
    use floating_point_parameter_m
    use vtk_field_data_m
    implicit none
    private

    character(*),parameter :: HEADER_ = "# vtk DataFile Version 3.0"
    character(*),parameter :: COMMENT_ = "OUTPUT"
    character(*),parameter :: FORMAT_ASCII_ =  "ASCII"
    character(*),parameter :: FORMAT_BINARY_ =  "BINARY"

    character(*),parameter :: ERR_MSG_HEADER = "In writeout ::"
    character(*),parameter :: LF_ = new_line("")
    
    type :: base_vtk_legacy_str_grid_t
        integer(ip) :: dimensions(3)
        integer(ip) :: cell_count_
        logical :: is_binary = .false.

        contains
        procedure save_as_binary
        procedure set_dimensions
        procedure :: writeout
    end type 

    public base_vtk_legacy_str_grid_t, &
           writeout_common, &
           writeout_common_field_data

contains
    subroutine save_as_binary(this)
        !!ファイルをバイナリ形式で保存する.
        !!@note 呼び出されない場合はASCIIで保存される.
        class(base_vtk_legacy_str_grid_t),intent(inout) :: this

        this%is_binary = .true.

    end subroutine

    subroutine set_dimensions(this, dims)
        !!次元を設定する.
        class(base_vtk_legacy_str_grid_t),intent(inout) :: this
        integer(ip),intent(in) :: dims(3)

        this%dimensions = dims
        this%cell_count_ = product(dims - 1)

    end subroutine

    subroutine writeout(this, path, holder)
        implicit none
        class(base_vtk_legacy_str_grid_t),intent(in) :: this
        type(attrib_data_holder_t),optional,allocatable,intent(in) :: holder(:)
        character(*),intent(in) :: path

        error stop ERR_MSG_HEADER//"This proceudure should be implemented in extended types."

    end subroutine

    subroutine writeout_common(this, unit, path)
        !!共通部分の書き出しヘルパー. ファイルを開き, 共通部分まで出力する.
        class(base_vtk_legacy_str_grid_t),intent(in) :: this
        integer(ip),intent(out) :: unit
            !!ファイルへの装置番号.
        character(*),intent(in) :: path
            !!ファイル名.

        character(128) errmsg
        integer(ip) iost
        character(:),allocatable :: form_, access_

        if ( this%is_binary ) then
            form_ = "unformatted"
            access_ = "stream"
        else
            form_ = "formatted"
            access_ = "sequential"
        end if

        open(newunit = unit, file = path, form = form_, access = access_, iostat = iost, iomsg = errmsg)
        if ( iost /= 0 ) then
            print"('code = ',i0, 1x, A)", iost, errmsg 
        end if

        if ( this%is_binary ) then
            
            write(unit) HEADER_//LF_
            write(unit) COMMENT_//LF_
            write(unit) FORMAT_BINARY_//LF_

        else 

            write(unit,"(A)") HEADER_
            write(unit,"(A)") COMMENT_
            write(unit,"(A)") FORMAT_ASCII_

        end if

    end subroutine

    subroutine writeout_common_field_data(this, unit, holders)
        !!フィールドデータの書き出し.
        class(base_vtk_legacy_str_grid_t),intent(in) :: this
        integer(ip),intent(in) :: unit
        type(attrib_data_holder_t),intent(in) :: holders(:)

        integer n_
        if ( this%is_binary ) then
            write(unit) "CELL_DATA "//int32_to_string(this%cell_count_)//LF_
            do n_ = 1, size(holders)
                call holders(n_)%write_data_binary(unit, this%cell_count_)
            end do
        else
            write(unit,"(A)")
            write(unit,"(A,1x,i0)") "CELL_DATA",this%cell_count_
            do n_ = 1, size(holders)
                call holders(n_)%write_data(unit, this%cell_count_)
                write(unit,"(A)")
            end do
        end if



    end subroutine

    function int32_to_string(var) result(char)
        !!整数を文字列へ変換する.
        integer(ip),intent(in) :: var
        character(:),allocatable :: char
        character(32) temp_

        write(temp_, "(i0)") var

        char = trim(adjustl(temp_))

    end function

    function int32_array_to_string(var) result(char)
        !!整数配列を文字列へ変換する.
        integer(ip),intent(in) :: var(:)
        character(:),allocatable :: char
        character(52) temp_
    
        write(temp_, "(*(i0,1x))") var(:)

        char = trim(adjustl(temp_))

    end function
end module abst_vtk_legacy_writer_m