module IO_operator_m
    use,intrinsic :: iso_fortran_env, only : int32, real64
    implicit none
    
contains

logical function open_text_file(unit, filename, title, unformatted) result(is_already_opened)
    !! ファイルを開くが, すでに開いてる場合はその装置番号を返す. 
    !! 開いていない場合で、すでに対象のファイルが存在する場合はappendする. 
    implicit none
    integer,intent(inout) :: unit
    character(*),intent(in) :: filename
    character(*),intent(in) :: title
        !!1行目の内容. 
    logical, intent(in), optional :: unformatted
        !!format style. if true, the file is opened as binary file.
        !![default] false

    integer ios
    character(256) errmsg
    character(:),allocatable :: form_
    logical exists

    INQUIRE(file = filename, opened = is_already_opened, number = unit)

    if(is_already_opened) return

    option_arg:block
        logical is_unformatted_

        if ( .not.present(unformatted) ) then
            is_unformatted_ = .false.
        else if(unformatted) then
            is_unformatted_ = .true.
        else
            is_unformatted_ = .false.
        end if
        
        if ( is_unformatted_ ) then
            form_ = "unformatted"
        else
            form_ = "formatted"
        end if
        
    end block option_arg

    INQUIRE(file = filename, exist = exists)
    if ( exists ) then
        open(newunit = unit, file = filename, position = "append", action = "write", iostat = ios, iomsg = errmsg, form=form_)
    else
        open(newunit = unit, file = filename, action = "write", iostat = ios, iomsg = errmsg, form=form_)
    end if

    is_already_opened = .true.

    if(ios /= 0) then
        print "(i0, 1x, A)", ios, errmsg
        error stop
    end if

    if(.not. exists) then
        select case(form_)
        case("formatted")
            write(unit, "(A)") trim(title)
        case("unformatted")
            write(unit) trim(title)
        case default
            error stop "Unexpected command."
        end select
    endif

end function

logical function file_exists(path) result(f_exists)
    !!ファイルが存在するか判定する.
    character(*),intent(in) :: path

    inquire(file=path, exist=f_exists)

end function

!>ファイル名と対応する装置番号を返す．
!>ファイルが存在しない場合は新しく装置を開く．
function get_unit_number(filename) result(unit_number)
    implicit none
    character(*), intent(in) :: filename
        !! ファイル名
    integer(int32) :: unit_number
        !! 戻り値
        !! ファイル名に紐付いた装置番号

    logical :: is_unit_opened

    ! 引数で指定したファイル名のファイルが，装置として開かれているかを問い合わせる
    inquire (file=filename, opened=is_unit_opened, number=unit_number)

    ! 開かれていたらunit_numberが取得されているのでreturnする
    if (is_unit_opened) return

    ! 開かれていなかったらファイルを開き，その装置番号を返す
    open (newunit=unit_number, file=filename)
end function get_unit_number

logical function open_as_binary(unit, path, access, big_endian) result(opened)
    !!ファイルをバイナリで開く. 
    integer,intent(out) :: unit
    character(*),intent(in) :: path
    character(*),intent(in) :: access
        !!`sequential`, `direct`, `stream` のいずれか.
    logical,intent(in),optional :: big_endian
        !!ビッグエンディアンで開く場合. デフォルトは.false.
    
    integer iost
    character(:),allocatable :: endian_
    character(256) errmsg

    if ( present(big_endian) .and. big_endian .eqv. .true.) then
        endian_ = "big_endian"
    else
        endian_ = "little_endian"
    end if

    open(newunit=unit, file=path, access=access, form="unformatted", convert=endian_, iostat=iost, iomsg=errmsg)

    if ( iost /= 0 ) then
        print "('error(',i0') :: ',A)",iost, trim(errmsg)
        opened = .false.
    else
        opened = .true.
    end if

end function


end module IO_operator_m