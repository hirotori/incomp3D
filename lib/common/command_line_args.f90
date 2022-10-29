module command_line_args_m
    !/////////////////////////////////////////////////////////////////
    !! version: 1.0.0
    !! author:  T.Ikeda
    !! summary:
    !! - コマンドライン引数を取得する関数
    !! @note
    !! 構造体`args_t`の配列に格納して返す. インデックスは1から. 1番目は実行ファイル名.
    !! @endnote
    !! @warning
    !! Fortran2018の機能を用いるため, 一部コンパイラでは動作しない. `ifort`でのみ動作する. gfortranではコンパイルエラー.
    !! @endwarning
    !/////////////////////////////////////////////////////////////////
    implicit none

    type args_t
        character(:),allocatable,public :: v
    end type

contains

function argv() result(args)
    !!コマンドライン引数を返す.
    !!インデックスは1から.
    type(args_t),allocatable :: args(:)
    integer :: argc
    integer i_, length_, stat_
    character(255) err_msg_

    !|## 例
    !```Fortran
    !program test
    !use command_line_args_m
    !implicit none
    !integer n
    !type(args_t),allocatable :: args(:)
    !
    !args = argv()
    !
    !do n = 1, size(args)
    !   print*, args(n)%v
    !end do
    !end program test
    !```

    argc = command_argument_count() 
    allocate(args(0:argc))
    do i_ = 0, argc
#ifdef __GFORTRAN__
        call get_command_argument(number=i_, length=length_, status=stat_)
        if (stat_ /= 0) error stop "Some error occured in counting argument."
#else
        call get_command_argument(number=i_, length=length_, status=stat_, errmsg=err_msg_)
        if ( stat_ /= 0 ) then
                print*, trim(err_msg_)
                error stop "Some error occured in counting argument."
        end if
#endif
        allocate(character(length_) :: args(i_)%v)
#ifdef __GFORTRAN__
        call get_command_argument(number=i_, value=args(i_)%v, status=stat_)
        if (stat_ /= 0) error stop "Some error occured in counting argument."
#else
        call get_command_argument(number=i_, value=args(i_)%v, status=stat_, errmsg=err_msg_)
        if ( stat_ /= 0 ) then
                print*, trim(err_msg_)
                error stop "Some error occured in getting command argument."
        end if
#endif

    end do
end function
    
end module command_line_args_m

! program test
!     use command_line_args_m
!     implicit none
!     integer n
!     type(args_t),allocatable :: args(:)

!     args = argv()

!     do n = 1, size(args)
!         print*, args(n)%v
!     end do

! end program test