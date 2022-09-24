module time_monitor_m
    !/////////////////////////////////////////////////////////////////
    !! version: 1.0.0
    !! author: T.Ikeda
    !! summary:
    !! - 時刻操作を取り扱うモジュール
    !! - アプリケーションが起動する国での標準時で値を取り扱う.
    !/////////////////////////////////////////////////////////////////
    implicit none
    type date_t
        !! 現在時刻を取り扱うクラス.
        integer Year
        integer Month
        integer Day
        integer Hour
        integer Minute
        real(8) Second
        contains
        procedure print_self
        procedure is_leap_year
    end type date_t
contains
    !--------------------------------------------------------------------------
    !class date_t
    !--------------------------------------------------------------------------

    subroutine print_self(this)
        implicit none
        class(date_t),intent(in) :: this

        write(*,"(i0, '/',i0, 1x, i0, ':', i0, ':', f8.4 )") &
        this%Month, this%Day, this%Hour, this%Minute, this%Second

    end subroutine print_self

    logical function is_leap_year(this)
        class(date_t),intent(in) :: this

        if (mod(this%Year,4)==0 .or. (mod(this%Year,100)==0 .and. mod(this%Year,400)/=0)) then
            is_leap_year = .true.
        else
            is_leap_year = .false.
        end if
    end function is_leap_year
    
    !> dateクラスのオブジェクトを返す．呼び出した時点での月日と時刻（秒単位）を持つ
    function get_current_date_and_time() result(current_time)
        use,intrinsic :: iso_fortran_env
        implicit none
        type(date_t) :: current_time

        enum, bind(c)
            enumerator :: Year = 1
            enumerator :: Month
            enumerator :: Day
            enumerator :: TimeDifference_min
            enumerator :: Hour
            enumerator :: Minute
            enumerator :: Second
            enumerator :: Millisecond
        end enum

        INTEGER T0(8)

        CALL DATE_AND_TIME(VALUES=T0)
        current_time%Year = T0(Year)
        current_time%Month = T0(Month)
        current_time%Day = T0(Day)
        ! current_time%timeHMS_sec = 3600.0 * T0(Hour) + 60.0 * T0(Minute) + 1.0 * T0(Second)+0.001 * T0(Millisecond)
        current_time%Hour = T0(Hour)
        current_time%Minute = T0(Minute)
        current_time%Second = T0(Second) + 0.001 * T0(Millisecond)
    
    end function get_current_date_and_time

    !--------------------------------------------------------------------------
    ! operater
    !--------------------------------------------------------------------------

    !> 経過時間の計算.
    !! 計算が1年以内で終わることを想定. 
    function calc_elapsed_time(date_to, date_from) result(res)
        implicit none
        type(date_t), intent(in):: date_to, date_from
        type(date_t) res

        !月ごとの日数. 
        integer :: days(12) = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

        res%Month = date_to%Month - date_from%Month

        ! 現時点の日付が開始時より小さいとき，開始月の日数を足してから引く. 同時にres%Monthから1引く
        if(date_to%Day < date_from%Day) then

            ! 開始年月が閏年の2月の場合の処理
            if ( date_from%is_leap_year() .and. date_from%Month == 2) then
                res%Day = date_to%Day + 29 - date_from%Day
            else
                res%Day = date_to%Day + days(date_from%Month) - date_from%Day
            end if
            
            res%Month = res%Month - 1
        
        else
            res%Day = date_to%Day - date_from%Day
        endif

        ! 現時刻の時間が開始時より小さい場合24足してから引く. このときres%Dayから1引く
        if (date_to%Hour < date_from%Hour ) then
            res%Hour = date_to%Hour + 24 - date_from%Hour
            res%Day = res%Day - 1
        else
            res%Hour = date_to%Hour - date_from%Hour
        end if

        ! 現時刻の分が開始時の分より小さい場合60を足して引く．このときres%hourから１引く
        if(date_to%Minute < date_from%Minute) then
            res%Minute = date_to%Minute + 60 - date_from%Minute
            res%Hour   = res%Hour - 1
        else
            res%Minute = date_to%Minute - date_from%Minute
        endif
        
        ! 現時刻の秒が開始時の秒より小さい場合60を足して引く。このときres%Minuteから１引く
        if(date_to%Second < date_from%Second) then
            res%Second = date_to%Second + 60 - date_from%Second
            res%Minute   = res%Minute - 1
        else
            res%Second = date_to%Second - date_from%Second
        endif

    end function calc_elapsed_time

end module time_monitor_m