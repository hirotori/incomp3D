module system_operator_m
    use floating_point_parameter_m
    implicit none
    character(*),parameter,private :: ERRMSG_HEADER_ = "system_operator :: "
contains
    logical function mkdir(dirname) result(success)
        !!ディレクトリを作成する. 既にある場合は端末にメッセージのみ表示される.
        implicit none
        character(*),intent(in) :: dirname

        integer(IP) stat, exitstat
        
        success = .false.
        call execute_command_line_core_("mkdir "//dirname, stat, exitstat)
        ! print*, stat
        if ( stat == 0 ) then
            print "(A,A)", ERRMSG_HEADER_, "directory generated. ", dirname
            success = .true. 
        endif
    end function

    subroutine execute_command_line_core_(command, cmdstat, exitstat)
        implicit none
        character(*),intent(in) :: command
        integer(IP),intent(out) :: cmdstat
        integer(IP),intent(out) :: exitstat
        character(256) :: cmdmsg

        call execute_command_line(command, exitstat=exitstat, cmdstat=cmdstat, cmdmsg=cmdmsg)
        ! print*, cmdstat
        if (cmdstat > 0) then
            print"(A,A)", ERRMSG_HEADER_, trim(cmdmsg)
        else if(cmdstat < 0) then
            print"(A,A)", ERRMSG_HEADER_, "command not supported."
        else
            print"(A,i0)", "Command completed with status ", exitstat
        end if

    end subroutine
end module system_operator_m