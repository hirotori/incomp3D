module setting_parameter_m
    ! use,intrinsic :: iso_fortran_env, only : IP => int32, DP => real64
    use floating_point_parameter_m
    use boundary_condition_m, only : bc_fix_grad, bc_fix_val, bc_outlet_unsteady, bc_periodic, bc_periodic_buffer
    implicit none
    private

    character(*),parameter :: current_file_version = "1.4"

    type case_setting
        character(128) :: grid_file_name
        integer(ip) :: nstart = 1
        integer(ip) :: nend = 10000
        integer(ip) :: nwrite = 100
        real(dp) :: dt = 0.001_dp
        real(dp) :: p_ref = 0.0_dp
        real(dp) :: u_ic(3) = 0.0_dp
        real(dp) :: body_force(3) = 0.0_dp
        real(dp) :: reynolds_number = 100.0_dp
    end type
    
    public case_setting, read_config

contains
subroutine read_config(fname, c_setting, bc_ids, bc_properties)
    character(*),intent(in) :: fname
    type(case_setting),intent(out) :: c_setting
    integer(ip),intent(out) :: bc_ids(2,6)
        !!速度/圧力境界条件. 速度, 圧力の順番.
    real(dp),intent(out) :: bc_properties(4,6)

    print "(A)", ">>> Loading config file .."


    call read_config_core(fname, c_setting, bc_ids, bc_properties)

    !TODO リファクタリング. namelistにして各構造体は別ファイル扱いとする??

end subroutine

subroutine read_config_core(fname, setting_c, bc_ids, bc_properties)
    character(*),intent(in) :: fname
    type(case_setting),intent(out) :: setting_c
    integer(ip),intent(out) :: bc_ids(2,6)
    real(dp),intent(out) :: bc_properties(4,6)

    integer(IP) unit, iost, l
    character(128) tmp_

    open(newunit = unit, file = fname, action = "read", status = "old", iostat = iost, iomsg = tmp_)
    if ( iost /= 0 ) then
        call create_sample_()
        error stop trim(tmp_)
    end if

    read(unit,*,err=99) tmp_
    !version判定.
    block
        character(:),allocatable :: ver_
        integer(ip) mj, mn
        ver_ = trim(adjustl(tmp_))
        if ( ver_ /= current_file_version ) then
            print "(A,A,A,A)", &
            "config file version wrong. ver.",current_file_version, " recquired, but current file is ", ver_
            error stop
        end if
    end block

    block
        integer pos_end_
        read(unit,"(A)",err=99) tmp_ 
        !コメント行も含まれるので, `!`以降を捨てる
        pos_end_ = index(tmp_, "!")
        if (pos_end_ == 0) pos_end_ = len(tmp_)
        setting_c%grid_file_name = trim(adjustl(tmp_(1:pos_end_-1)))
    end block
    
    read(unit,*,err=99) setting_c%nstart, setting_c%nend, setting_c%nwrite
    read(unit,*,err=99) setting_c%dt
    read(unit,*,err=99) setting_c%reynolds_number
    read(unit,*,err=99) setting_c%body_force(:)
    read(unit,*,err=99) setting_c%p_ref
    read(unit,*,err=99) setting_c%u_ic(:)
    do l = 1, 6
        read(unit,*,err=99) bc_ids(:,l), bc_properties(:,l)        
    end do
    print "('======= from config file =======')"
    print "('time step : ',i0, '<= n <=', i0)", setting_c%nstart, setting_c%nend
    print "('delta t   : ',g0.5)", setting_c%dt
    print "('Re        : ',g0)", setting_c%reynolds_number
    print "('Force  : (', g0.2, ', ', g0.2, ', ', g0.2, ') ')", setting_c%body_force(:)
    close(unit)

    return

99 continue

    print "('ERROR OCCURED IN READ_CONFIG')"
    call create_sample_()
    print "('Please check sample file.')"
    error stop

end subroutine

subroutine create_sample_()
    integer(ip) unit
    integer(ip) l
    character(8) :: label_(6) = ["i = 1   ", "i = imax", "j = 1   ", "j = jmax", "k = 1   ", "k = kmax"] 
    character(*),parameter :: dummy_flag = "T"
    open(newunit = unit, file = "config_sample.txt", status = "replace")
        write(unit, "(A, 3x, '! current config file version.')") current_file_version
        write(unit, "(A, 3x, '! grid file name. extention must be "".msh"".')") "path/to/your/grid_file_name.msh"
        write(unit, "(3(i0,1x), ' !time step (start,end,write_interval')") 0, 1000, 100 
        write(unit, "(g0.3      , ' !time spacing')") 0.001d0         
        write(unit, "(g0.2      , ' !Reynolds Number')") 100.0d0
        write(unit, "(3(g0.2,1x), ' !body force (x,y,z) ')") 0.0d0, 0.0d0, 0.0d0             
        write(unit, "(g0.2      , ' !Pressure reference (also for initial condition)')") 1.0d0
        write(unit, "(3(g0.2,1x), ' !initial condition for velocity (x,y,z)     ')") 0.0d0, 0.0d0, 0.0d0 
        do l = 1, 6
            write(unit, "(2(i0,1x),4(g0.2,1x),' !B.C. for ',A)") 1, 1, 0.0d0, 0.0d0, 0.0d0, 1.0d0, &
                                                label_(l)//" face. [bc_u, bc_p, u, v, w, p]"
        end do
        write(unit,"(A)") "!**Boundary Condition**"
        write(unit,"(A,1x,i0)") "! Fix value            ", bc_fix_val
        write(unit,"(A,1x,i0)") "! Fix gradient         ", bc_fix_grad
        write(unit,"(A,1x,i0)") "! Periodic             ", bc_periodic
        write(unit,"(A,1x,i0)") "! Outlet unsteady      ", bc_outlet_unsteady
        write(unit,"(A,1x,i0)") "! Periodic with buffer ", bc_periodic_buffer
        write(unit,"(A)") "!bc velocity(u,v,or w) interplated as a pair cell index. "

    close(unit)
end subroutine
    
end module setting_parameter_m