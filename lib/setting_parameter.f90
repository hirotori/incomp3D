module setting_parameter_m
    ! use,intrinsic :: iso_fortran_env, only : IP => int32, DP => real64
    use floating_point_parameter_m
    use boundary_condition_m, only : bc_inlet, bc_outlet, bc_wall, bc_slip, &
                                     bc_periodic, bc_outlet_unsteady, bc_periodic_buffer
    implicit none
    private
    type slv_setting
        !!poisson方程式のソルバ設定.
        real(DP) tolerance
        !!反復終了のための残差の閾値.
        integer(IP) itr_max
        !!反復回数上限
        real(DP) relax_coeff
        !!加速係数.
        character(2) :: conv_type
        integer(ip) :: diff_type
        logical :: correct_face_flux
    end type

    type case_setting
        integer(ip) nstart, nend, nwrite
        real(dp) dt
        real(dp) p_ref
        real(dp) u_ic(3)
        real(dp) body_force(3)
        real(dp) reynolds_number
    end type
    
    public slv_setting, case_setting, read_config

contains
subroutine read_config(fname, imx, jmx, kmx, l_xyz, c_setting, s_setting, bc_ids, bc_properties)
    character(*),intent(in) :: fname
    integer(IP),intent(out) :: imx, jmx, kmx
    real(DP),intent(out) :: l_xyz(3)
    type(case_setting),intent(out) :: c_setting
    type(slv_setting),intent(out) :: s_setting
    integer(ip),intent(out) :: bc_ids(6)
    real(dp),intent(out) :: bc_properties(4,6)



    call read_config_core(fname, imx, jmx, kmx, l_xyz, &
                          c_setting, s_setting, &
                          bc_ids, bc_properties)

    !TODO リファクタリング. namelistにして各構造体は別ファイル扱いとする??

end subroutine

subroutine read_config_core(fname, imx, jmx, kmx, l_xyz, setting_c ,setting, bc_ids, bc_properties)
    character(*),intent(in) :: fname
    integer(IP),intent(out) :: imx, jmx, kmx
    real(DP),intent(out) :: l_xyz(3)
    type(case_setting),intent(out) :: setting_c
    type(slv_setting),intent(out) :: setting
    integer(ip),intent(out) :: bc_ids(6)
    real(dp),intent(out) :: bc_properties(4,6)

    integer(IP) unit, iost, l
    character(128) tmp_

    open(newunit = unit, file = fname, action = "read", status = "old", iostat = iost, iomsg = tmp_)
    if ( iost /= 0 ) then
        call create_sample_()
        error stop trim(tmp_)
    end if

    read(unit,*,err=99) imx, jmx, kmx 
    read(unit,*,err=99) l_xyz(:)      
    read(unit,*,err=99) setting_c%nstart, setting_c%nend, setting_c%nwrite
    read(unit,*,err=99) setting_c%dt
    read(unit,*,err=99) setting_c%reynolds_number
    read(unit,*,err=99) setting_c%body_force(:)
    read(unit,*,err=99) setting_c%p_ref
    read(unit,*,err=99) setting_c%u_ic(:)
    read(unit,*,err=99) setting%itr_max
    read(unit,*,err=99) setting%tolerance
    read(unit,*,err=99) setting%relax_coeff
    do l = 1, 6
        read(unit,*,err=99) bc_ids(l), bc_properties(:,l)        
    end do
    read(unit,*,err=99) tmp_ ; setting%conv_type = trim(adjustl(tmp_))
    read(unit,*,err=99) setting%diff_type
    read(unit,*,err=99) setting%correct_face_flux
    print "('extent : (', i0, ', ', i0, ', ', i0, ') ')", imx, jmx, kmx
    print "('Length : (', g0.2, ', ', g0.2, ', ', g0.2, ') ')", l_xyz(:)
    print "('time step : ',i0, '<= n <=', i0)", setting_c%nstart, setting_c%nend
    print "('delta t   : ',g0.5)", setting_c%dt
    print "('Re        : ',g0)", setting_c%reynolds_number
    print "('Force  : (', g0.2, ', ', g0.2, ', ', g0.2, ') ')", setting_c%body_force(:)
    print "('Poisson max iteration = ',i0)", setting%itr_max
    print "('Poisson tolerance     = ',g0.6)", setting%tolerance
    print "('Poisson alpha         = ',g0.2)", setting%relax_coeff 
    close(unit)

    return

99 continue

    print "('ERROR OCCURED IN READ_CONFIG')"

    error stop

end subroutine

subroutine create_sample_()
    integer(ip) unit
    integer(ip) l
    character(8) :: label_(6) = ["i = 1   ", "i = imax", "j = 1   ", "j = jmax", "k = 1   ", "k = kmax"] 
    character(*),parameter :: dummy_flag = "T"
    open(newunit = unit, file = "config_sample.txt", status = "replace")
        write(unit, "(3(i0,1x), ' !grid points (imax,jmax,kmax)')") 21, 21, 21
        write(unit, "(3(g0.3,1x), ' !length (x,y,z)     ')") 10.0d0, 4.0d0, 2.0d0 
        write(unit, "(3(i0,1x), ' !time step (start,end,write_interval')") 0, 1000, 100 
        write(unit, "(g0.3      , ' !time spacing')") 0.001d0         
        write(unit, "(g0.2      , ' !Reynolds Number')") 100.0d0
        write(unit, "(3(g0.2,1x), ' !body force (x,y,z) ')") 0.0d0, 0.0d0, 0.0d0             
        write(unit, "(g0.2      , ' !Pressure reference (also for initial condition)')") 1.0d0
        write(unit, "(3(g0.2,1x), ' !initial condition for velocity (x,y,z)     ')") 0.0d0, 0.0d0, 0.0d0 
        write(unit, "(i0      , ' !max iteration no. for poisson eq.')") 1000
        write(unit, "(g0.2      , ' !tolerance for poisson eq.')") 0.0001d0
        write(unit, "(g0.2      , ' !accel. coefficient for poisson eq.')") 1.3d0
        do l = 1, 6
            write(unit, "(i0,1x,4(g0.2,1x),' !B.C. for ',A)") 1, 0.0d0, 0.0d0, 0.0d0, 1.0d0, label_(l)//" face. [bc_id, u, v, w, p]"
        end do
        write(unit, "(A         , ' !convection term. ""ud"": 1st order upwind, ""cd"": central diff')") "ud"
        write(unit, "(i0        , ' !diffusion  term.   1 : compact stencil ,   2 : large stencil')") 1
        write(unit, "(A         , ' ![T/F] if T, face fluxes are corrected by pressure.')") dummy_flag
        write(unit,"(A)") "!**Boundary Condition**"
        write(unit,"(A,1x,i0)") "! Inlet          ", bc_inlet
        write(unit,"(A,1x,i0)") "! Outlet         ", bc_outlet
        write(unit,"(A,1x,i0)") "! Wall           ", bc_wall
        write(unit,"(A,1x,i0)") "! Slip           ", bc_slip
        write(unit,"(A,1x,i0)") "! Periodic       ", bc_periodic
        write(unit,"(A,1x,i0)") "! Outlet unsteady ", bc_outlet_unsteady
        write(unit,"(A,1x,i0)") "! Periodic with buffer ", bc_periodic_buffer
        write(unit,"(A)") "!bc velocity(u,v,or w) interplated as a pair cell index. "

    close(unit)
end subroutine
    
end module setting_parameter_m