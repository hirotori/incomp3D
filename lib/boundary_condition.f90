module boundary_condition_m
    ! use,intrinsic :: iso_fortran_env, only : IP => int32, DP => real64
    use floating_point_parameter_m
    implicit none

    !!@TODO 境界条件は速度/圧力それぞれで固定値or固定勾配or周期境界を指定する.

    integer(ip),parameter :: bc_fix_val = 1
    integer(ip),parameter :: bc_fix_grad = 2
    integer(ip),parameter :: bc_periodic = 3

    integer(ip),parameter :: bc_outlet_unsteady = 4
    integer(ip),parameter :: bc_periodic_buffer = 5
    
    type bc_t
        !!version : 1.2
        !!境界条件構造体. 1.2から境界条件は固定勾配か固定値かで分岐させる.
        integer(ip) :: type(2) = -99
            !!境界条件の種類. type(1) for velocity, type(2) for pressure.
        real(dp) v_bc(3)
            !!境界面の速度.
        real(dp) p_bc
            !!境界面の圧力.
        integer(ip) :: index_ = -999
            !!仮想セルのインデックス. 
        integer(ip) :: neighbor_index_(2) = -999
            !!index_と対応するペアのセルのインデックス.
            !!圧力と速度とで違う場合がある.
            !!index_(1) for velocity, index_(2) for pressure.
    end type

    public :: bc_t, &
              bc_fix_grad, bc_fix_val, bc_periodic, bc_outlet_unsteady, bc_periodic_buffer, &
              boundary_condition_pressure, &
              boundary_condition_velocity
    
contains
subroutine set_bc_type(types, type_ids, properties, extents)
    type(bc_t),intent(out) :: types(6)
    integer(ip),intent(in) :: type_ids(2,6)
    real(dp),intent(in) :: properties(3+1,6)
    integer(ip),intent(in) :: extents(3)

    integer(ip) l

    print "(A)", ">>> Setting Boundary conditions ..."

    do l = 1, 6
        types(l)%type(:) = type_ids(:,l)
        types(l)%v_bc(:) = properties(1:3,l)
        types(l)%p_bc    = properties(4,l)
    end do

    !i=1 or j=1 or k=1
    do l = 1, 5, 2
        types(l)%index_ = 1
        types(l)%neighbor_index_(:) = merge(extents((l+1)/2), 2, type_ids(:,l) == bc_periodic) !周期境界ならimx(jmx,kmx), それ以外2
    end do

    !i=imax+1 or j=jmax+1 or k=kmax+1
    do l = 2, 6, 2
        types(l)%index_ = extents(l/2) + 1
        types(l)%neighbor_index_(:) = merge(2, extents(l/2), type_ids(:,l) == bc_periodic) !周期境界なら2, それ以外imx+1(jmx+1,kmx+1)
    end do

    !例外. for bc_periodic_buffer
    !この境界条件は, 境界セル以外の内部セルの**速度**を参照する. 
    !参照するインデックスはint(v_bc())から読み取る.
    !この条件は速度にのみ有効.
    block
        integer pair_ids(3), id, pos_
        do l = 1, 6
            if ( types(l)%type(2) == bc_periodic_buffer ) error stop "bc_periodic_buffer is valid for only velocity." 
            if ( types(l)%type(1) == bc_periodic_buffer ) then
                pair_ids(:) = int(types(l)%v_bc, ip) ![u,v,w]を[i,j,k]と対応させる.

                if ( mod(l,2) /= 0 ) then !i=1 or j=1 or k=1
                    pos_ = (l+1)/2
                else                      !i=imax+1 or j=jmax+1 or k=kmax+1
                    pos_ = l/2
                end if

                id = pair_ids(pos_)
                if ( id < 0 .or. extents(pos_) < id) then
                    print "(2(i0,1x,:,'<-->'))", id, extents(pos_)
                    error stop "id is out of range."
                end if
                types(l)%neighbor_index_(1) = pair_ids(pos_)

            end if 
        end do
    end block

    print "('======= Boundary condition Setup =======')"

    block
        character(1) :: label_(6) = ["i","i","j","j","k","k"]
        character(:),allocatable :: format_
        format_ = "(A,' = ',i0, ' :: type (vel,pres) = (', 2(i0,:,','), '), neighbor (vel,pres)= (', 2(i0,:,','), ')')"
        do l = 1, 6
            print format_, label_(l) ,types(l)%index_, types(l)%type(1:2), types(l)%neighbor_index_(1:2)
        end do
    end block

end subroutine


subroutine boundary_condition_velocity(extents, v, bc_type)
    !!境界条件の適用. ここでは仮想セルに値を入れる.
    !!@note 勾配規定条件は全て勾配ゼロ条件のみとする.
    integer(IP),intent(in) :: extents(3)
        !!格子点の個数. 1:x, 2:y, 3:z
    real(DP),intent(inout) :: v(:,:,:,:)
        !!速度
    ! real(dp),intent(in) :: dx(3)
    !     !!格子幅.
    type(bc_t),intent(in) :: bc_type(6)
        !!各面に設定された境界条件. これに対応して仮想セルの値を設定する.

    integer(IP) imx, jmx, kmx
    integer(IP) inb, jnb, knb
    integer(IP) i, j, k, l, m

    imx = extents(1)
    jmx = extents(2)
    kmx = extents(3)
    
    !i-dir.
    do m = 1, 2
        i = bc_type(m)%index_
        inb = bc_type(m)%neighbor_index_(1)

        select case (bc_type(m)%type(1))
        case (bc_fix_val)
            
            do l = 1, 3
                v(l,i,2:jmx,2:kmx) = 2.0_dp*bc_type(m)%v_bc(l) - v(l,inb,2:jmx,2:kmx)
            end do   
    
        !勾配0の条件しか使われないので廃止.  
        ! case (bc_fix_grad)
            
            ! !!勾配規定条件は速度境界条件をそのまま利用する. ただし面に垂直な勾配として与えられる.
            ! !!i方向であれば有効なのはv_bc(1).
            ! do l = 1, 3
            !     v(:,i,2:jmx,2:kmx) = dx(1)*bc_type(m)%v_bc(l) + v(:,inb,2:jmx,2:kmx)              
            ! end do
            
        case (bc_periodic, bc_periodic_buffer, bc_fix_grad)
            
            v(:,i,2:jmx,2:kmx) = v(:,inb,2:jmx,2:kmx)

        case (bc_outlet_unsteady)
            error stop "Sorry, Not Implemented."

        case default
            print "('bc(u,i) = ', i0, ', m = ', i0)", bc_type(m)%type(1), m
            error stop "Invalid BC type."
        end select    

    end do
    
    !j-dir.
    do m = 3, 4
        j = bc_type(m)%index_
        jnb = bc_type(m)%neighbor_index_(1)

        select case (bc_type(m)%type(1))
        case (bc_fix_val)
            
            do l = 1, 3
                v(l,2:imx,j,2:kmx) = 2.0_dp*bc_type(m)%v_bc(l) - v(l,2:imx,jnb,2:kmx)
            end do   
    
        ! case (bc_fix_grad)
            
        !     do l = 1, 3
        !         v(l,2:imx,j,2:kmx) = dx(2)*bc_type(m)%v_bc(l) + v(l,2:imx,jnb,2:kmx)
        !     end do

        case (bc_periodic, bc_periodic_buffer, bc_fix_grad)
            
            v(:,2:imx,j,2:kmx) = v(:,2:imx,jnb,2:kmx)
    
        case (bc_outlet_unsteady)
            error stop "Sorry, Not Implemented."
        case default
            print "('bc(u,j) = ', i0, ', m = ', i0)", bc_type(m)%type(1), m
            error stop "Invalid BC type."
        end select    

    end do
    !k-dir.
    do m = 5, 6
        k = bc_type(m)%index_
        knb = bc_type(m)%neighbor_index_(1)

        select case (bc_type(m)%type(1))
        case (bc_fix_val)
            
            do l = 1, 3
                v(l,2:imx,2:jmx,k) = 2.0_dp*bc_type(m)%v_bc(l) - v(l,2:imx,2:jmx,knb)
            end do   
    
        ! case (bc_fix_grad)
            
        !     do l = 1, 3
        !         v(l,2:imx,2:jmx,k) = dx(3)*bc_type(m)%v_bc(l) + v(l,2:imx,2:jmx,knb)
        !     end do

        case (bc_periodic, bc_periodic_buffer, bc_fix_grad)
            
            v(:,2:imx,2:jmx,k) = v(:,2:imx,2:jmx,knb)
    
        case (bc_outlet_unsteady)
            error stop "Sorry, Not Implemented."
        case default
            print "('bc(u,j) = ', i0, ', m = ', i0)", bc_type(m)%type(1), m
            error stop "Invalid BC type."
        end select    

    end do

end subroutine

subroutine boundary_condition_pressure(extents, p, bc_type)
    !!境界条件の適用. ここでは仮想セルに値を入れる.
    !!@note 勾配規定条件は全て勾配ゼロ条件のみとする.
    integer(IP),intent(in) :: extents(3)
    real(DP),intent(inout) :: p(:,:,:)
    ! real(dp),intent(in) :: dx(3)
    type(bc_t),intent(in) :: bc_type(6)

    integer(IP) inb, jnb, knb
    integer(IP) i, j, k, m
    integer(IP) imx, jmx, kmx

    imx = extents(1)
    jmx = extents(2)
    kmx = extents(3)

    !i-dir.
    do m = 1, 2
        i = bc_type(m)%index_
        inb = bc_type(m)%neighbor_index_(2)

        select case (bc_type(m)%type(2))
        case (bc_fix_val)

            p(i,2:jmx,2:kmx) = 2.0_dp*bc_type(m)%p_bc - p(inb,2:jmx,2:kmx)
           
        ! case (bc_fix_grad)

        !     !!固定勾配の場合は境界の圧力値を面に垂直な勾配として扱う.
        !     p(i,2:jmx,2:kmx) = dx(1)*bc_type(m)%p_bc + p(inb,2:jmx,2:kmx)
        
        case (bc_periodic, bc_fix_grad)
            p(i,2:jmx,2:kmx) = p(inb,2:jmx,2:kmx)
            
        case (bc_periodic_buffer)
            error stop "bc_periodic_buffer is only valid for velocity. Check your config file."

        case (bc_outlet_unsteady)
            error stop "Sorry, Not Implemented."

        case default
            print "('bc(p,i) = ', i0)", bc_type(m)%type(2)
            error stop "Invalid BC type."
        end select    

    end do
    
    !j-dir.
    do m = 3, 4
        j = bc_type(m)%index_
        jnb = bc_type(m)%neighbor_index_(2)

        select case (bc_type(m)%type(2))
        case (bc_fix_val)

            p(2:imx,j,2:kmx) = 2.0_dp*bc_type(m)%p_bc - p(2:imx,jnb,2:kmx)
           
        ! case (bc_fix_grad)

        !     p(2:imx,j,2:kmx) = dx(2)*bc_type(m)%p_bc + p(2:imx,jnb,2:kmx)

        case (bc_periodic, bc_fix_grad)
            p(2:imx,j,2:kmx) = p(2:imx,jnb,2:kmx)

        case (bc_periodic_buffer)
            error stop "bc_periodic_buffer is only valid for velocity. Check your config file."

        case (bc_outlet_unsteady)
            error stop "Sorry, Not Implemented."

        case default
            print "('bc(p,j) = ', i0)", bc_type(m)%type(2)
            error stop "Invalid BC type."
        end select

    end do
    !k-dir.
    do m = 5, 6
        k = bc_type(m)%index_
        knb = bc_type(m)%neighbor_index_(2)

        select case (bc_type(m)%type(2))
        case (bc_fix_val)

            p(2:imx,2:jmx,k) = 2.0_dp*bc_type(m)%p_bc - p(2:imx,2:jmx,knb)
           
        ! case (bc_fix_grad)

        !     p(2:imx,2:jmx,k) = dx(3)*bc_type(m)%p_bc + p(2:imx,2:jmx,knb)

        case (bc_periodic, bc_fix_grad)
            
            p(2:imx,2:jmx,k) = p(2:imx,2:jmx,knb)

        case (bc_periodic_buffer)
            error stop "bc_periodic_buffer is only valid for velocity. Check your config file."

        case (bc_outlet_unsteady)
            error stop "Sorry, Not Implemented."

        case default
            print "('bc(p,k) = ', i0)", bc_type(m)%type(2)
            error stop "Invalid BC type."
        end select

    end do
    !fix pressure
    ! p(1,2,2) = 0.0_dp

end subroutine
end module boundary_condition_m