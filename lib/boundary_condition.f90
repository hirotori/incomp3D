module boundary_condition_m
    ! use,intrinsic :: iso_fortran_env, only : IP => int32, DP => real64
    use floating_point_parameter_m
    implicit none

    integer(ip),parameter :: bc_inlet = 1
    integer(ip),parameter :: bc_outlet = 2
    integer(ip),parameter :: bc_wall = 3
    integer(ip),parameter :: bc_slip = 4
    integer(ip),parameter :: bc_periodic = 5
    integer(ip),parameter :: bc_outlet_unsteady = 6
    integer(ip),parameter :: bc_periodic_buffer = 7
    
    type bc_t
        integer(ip) type
        real(dp) v_bc(3)
        real(dp) p_bc
        integer(ip) index
        integer(ip) neighbor_index
    end type

    public :: bc_t, &
              bc_inlet, bc_outlet, bc_wall, bc_slip, bc_periodic, bc_outlet_unsteady, bc_periodic_buffer, &
              boundary_condition_pressure, &
              boundary_condition_velocity
    
contains
subroutine set_bc_type(types, type_ids, properties, extents)
    type(bc_t),intent(out) :: types(6)
    integer(ip),intent(in) :: type_ids(6)
    real(dp),intent(in) :: properties(3+1,6)
    integer(ip),intent(in) :: extents(3)

    integer(ip) l

    do l = 1, 6
        types(l)%type = type_ids(l)
        types(l)%v_bc(:) = properties(1:3,l)
        types(l)%p_bc    = properties(4,l)
    end do

    !i=1 or j=1 or k=1
    do l = 1, 5, 2
        types(l)%index = 1
        types(l)%neighbor_index = merge(extents((l+1)/2), 2, type_ids(l) == bc_periodic) !周期境界ならimx(jmx,kmx), それ以外2
    end do

    !i=imax+1 or j=jmax+1 or k=kmax+1
    do l = 2, 6, 2
        types(l)%index = extents(l/2) + 1
        types(l)%neighbor_index = merge(2, extents(l/2), type_ids(l) == bc_periodic) !周期境界なら2, それ以外imx+1(jmx+1,kmx+1)
    end do

    !例外. for bc_periodic_buffer
    !この境界条件は, 境界セル以外の内部セルの速度を参照する. 
    !参照するインデックスはint(v_bc())から読み取る.
    block
        integer pair_ids(3), id, pos_
        do l = 1, 6
            if ( types(l)%type == bc_periodic_buffer ) then
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
                types(l)%neighbor_index = pair_ids(pos_)

            end if 
        end do
    end block

    print "(' -- Boundary condition Setup -- ')"

    block
        character(1) :: label_(6) = ["i","i","j","j","k","k"]
        do l = 1, 6
            print "(A,' = ',i0, ' :: type = ',i0, ' (neighbor_index = ', i0, ' )')", &
            label_(l) ,types(l)%index, types(l)%type, types(l)%neighbor_index
        end do
    end block

end subroutine


subroutine boundary_condition_velocity(extents, v, bc_type)
    !!境界条件の適用. ここでは仮想セルに値を入れる.
    integer(IP),intent(in) :: extents(3)
        !!格子点の個数. 1:x, 2:y, 3:z
    real(DP),intent(inout) :: v(:,:,:,:)
        !!速度
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
        i = bc_type(m)%index
        inb = bc_type(m)%neighbor_index

        select case (bc_type(m)%type)
        case (bc_inlet, bc_wall)
            do l = 1, 3
                v(l,i,2:jmx,2:kmx) = 2.0_dp*bc_type(m)%v_bc(l) - v(l,inb,2:jmx,2:kmx)
            end do   
    
        case (bc_slip, bc_periodic, bc_outlet, bc_periodic_buffer)

            v(:,i,2:jmx,2:kmx) = v(:,inb,2:jmx,2:kmx)
    
        case (bc_outlet_unsteady)
            error stop "Sorry, Not Implemented."

        case default
            error stop "Invalid BC type."
        end select    

    end do
    
    !j-dir.
    do m = 3, 4
        j = bc_type(m)%index
        jnb = bc_type(m)%neighbor_index

        select case (bc_type(m)%type)
        case (bc_inlet, bc_wall)
            do l = 1, 3
                v(l,2:imx,j,2:kmx) = 2.0_dp*bc_type(m)%v_bc(l) - v(l,2:imx,jnb,2:kmx)
            end do   
    
        case (bc_slip, bc_periodic, bc_outlet, bc_periodic_buffer)

            v(:,2:imx,j,2:kmx) = v(:,2:imx,jnb,2:kmx)
    
        case (bc_outlet_unsteady)
            error stop "Sorry, Not Implemented."

        case default
            error stop "Invalid BC type."
        end select    

    end do
    !k-dir.
    do m = 5, 6
        k = bc_type(m)%index
        knb = bc_type(m)%neighbor_index

        select case (bc_type(m)%type)
        case (bc_inlet, bc_wall)
            do l = 1, 3
                v(l,2:imx,2:jmx,k) = 2.0_dp*bc_type(m)%v_bc(l) - v(l,2:imx,2:jmx,knb)
            end do   
    
        case (bc_slip, bc_periodic, bc_outlet, bc_periodic_buffer)

            v(:,2:imx,2:jmx,k) = v(:,2:imx,2:jmx,knb)
    
        case (bc_outlet_unsteady)
            error stop "Sorry, Not Implemented."

        case default
            error stop "Invalid BC type."
        end select    

    end do

end subroutine

subroutine boundary_condition_pressure(extents, p, bc_type)
    !!境界条件の適用. ここでは仮想セルに値を入れる.
    integer(IP),intent(in) :: extents(3)
    real(DP),intent(inout) :: p(:,:,:)
    type(bc_t),intent(in) :: bc_type(6)

    integer(IP) inb, jnb, knb
    integer(IP) i, j, k, m
    integer(IP) imx, jmx, kmx

    imx = extents(1)
    jmx = extents(2)
    kmx = extents(3)

    !i-dir.
    do m = 1, 2
        i = bc_type(m)%index
        inb = bc_type(m)%neighbor_index

        select case (bc_type(m)%type)
        case (bc_inlet, bc_wall, bc_slip, bc_periodic)
           
            p(i,2:jmx,2:kmx) = p(inb,2:jmx,2:kmx)
           
        case (bc_outlet)

            p(i,2:jmx,2:kmx) = 2.0_dp*bc_type(m)%p_bc - p(inb,2:jmx,2:kmx)

        case (bc_periodic_buffer)
            !これだけNeumann境界を別途実装する必要があった. インデックスで判定出来ないので, 
            !オフセットをm=1なら1, m=2なら-1となるように設定.
            p(i,2:jmx,2:kmx) = p(i-2*m+3,2:jmx,2:kmx)

        case (bc_outlet_unsteady)
            error stop "Sorry, Not Implemented."

        case default
            error stop "Invalid BC type."
        end select    

    end do
    
    !j-dir.
    do m = 3, 4
        j = bc_type(m)%index
        jnb = bc_type(m)%neighbor_index

        select case (bc_type(m)%type)
        case (bc_inlet, bc_wall, bc_slip, bc_periodic)
           
            p(2:imx,j,2:kmx) = p(2:imx,jnb,2:kmx)
               
        case (bc_outlet)

            p(2:imx,j,2:kmx) = 2.0_dp*bc_type(m)%p_bc - p(2:imx,jnb,2:kmx)

        case (bc_periodic_buffer)
            !これだけNeumann境界を別途実装する必要があった. インデックスで判定出来ないので, 
            !オフセットをm=3なら1, m=4なら-1となるように設定.
            p(2:imx,j,2:kmx) = p(2:imx,j-2*(m-2)+3,2:kmx)

        case (bc_outlet_unsteady)
            error stop "Sorry, Not Implemented."

        case default
            error stop "Invalid BC type."
        end select    

    end do
    !k-dir.
    do m = 5, 6
        k = bc_type(m)%index
        knb = bc_type(m)%neighbor_index

        select case (bc_type(m)%type)
        case (bc_inlet, bc_wall, bc_slip, bc_periodic)
           
            p(2:imx,2:jmx,k) = p(2:imx,2:jmx,knb)
               
        case (bc_outlet)

            p(2:imx,2:jmx,k) = 2.0_dp*bc_type(m)%p_bc - p(2:imx,2:jmx,knb)

        case (bc_periodic_buffer)
            !これだけNeumann境界を別途実装する必要があった. インデックスで判定出来ないので, 
            !オフセットをm=5なら1, m=6なら-1となるように設定.
            p(2:imx,2:jmx,k) = p(2:imx,2:jmx,k-2*(m-4)+3)
        
        case (bc_outlet_unsteady)
            error stop "Sorry, Not Implemented."

        case default
            error stop "Invalid BC type."
        end select    

    end do
    !fix pressure
    ! p(1,2,2) = 0.0_dp

end subroutine
end module boundary_condition_m