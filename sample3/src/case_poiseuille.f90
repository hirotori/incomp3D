module case_poiseuille_m
    use floating_point_parameter_m
    use case_common_m
    use mesh_m
    use fluid_field_m
    use IO_operator_m
    implicit none

    type, extends(case_common_t),public :: case_poiseuille_t
        real(dp),allocatable :: exact(:,:,:,:)        
            !!厳密解.
        
        contains
        procedure add_on_pre_process
        procedure phase_post_process
    
    end type 
    
contains

subroutine add_on_pre_process(this, grid, fld)
    !!追加で初期状態を管理したい場合に呼び出す.
    class(case_poiseuille_t),intent(inout) :: this
    class(equil_mesh_t),intent(inout) :: grid
    type(fluid_field_t),intent(inout) :: fld

    integer(ip) i, j, k, imx, jmx, kmx
    real(dp) y_

    call grid%get_extents_sub(imx, jmx, kmx)

    allocate(this%exact(1:3,2:imx,2:jmx,2:kmx),source = 0.0_dp)

    !! 初期条件としてポアズイユ流の厳密解を与える. ここでは比較用に厳密解を配列に持たせておく.
    do k = 2, kmx
    do j = 2, jmx
    do i = 2, imx
        y_ = grid%yc(j)
        this%exact(1,i,j,k) = 1.5_dp*(y_*(2.0_dp - y_))
    end do
    end do    
    end do

    fld%velocity(:,2:imx,2:jmx,2:kmx) = this%exact(:,2:imx,2:jmx,2:kmx)


end subroutine

subroutine phase_post_process(this, grid, fld, step, time)
    !!現時間段階の計算が終わった後に呼び出される. 
    class(case_poiseuille_t),intent(inout) :: this
    class(equil_mesh_t),intent(in) :: grid
    type(fluid_field_t),intent(in) :: fld
    integer(ip),intent(in) :: step
    real(dp),intent(in) :: time

    real(dp) resid_, den_
    integer(ip) i, j, k, imx, jmx, kmx
    integer(ip) unit

    !!厳密解との誤差ノルムの計算.
    call grid%get_extents_sub(imx, jmx, kmx)

    resid_ = 0.0_dp
    den_ = 0.0_dp

    do k = 2, kmx
    do j = 2, jmx
    do i = 2, imx
        resid_ = resid_ + norm2(fld%velocity(1:3,i,j,k) - this%exact(1:3,i,j,k))**2.0_dp
        den_ = den_ + norm2(this%exact(1:3,i,j,k))**2.0_dp
    end do
    end do        
    end do

    resid_ = sqrt(resid_/den_)

    print "('L2-norm (velocity - exact) = ',g0)", resid_

    !!l2-normはファイルにストレージされる.
    if ( open_text_file(unit, "result.txt", "# time, l2-err") ) then
        write(unit,"(*(g0,1x))") time, resid_
    end if

end subroutine
    
end module case_poiseuille_m