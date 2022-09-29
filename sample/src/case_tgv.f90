module case_tgv_m
    use floating_point_parameter_m
    use case_common_m
    use mesh_m
    use fluid_field_m
    implicit none

    type, extends(case_common_t), public :: case_tgv_t
        real(dp),allocatable :: exact(:,:)    
    contains
        procedure add_on_pre_process
        procedure phase_post_process
    
    end type 
    
contains
subroutine add_on_pre_process(this, grid, fld)
    !!追加で初期状態を管理したい場合に呼び出す.
    class(case_tgv_t),intent(inout) :: this
    type(rectilinear_mesh_t),intent(inout) :: grid
    type(fluid_field_t),intent(inout) :: fld

    !!初期条件を設定する.
    real(dp) Re
    integer(ip) i, j, k, imx, jmx, kmx
    real(dp) rc_(3)

    call grid%get_extents_sub(imx, jmx, kmx)
    Re = this%settings_case%reynolds_number
    do k = 2, kmx
    do j = 2, jmx
    do i = 2, imx
        rc_(:) = grid%rc(:,i,j,k)
        fld%velocity(1,i,j,k) = -cos(rc_(1))*sin(rc_(2))
        fld%velocity(2,i,j,k) =  sin(rc_(2))*cos(rc_(1))
        fld%pressure(i,j,k) = -0.25_dp*Re*(cos(2.0_dp*rc_(1)) + cos(2.0_dp*rc_(2)))
    end do
    end do        
    end do

end subroutine

subroutine phase_post_process(this, grid, fld)
    !!現時間段階の計算が終わった後に呼び出される. 
    class(case_tgv_t),intent(inout) :: this
    type(rectilinear_mesh_t),intent(in) :: grid
    type(fluid_field_t),intent(in) :: fld

    !HELLO

end subroutine
    
end module case_tgv_m