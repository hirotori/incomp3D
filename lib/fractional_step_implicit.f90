module fractional_step_implicit_m
    use floating_point_parameter_m
    use fluid_field_m
    use mesh_m
    use fractional_step_m, only : solver_fs, mat_a, slvr_init_common, calc_convective_and_diffusive_flux
    use setting_parameter_m, only : slv_setting, case_setting
    use boundary_condition_m, only : bc_t, boundary_condition_velocity
    implicit none
    
    private

    type,extends(solver_fs) :: solver_fs_imp_t
        !!陰解法Fractional stepクラス. 粘性項のみをクランク･ニコルソン法で評価する半陰解法.
        !!@warning 現状, diffus_type = 2には対応出来ていない. diffus_type = 2を選ぶと恐らく破綻する(整合性がないので).
        type(mat_a) :: matrix_v
        real(dp),allocatable :: rhs_(:,:,:,:)
        integer(ip) :: iter_max = 1000
        real(dp) :: tolerance = 0.001_dp
        real(dp) :: alpha = 1.0_dp
        contains
        procedure :: init => init_solver_2
        procedure predict_pseudo_velocity

    end type

    public solver_fs_imp_t, calc_pseudo_velocity_common_core
contains
subroutine init_solver_2(this, fld, grd, settings_slv, setting_case)
    !!ソルバを初期化する.
    class(solver_fs_imp_t),intent(inout) :: this
    type(fluid_field_t),intent(in) :: fld
    type(rectilinear_mesh_t),intent(in) :: grd
    type(slv_setting),intent(in) :: settings_slv
    type(case_setting),intent(in) :: setting_case

    integer(ip) mx(3)

    mx(:) = grd%get_extents()

    call slvr_init_common(this, fld, grd, settings_slv, setting_case)
    
    call set_matrix_(this%matrix_v, grd%ds, grd%dv, grd%dx, setting_case%reynolds_number, setting_case%dt, this%diffus_type)
    allocate(this%rhs_(1:3,2:mx(1),2:mx(2),2:mx(3)))

end subroutine

subroutine predict_pseudo_velocity(this, extents, ds, dv, dx, del_t, re, force, v0, v, mi, mj, mk, dudr, bc_types)
    !!中間速度を計算する.
    class(solver_fs_imp_t),intent(inout) :: this
    integer(ip),intent(in) :: extents(3)
    real(DP),intent(in) :: ds(3)
        !!検査面の大きさ
    real(DP),intent(in) :: dv
        !!検査体積の大きさ
    real(DP),intent(in) :: dx(3)
        !!検査体積の幅
    real(DP),intent(in) :: del_t
        !!時間刻み
    real(DP),intent(in) :: re
        !!レイノルズ数
    real(DP),intent(in) :: force(3)
        !!体積力
    real(DP),intent(in) :: v0(:,:,:,:)
        !!既知時間段階の速度
    real(DP),intent(out) :: v(:,:,:,:)
        !!更新される速度(中間段階)
    real(DP),intent(inout) :: mi(:,2:,2:), mj(2:,:,2:), mk(2:,2:,:)
        !!面に垂直な流速成分.
    real(dp),intent(in) :: dudr(:,:,:,:,:)
        !!速度勾配テンソル.
    type(bc_t),intent(in) :: bc_types(:)
    
    real(dp),allocatable :: conv(:,:,:,:), diff(:,:,:,:)
    integer(ip) i, j, k, l
    real(dp) rei

    allocate(conv(3,2:extents(1),2:extents(2),2:extents(3)))
    allocate(diff(3,2:extents(1),2:extents(2),2:extents(3)))    
        
    call calc_convective_and_diffusive_flux(this, extents, ds, dv, dx, re, v0, mi, mj, mk, dudr, conv, diff)

    rei = 1.0_dp/re
    do k = 2, extents(3)
    do j = 2, extents(2)
    do i = 2, extents(1)
        this%rhs_(:,i,j,k) = v0(:,i,j,k)*dv - del_t*(conv(:,i,j,k) - 0.5_dp*rei*diff(:,i,j,k) - force(:)*dv)
    end do
    end do        
    end do

    call calc_pseudo_velocity_common_core(this, extents, v, bc_types)

    !@TODO 発散した後の処理. 


end subroutine

subroutine calc_pseudo_velocity_common_core(this, extents, v, bc_types)
    !!時間陰解法に共通の処理. 連立方程式を解いて新しい中間速度を求める.
    !!係数行列と右辺は外部で計算されている前提.
    use,intrinsic :: ieee_arithmetic, only : ieee_is_nan
    class(solver_fs_imp_t),intent(inout) :: this
    integer(ip),intent(in) :: extents(3)
    real(dp),intent(inout) :: v(:,:,:,:)
    type(bc_t),intent(in) :: bc_types(:)

    integer(ip) i, j, k, l
    integer(ip) iter
    real(dp) den_, resid_
    real(dp),allocatable :: v_tmp(:,:,:,:)
    character(3) :: label_(3) = ["[u]", "[v]", "[w]"]
    integer(ip) imx, jmx, kmx, icmx

    imx = extents(1)
    jmx = extents(2)
    kmx = extents(3)
    icmx = product(extents(:) - 1)
    allocate(v_tmp, mold = v)
    v_tmp(:,:,:,:) = 0.0_dp

    !各速度成分毎に方程式を解く. 
    !3成分まとめて解いても良いのか(l=1, 3に関するループではなく, v(1:3,i,j,k)とするようなループ)?
    do l = 1, 3
        print "(A)", '------- velocity '//label_(l)//' loop -------'
    do iter = 1, this%iter_max
        do k = 2, kmx
        do j = 2, jmx
        do i = 2, imx
            !@TODO large stencilの場合の対応. ステンシルが一つ飛びになる.
            v(l,i,j,k) = (1.0_dp - this%alpha)*v(l,i,j,k) &
                                 + this%alpha*(this%rhs_(l,i,j,k) &
                                 - this%matrix_v%aw*v(l,i-1,j  ,k  ) - this%matrix_v%ae*v(l,i+1,j  ,k  ) &
                                 - this%matrix_v%as*v(l,i  ,j-1,k  ) - this%matrix_v%an*v(l,i  ,j+1,k  ) &
                                 - this%matrix_v%ab*v(l,i  ,j  ,k-1) - this%matrix_v%at*v(l,i  ,j  ,k+1))/this%matrix_v%ap
        end do
        end do        
        end do
    
        call boundary_condition_velocity(extents, v, bc_types)

        resid_ = 0.0_dp
        den_ = 0.0_dp
        do k = 2, kmx
        do j = 2, jmx
        do i = 2, imx
            resid_ = resid_ + (v(l,i,j,k) - v_tmp(l,i,j,k))**2.0_dp
            den_ = den_ + (v(l,i,j,k))**2.0_dp
        end do
        end do
        end do

        if ( abs(den_) == 0.0_dp ) exit !流れによっては速度成分がずっとゼロの場合もありうる. NANを防ぐ.

        resid_ = sqrt(resid_/icmx)/sqrt(den_/icmx)

        if ( resid_ <= this%tolerance .or. iter >= this%iter_max) then
            print "('velocity',A,'(',i0,') converged, resid = ',g0, ' den = ',g0)", label_(l), iter, resid_, den_
            ! print "('exit.')"
            exit
        end if
        
        if ( resid_ > huge(0.0_dp)) then
            print "(A)", "velocity "//label_(l)//" diverged."
            ! diverged = .true.
            return
        end if

        if ( any(ieee_is_nan(v(l,:,:,:)))) then
            print"(A)", "NAN detected in velocity "//label_(l)
            ! diverged = .true.
            return
        end if

        if ( (iter > 1000 .and. mod(iter,1000) == 0) ) then
            print "('velocity',A,'(',i0,'), resid = ',g0, ' den = ',g0)", label_(l), iter, resid_, den_
        end if
        
        v_tmp(l,:,:,:) = v(l,:,:,:)

    end do
    end do



end subroutine

subroutine set_matrix_(mat, ds, dv, dx, re, dt, stencil_type)
    !!係数行列を設定する.
    type(mat_a),intent(inout) :: mat
    real(DP),intent(in) :: ds(3)
        !!検査面の大きさ
    real(DP),intent(in) :: dv
        !!検査体積の大きさ
    real(DP),intent(in) :: dx(3)
        !!検査体積の幅
    real(DP),intent(in) :: re
        !!レイノルズ数
    real(DP),intent(in) :: dt
        !!時間刻み
    integer(ip),intent(in) :: stencil_type
        !!large stencilかそうで無いかを区別する.

    !@TODO large stencilの場合係数が微妙に変わる.それの実装.
    mat%ae = -0.5_dp*dt*ds(1)/(re*dx(1))
    mat%aw = mat%ae
    mat%an = -0.5_dp*dt*ds(2)/(re*dx(2))
    mat%as = mat%an
    mat%at = -0.5_dp*dt*ds(3)/(re*dx(3))
    mat%ab = mat%at

    mat%ap = dv - (mat%ae + mat%aw + mat%an + mat%as + mat%at + mat%ab)


end subroutine

! subroutine sor(extents, A, x, b, bc_types, tolerance, max_iter, accl, diverged)
!     use,intrinsic :: ieee_arithmetic, only : ieee_is_nan
!     !!SOR法. 速度成分をそれぞれループにかけて解く.
!     integer(ip),intent(in) :: extents(3)
!     type(mat_a),intent(in) :: A
!     real(dp),intent(inout) :: x(:,:,:)
!     real(dp),intent(in) :: b(:,:,:)
!     type(bc_t),intent(in) :: bc_types(:)
!     real(dp),intent(in) :: tolerance
!     integer(ip),intent(in) :: max_iter
!     real(dp),intent(in) :: accl
!     logical,intent(out) :: diverged

!     integer(ip) iter, i, j, k, imx, jmx, kmx, icmx
!     real(dp) den_, resid_
!     real(dp),allocatable :: x_tmp(:,:,:,:)
    
!     imx = extents(1)
!     jmx = extents(2)
!     kmx = extents(3)
!     icmx = product(extents(:) - 1)

!     allocate(x_tmp, mold = x)
!     x_tmp(:,:,:,:) = 0.0_dp

!     do iter = 1, max_iter
!         do k = 2, kmx
!         do j = 2, jmx
!         do i = 2, imx
!             x(i,j,k) = (1.0_dp - accl)*x(i,j,k) &
!                                  + accl*(b(i,j,k) &
!                                  - A%aw*x(i-1,j  ,k  ) - A%ae*x(i+1,j  ,k  ) &
!                                  - A%as*x(i  ,j-1,k  ) - A%an*x(i  ,j+1,k  ) &
!                                  - A%ab*x(i  ,j  ,k-1) - A%at*x(i  ,j  ,k+1))/A%ap
!         end do
!         end do        
!         end do
    
!         call boundary_condition_velocity(extents, v, bc_types)

!         resid_ = 0.0_dp
!         den_ = 0.0_dp
!         do k = 2, kmx
!         do j = 2, jmx
!         do i = 2, imx
!             resid_ = resid_ + norm2(v(:,i,j,k) - v_tmp(:,i,j,k))**2.0_dp
!             den_ = den_ + norm2(v(:,i,j,k))**2.0_dp
!         end do
!         end do
!         end do

!         resid_ = sqrt(resid_/icmx)/sqrt(den_/icmx)

!         if ( resid_ <= this%tolerance .or. iter >= this%iter_max) then
!             print "('velocity(',i0,') converged, resid = ',g0)", iter, resid_
!             ! print "('exit.')"
!             exit
!         end if
        
!         if ( resid_ > huge(0.0_dp)) then
!             print "(A)", "velocity diverged."
!             ! diverged = .true.
!             return
!         end if

!         if ( any(ieee_is_nan(v))) then
!             print"(A)", "NAN detected in pressure."
!             ! diverged = .true.
!             return
!         end if

!         if ( (iter > 1000 .and. mod(iter,1000) == 0) ) then
!             print "('velocity(',i0,'), resid = ',g0)", iter, resid_
!         end if
        
!         v_tmp(:,:,:,:) = v(:,:,:,:)

!     end do

! end subroutine


end module fractional_step_implicit_m
