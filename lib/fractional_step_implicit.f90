module fractional_step_implicit_m
    use floating_point_parameter_m
    use fluid_field_m
    use mesh_m
    use fractional_step_m, only : solver_fs, mat_a, slvr_init_common, calc_convective_and_diffusive_flux
    use setting_parameter_m, only : case_setting
    use boundary_condition_m, only : bc_t, boundary_condition_velocity
    implicit none
    
    private

    type,extends(solver_fs) :: solver_fs_imp_t
        !!陰解法Fractional stepクラス. 粘性項のみをクランク･ニコルソン法で評価する半陰解法.
        !!@note おそらく, 第二段階の離散化が厳密ではない.
        !!@warning 現状, diffus_type = 2には対応出来ていない. diffus_type = 2を選ぶと恐らく破綻する(整合性がないので).
        type(mat_a) :: matrix_v
        real(dp),allocatable :: rhs_(:,:,:,:)
        integer(ip),private :: iter_max = 1000
        real(dp),private :: tolerance = 1.0d-16
        real(dp),private :: alpha = 1.0_dp
        contains
        procedure :: init => init_solver_2
        procedure predict_pseudo_velocity
        procedure set_parameter

    end type

    public solver_fs_imp_t, calc_pseudo_velocity_common_core
contains
subroutine set_parameter(this, iter_max, tolerance, alpha)
    !!中間速度を計算する反復法のためのパラメータを設定する.
    !!解法は現在SOR法を適用している. 
    class(solver_fs_imp_t),intent(inout) :: this
    integer(ip),intent(in) :: iter_max
        !!反復回数上限.
    real(dp),intent(in) :: tolerance
        !!反復の閾値.
    real(dp),intent(in) :: alpha
        !!SOR法の加速係数.

    this%iter_max = iter_max
    this%tolerance = tolerance
    this%alpha = alpha
end subroutine

subroutine init_solver_2(this, fld, grd, setting_case)
    !!ソルバを初期化する.
    class(solver_fs_imp_t),intent(inout) :: this
    type(fluid_field_t),intent(in) :: fld
    class(equil_mesh_t),intent(in) :: grd
    type(case_setting),intent(in) :: setting_case

    integer(ip) mx(3)

    mx(:) = grd%get_extents()
    this%method_name = "Fractional Step Implicit(CN)"

    call slvr_init_common(this, fld, grd, setting_case)
    
    call set_matrix_(this%matrix_v, grd%dsx, grd%dsy, grd%dsz, &
        grd%dv, grd%xc, grd%yc, grd%zc, setting_case%reynolds_number, setting_case%dt)
    allocate(this%rhs_(1:3,2:mx(1),2:mx(2),2:mx(3)))

    print "('   - velocity linear solver : SOR')"
    print "('   - tolerance              : ', g0)", this%tolerance
    print "('   - accelaration           : ', g0)", this%alpha

end subroutine

subroutine predict_pseudo_velocity(this, grid, fluid, del_t, bc_types)
    !!中間速度を計算する.
    class(solver_fs_imp_t),intent(inout) :: this
    class(equil_mesh_t),intent(in) :: grid
    type(fluid_field_t),intent(inout) :: fluid
    real(DP),intent(in) :: del_t
        !!時間刻み
    type(bc_t),intent(in) :: bc_types(:)
    
    real(dp),allocatable :: conv(:,:,:,:), diff(:,:,:,:)
    integer(ip) i, j, k, l, imx, jmx, kmx
    real(dp) rei

    call grid%get_extents_sub(imx, jmx, kmx)

    allocate(conv(3,2:imx,2:jmx,2:kmx))
    allocate(diff(3,2:imx,2:jmx,2:kmx))    
        
    call calc_convective_and_diffusive_flux(this, grid, fluid, this%v0, conv, diff)
    do k = 2, kmx
    do j = 2, jmx
    do i = 2, imx
        rei = 1.0_dp/fluid%kinetic_viscosity(i,j,k)
        this%rhs_(:,i,j,k) = this%v0(:,i,j,k)*grid%dv(i,j,k) & 
        - del_t*(conv(:,i,j,k) - 0.5_dp*rei*diff(:,i,j,k) - fluid%vol_force(:)*grid%dv(i,j,k))
    end do
    end do        
    end do

    call calc_pseudo_velocity_common_core(this, grid%get_extents(), fluid%velocity, bc_types)

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
              - this%matrix_v%a_nb(1,i,j,k)*v(l,i-1,j,k) - this%matrix_v%a_nb(2,i,j,k)*v(l,i+1,j,k) &
              - this%matrix_v%a_nb(3,i,j,k)*v(l,i,j-1,k) - this%matrix_v%a_nb(4,i,j,k)*v(l,i,j+1,k) &
              - this%matrix_v%a_nb(5,i,j,k)*v(l,i,j,k-1) - this%matrix_v%a_nb(6,i,j,k)*v(l,i,j,k+1))/this%matrix_v%a_p(i,j,k)
            ! point jacobi
            ! v(l,i,j,k) = (this%rhs_(l,i,j,k) &
            ! - this%matrix_v%aw*v_tmp(l,i-1,j  ,k  ) - this%matrix_v%ae*v_tmp(l,i+1,j  ,k  ) &
            ! - this%matrix_v%as*v_tmp(l,i  ,j-1,k  ) - this%matrix_v%an*v_tmp(l,i  ,j+1,k  ) &
            ! - this%matrix_v%ab*v_tmp(l,i  ,j  ,k-1) - this%matrix_v%at*v_tmp(l,i  ,j  ,k+1))/this%matrix_v%ap
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

! subroutine set_matrix_(mat, ds, dv, dx, re, dt, stencil_type)
subroutine set_matrix_(mat, dsx, dsy, dsz, dv, xc, yc, zc, re, dt)
    !!係数行列を設定する.
    type(mat_a),intent(inout) :: mat
    real(DP),intent(in) :: dsx(2:,2:), dsy(2:,2:), dsz(2:,2:)
        !!検査面の大きさ
    real(DP),intent(in) :: dv(2:,2:,2:)
        !!検査体積の大きさ
    real(DP),intent(in) :: xc(:)
        !!セル中心のx座標.
    real(DP),intent(in) :: yc(:)
        !!セル中心のy座標.
    real(DP),intent(in) :: zc(:)
        !!セル中心のz座標.
    real(DP),intent(in) :: re
        !!レイノルズ数
    real(DP),intent(in) :: dt
        !!時間刻み

    integer i, j, k
    integer imx, jmx, kmx
    real(dp) anb_, dx_

    !間接的に取得する.
    imx = ubound(dv, dim = 1)
    jmx = ubound(dv, dim = 2)
    kmx = ubound(dv, dim = 3)
    
    !@TODO large stencilの場合係数が微妙に変わる.それの実装.
    do k = 2, kmx
    do j = 2, jmx
    do i = 1, imx
        dx_ = -xc(i) + xc(i+1)
        anb_ = -0.5_dp*dt*dsx(j,k)/(re*dx_)
        if (i /= 1  ) mat%a_nb(2,i  ,j,k) = anb_
        if (i /= imx) mat%a_nb(1,i+1,j,k) = anb_
    end do
    end do        
    end do

    do k = 2, kmx
    do j = 1, jmx
    do i = 2, imx
        dx_ = -yc(j) + yc(j+1)
        anb_ = -0.5_dp*dt*dsy(i,k)/(re*dx_)
        if (j /= 1  ) mat%a_nb(4,i,j  ,k) = anb_
        if (j /= jmx) mat%a_nb(3,i,j+1,k) = anb_
    end do
    end do        
    end do

    do k = 1, kmx
    do j = 2, jmx
    do i = 2, imx
        dx_ = -zc(k) + zc(k+1)
        anb_ = -0.5_dp*dt*dsz(i,j)/(re*dx_)
        if (k /= 1  ) mat%a_nb(6,i,j,k  ) = anb_
        if (k /= kmx) mat%a_nb(5,i,j,k+1) = anb_
    end do
    end do        
    end do

    do k = 2, kmx
    do j = 2, jmx
    do i = 2, imx
        mat%a_p(i,j,k) = dv(i,j,k) - sum(mat%a_nb(:,i,j,k))
    end do
    end do        
    end do


end subroutine

end module fractional_step_implicit_m
