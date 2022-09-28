module fractional_step_m
    ! use,intrinsic :: iso_fortran_env, only : IP => int32, DP => real64
    use floating_point_parameter_m
    use fluid_field_m
    use mesh_m
    use setting_parameter_m, only : slv_setting, case_setting
    use boundary_condition_m, only : bc_t, boundary_condition_pressure
    implicit none
    private
    type :: mat_a
        !!線形方程式の係数.セルから見て順に東, 西, 北, 南, 上, 下側の係数
        real(DP),allocatable :: ae, aw, an, as, at, ab, ap
    end type

    type solver_fs
        !!フラクショナルステップソルバ. 
        !!時間精度はEuler陽解法.
        type(mat_a) :: coeffs_p
        character(:),allocatable :: convec_type
            !!移流項の評価方法. "ud1":1次精度風上法, "cd":2次精度中心差分, "qu":QUICK
        integer(ip),allocatable :: diffus_type
            !!拡散項の評価方法. 1:通常の中心差分, 2:ステンシルが一つ飛びの中心差分. ハウスコードに対応.
        logical :: correct_face_flux = .true.
        contains
        procedure :: init => init_solver 
        procedure :: predict_pseudo_velocity
        procedure calc_corrected_velocity
    end type

    public solver_fs, slvr_init_common, calc_convective_and_diffusive_flux, mat_a, cal_face_velocity
contains
subroutine init_solver(this, fld, grd, settings_slv, setting_case)
    !!ソルバを初期化する.
    class(solver_fs),intent(inout) :: this
    type(fluid_field_t),intent(in) :: fld
    type(rectilinear_mesh_t),intent(in) :: grd
    type(slv_setting),intent(in) :: settings_slv
    type(case_setting),intent(in) :: setting_case

    call slvr_init_common(this, fld, grd, settings_slv, setting_case)

end subroutine

subroutine slvr_init_common(this, fld, grd, settings_slv, setting_case)
    !!ソルバの初期化を行う. 圧力の係数行列設定, 移流/拡散スキームの設定.
    class(solver_fs),intent(inout) :: this
    type(fluid_field_t),intent(in) :: fld
    type(rectilinear_mesh_t),intent(in) :: grd
    type(slv_setting),intent(in) :: settings_slv
    type(case_setting),intent(in) :: setting_case

    call set_matrix_p(this%coeffs_p, grd%ds, grd%dv) !係数行列は1回だけ更新すれば良い.
    this%convec_type = settings_slv%conv_type
    this%diffus_type = settings_slv%diff_type
    this%correct_face_flux = settings_slv%correct_face_flux

    print "('   -------- Fractional Step Method ---------    ')"
    print "('      - 3D unsteady problem                     ')"
    print "('      - Convection : ', A)", merge("Up-wind", "Central", this%convec_type=="ud")
    print "('      - Diffusion  : ', A)", merge("compact", " large ", this%diffus_type==1)  
    print "('      - face flux correction : ', A)", merge("Yes","No ",this%correct_face_flux)
    print "('      - SOR linear solver for Poisson           ')"
    print "('   -----------------------------------------    ')"

end subroutine

subroutine predict_pseudo_velocity(this, extents, ds, dv, dx, del_t, re, force, v0, v, mi, mj, mk, dudr, bc_types)
    !!中間速度を計算する.
    class(solver_fs),intent(inout) :: this
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
        !!@note 初期状態から計算の場合は, ループ開始前に速度の初期値で更新しておく必要がある.
    real(dp),intent(in) :: dudr(:,:,:,:,:)
        !!速度勾配テンソル.
    type(bc_t),intent(in) :: bc_types(:)
    
    integer(ip) imx, jmx, kmx, i, j, k
    real(dp),allocatable :: conv(:,:,:,:), diff(:,:,:,:), rei

    imx = extents(1)
    jmx = extents(2)
    kmx = extents(3)

    allocate(conv(3,2:imx,2:jmx,2:kmx))
    allocate(diff(3,2:imx,2:jmx,2:kmx))

    !元々のコード. 移流/粘性項の評価が選択式に対応していないことに注意. 
    ! call fs_prediction_v2(imx, jmx, kmx, ds, dv, dx, del_t, re, force, v0, v, mi, mj, mk)

    !フラックスを作ってから後で右辺を計算する. メモリ量に不安がある場合は今まで通り上の方法で. 
    call calc_convective_and_diffusive_flux(this, extents, ds, dv, dx, re, v0, mi, mj, mk, dudr, conv, diff)

    rei = 1.0_dp/re
    do k = 2, kmx
    do j = 2, jmx
    do i = 2, imx
       
        v(:,i,j,k) = v0(:,i,j,k) - del_t*(conv(:,i,j,k) - rei*diff(:,i,j,k) - force(:)*dv)/dv

    end do
    end do
    end do



end subroutine

subroutine calc_corrected_velocity(this, extents, ds, dv, dx, del_t, mi, mj, mk, p, v,setting, p_ic, bc_types, diverged)
    !!圧力Poisson方程式を解いて速度とフラックスを修正する.
    class(solver_fs),intent(in) :: this
    integer(IP),intent(in) :: extents(3)
    real(DP),intent(in) :: ds(3), dv, dx(3)
    real(DP),intent(in) :: del_t
    real(DP),intent(inout) :: mi(:,2:,2:), mj(2:,:,2:), mk(2:,2:,:)
    real(DP),intent(inout) :: p(:,:,:)
    real(DP),intent(inout) :: v(:,:,:,:)
    type(slv_setting),intent(in) :: setting
    real(DP),intent(in) :: p_ic
    type(bc_t),intent(in) :: bc_types(6)
    logical,intent(out) :: diverged

    integer(IP) :: imx, jmx, kmx

    imx = extents(1)
    jmx = extents(2)
    kmx = extents(3)

    call cal_face_velocity(imx, jmx, kmx, v, mi, mj, mk)

    call fs_poisson_v2(imx, jmx, kmx, ds, dv, dx, del_t, this%coeffs_p, mi, mj, mk, p, &
                       setting, p_ic, bc_types, diverged)
    if(diverged) return
    !Poisson方程式の収束判定をパスした時点で圧力境界条件は満たされている.

    call fs_correction_v2(imx, jmx, kmx, ds, dv, del_t, p, v)

    !面の流束の速度補正. 
    if ( this%correct_face_flux ) then
        call fs_correction_face_flux(imx, jmx, kmx, ds, dv, del_t, p, mi, mj, mk)
    else
        !圧力で補正しない場合. mi, mj, mkは単に速度の線形補間として取り扱う. 
        !その場合ここで更新する必要がある. vはこの時点で新しい段階の速度なので, mi, mj, mkも新しい速度の線形補間となる.
        call cal_face_velocity(imx, jmx, kmx, v, mi, mj, mk)
    end if

end subroutine


!======================================================================

subroutine cal_face_velocity(imx, jmx, kmx, v, mi, mj, mk)
    !!部分段階速度から面の部分段階速度へ内挿する.
    implicit none
    integer(ip),intent(in) :: imx, jmx, kmx
    real(DP),intent(in) :: v(:,:,:,:)
    real(DP),intent(out) :: mi(:,2:,2:), mj(2:,:,2:), mk(2:,2:,:)
    integer(ip) i, j, k
    
    do K = 2, kmx
    do J = 2, jmx
    do I = 1, imx
        mi(i,j,k) = 0.5_dp*(v(1,i,j,k) + v(1,i+1,j,k))
    end do
    end do
    end do

    do K = 2, kmx
    do J = 1, jmx
    do I = 2, imx
        mj(i,j,k) = 0.5_dp*(v(2,i,j,k) + v(2,i,j+1,k))
    end do
    end do
    end do

    do K = 1, kmx
    do J = 2, jmx
    do I = 2, imx
        mk(i,j,k) = 0.5_dp*(v(3,i,j,k) + v(3,i,j,k+1))
    end do
    end do
    end do


end subroutine

subroutine fs_prediction_v2(imx, jmx, kmx, ds, dv, dx, del_t, re, force, v0, v, mi, mj, mk)
    !!フラクショナルステップの第一段階. 
    !!@deprecated : true
    !!@warning このプロセスはもう用いられない. 代わりに`calc_convective_and_diffusive_flux`でフラックスを作成してから計算する.
    !!時間精度: 1st order Euler (explicit)
    !!空間精度: 1st order upwind for convection, 2nd order central difference for diffusion
    integer(IP),intent(in) :: imx, jmx, kmx
    real(DP),intent(in) :: ds(3), dv, dx(3)
    real(DP),intent(in) :: del_t
    real(DP),intent(in) :: re
    real(DP),intent(in) :: force(3)
    real(DP),intent(in) :: v0(:,:,:,:)
    real(DP),intent(out) :: v(:,:,:,:)
    real(DP),intent(in) :: mi(:,2:,2:), mj(2:,:,2:), mk(2:,2:,:)

    real(DP) me, mw, mn, ms, mt, mb
    real(DP) conv(3), diff(3)
    real(DP) rei

    integer(IP) i, j, k

    rei = 1.0_dp/re
    do k = 2, kmx
    do j = 2, jmx
    do i = 2, imx
        !----convection----
        !::::mass flux
        me =  mi(i  ,j  ,k  )*ds(1)
        mw = -mi(i-1,j  ,k  )*ds(1)
        mn =  mj(i  ,j  ,k  )*ds(2)
        ms = -mj(i  ,j-1,k  )*ds(2)
        mt =  mk(i  ,j  ,k  )*ds(3)
        mb = -mk(i  ,j  ,k-1)*ds(3)

        !::::1st order upwind
        conv(:) = 0.5_dp*(me + abs(me))*v0(:,i  ,j  ,k  ) + 0.5_dp*(me - abs(me))*v0(:,i+1,j  ,k  ) &
                + 0.5_dp*(mw - abs(mw))*v0(:,i-1,j  ,k  ) + 0.5_dp*(mw + abs(mw))*v0(:,i  ,j  ,k  ) &
                + 0.5_dp*(mn + abs(mn))*v0(:,i  ,j  ,k  ) + 0.5_dp*(mn - abs(mn))*v0(:,i  ,j+1,k  ) &
                + 0.5_dp*(ms - abs(ms))*v0(:,i  ,j-1,k  ) + 0.5_dp*(ms + abs(ms))*v0(:,i  ,j  ,k  ) &
                + 0.5_dp*(mt + abs(mt))*v0(:,i  ,j  ,k  ) + 0.5_dp*(mt - abs(mt))*v0(:,i  ,j  ,k+1) &
                + 0.5_dp*(mb - abs(mb))*v0(:,i  ,j  ,k-1) + 0.5_dp*(mb + abs(mb))*v0(:,i  ,j  ,k  )

        !----diffusion----
        !::::2nd order central
        diff(:) = (v0(:,i-1,j  ,k  ) - 2.0_dp*v0(:,i  ,j  ,k  ) + v0(:,i+1,j  ,k  ))*ds(1)/dx(1) &
                + (v0(:,i  ,j-1,k  ) - 2.0_dp*v0(:,i  ,j  ,k  ) + v0(:,i  ,j+1,k  ))*ds(2)/dx(2) &
                + (v0(:,i  ,j  ,k-1) - 2.0_dp*v0(:,i  ,j  ,k  ) + v0(:,i  ,j  ,k+1))*ds(3)/dx(3)
        

        v(:,i,j,k) = v0(:,i,j,k) - del_t*(conv(:) - rei*diff(:) - force(:)*dv)/dv

    end do
    end do
    end do

end subroutine

subroutine calc_convective_and_diffusive_flux(this, extents, ds, dv, dx, re, v0, mi, mj, mk, dudr, convec, diff)
    !!右辺を計算するのに必要な移流と拡散流束を計算する.
    !!@note 陰解法に必要な右辺項は外部で計算する.
    class(solver_fs),intent(in) :: this
    integer(IP),intent(in) :: extents(3)
        !!内部セルの各方向の最大番号.
    real(DP),intent(in) :: ds(3), dv, dx(3)
        !!面積, 体積, 格子幅
    real(DP),intent(in) :: re
        !!レイノルズ数
    real(DP),intent(in) :: v0(:,:,:,:)
        !!既知段階速度.
    real(DP),intent(in) :: mi(:,2:,2:), mj(2:,:,2:), mk(2:,2:,:)
        !!面流束.
    real(dp),intent(in) :: dudr(:,:,:,:,:)
        !!速度勾配テンソル
    real(dp),intent(out) :: convec(:,2:,2:,2:)
        !!移流流束. 空間内部のみに存在するので, インデックスは2から始まる.
    real(dp),intent(out) :: diff(:,2:,2:,2:)
        !!拡散流束. 空間内部のみに存在するので, インデックスは2から始まる.
    real(DP) me, mw, mn, ms, mt, mb
    real(dp),dimension(3) :: dif_e, dif_w, dif_n, dif_s, dif_t, dif_b
    real(DP) rei


    integer(IP) i, j, k, ld, imx, jmx, kmx

    ld = this%diffus_type

    imx = extents(1)
    jmx = extents(2)
    kmx = extents(3)

    rei = 1.0_dp/re
    do k = 2, kmx
    do j = 2, jmx
    do i = 2, imx
        !----convection----
        !::::mass flux
        me =  mi(i  ,j  ,k  )*ds(1)
        mw = -mi(i-1,j  ,k  )*ds(1)
        mn =  mj(i  ,j  ,k  )*ds(2)
        ms = -mj(i  ,j-1,k  )*ds(2)
        mt =  mk(i  ,j  ,k  )*ds(3)
        mb = -mk(i  ,j  ,k-1)*ds(3)

        select case(this%convec_type) !計算効率的には良くないが, そこまで効率を求めていないので
        case("ud")
            !::::1st order upwind
            convec(:,i,j,k) = 0.5_dp*(me + abs(me))*v0(:,i  ,j  ,k  ) + 0.5_dp*(me - abs(me))*v0(:,i+1,j  ,k  ) &
                                + 0.5_dp*(mw - abs(mw))*v0(:,i-1,j  ,k  ) + 0.5_dp*(mw + abs(mw))*v0(:,i  ,j  ,k  ) &
                                + 0.5_dp*(mn + abs(mn))*v0(:,i  ,j  ,k  ) + 0.5_dp*(mn - abs(mn))*v0(:,i  ,j+1,k  ) &
                                + 0.5_dp*(ms - abs(ms))*v0(:,i  ,j-1,k  ) + 0.5_dp*(ms + abs(ms))*v0(:,i  ,j  ,k  ) &
                                + 0.5_dp*(mt + abs(mt))*v0(:,i  ,j  ,k  ) + 0.5_dp*(mt - abs(mt))*v0(:,i  ,j  ,k+1) &
                                + 0.5_dp*(mb - abs(mb))*v0(:,i  ,j  ,k-1) + 0.5_dp*(mb + abs(mb))*v0(:,i  ,j  ,k  )
        case("cd")
            !::::2nd order central difference
            convec(:,i,j,k) = 0.5_dp*mw*(v0(:,i-1,j  ,k  ) + v0(:,i,j,k)) + 0.5_dp*me*(v0(:,i,j,k) + v0(:,i+1,j  ,k  )) &
                                 + 0.5_dp*ms*(v0(:,i  ,j-1,k  ) + v0(:,i,j,k)) + 0.5_dp*mn*(v0(:,i,j,k) + v0(:,i,  j+1,k  )) &
                                 + 0.5_dp*mb*(v0(:,i  ,j  ,k-1) + v0(:,i,j,k)) + 0.5_dp*mt*(v0(:,i,j,k) + v0(:,i,  j  ,k+1))
        
        case default
            error stop this%convec_type//" not supported."
        end select
        
        !----diffusion----
        !::::2nd order central
        if ( ld == 1 ) then
            dif_w(:) = (  v0(:,i-1,j  ,k  ) - v0(:,i  ,j  ,k  ))/dx(1)
            dif_e(:) = (- v0(:,i  ,j  ,k  ) + v0(:,i+1,j  ,k  ))/dx(1)
            dif_s(:) = (  v0(:,i  ,j-1,k  ) - v0(:,i  ,j  ,k  ))/dx(2)
            dif_n(:) = (- v0(:,i  ,j  ,k  ) + v0(:,i  ,j+1,k  ))/dx(2)
            dif_b(:) = (  v0(:,i  ,j  ,k-1) - v0(:,i  ,j  ,k  ))/dx(3)
            dif_t(:) = (- v0(:,i  ,j  ,k  ) + v0(:,i  ,j  ,k+1))/dx(3)                
        else if ( ld == 2 ) then
            !面勾配*法線を計算する. セル勾配から面へ内挿する. 境界では片側差分となっている.
            dif_w(:) = -(dudr(:,1,i-1,j  ,k  ) + dudr(:,1,i  ,j  ,k  ))*0.5_dp
            dif_e(:) =  (dudr(:,1,i  ,j  ,k  ) + dudr(:,1,i+1,j  ,k  ))*0.5_dp
            dif_s(:) = -(dudr(:,2,i  ,j-1,k  ) + dudr(:,2,i  ,j  ,k  ))*0.5_dp
            dif_n(:) =  (dudr(:,2,i  ,j  ,k  ) + dudr(:,2,i  ,j+1,k  ))*0.5_dp
            dif_b(:) = -(dudr(:,3,i  ,j  ,k-1) + dudr(:,3,i  ,j  ,k  ))*0.5_dp
            dif_t(:) =  (dudr(:,3,i  ,j  ,k  ) + dudr(:,3,i  ,j  ,k+1))*0.5_dp

        end if

        diff(:,i,j,k) = (dif_e(:) + dif_w(:))*ds(1) + (dif_n(:) + dif_s(:))*ds(2) + (dif_t(:) + dif_b(:))*ds(3)
        
        !右辺のフラックスは例えば次のように作成する.
        ! rhs_flux(:,i,j,k) = v0(:,i,j,k)*dv - del_t*(conv(:) - rei*diff(:) - force(:)*dv)

    end do
    end do
    end do

end subroutine

subroutine set_matrix_p(coeffs, ds, dv)!, imx, jmx, kmx)
    !!圧力Poisson方程式の係数行列の設定.
    type(mat_a),intent(out) :: coeffs
    real(DP),intent(in) :: ds(3), dv
    ! integer(IP),intent(in) :: imx, jmx, kmx

    ! integer(IP) ic

    coeffs%ae = ds(1)*ds(1)/dv
    coeffs%aw = ds(1)*ds(1)/dv
    
    coeffs%an = ds(2)*ds(2)/dv
    coeffs%as = ds(2)*ds(2)/dv
    
    coeffs%at = ds(3)*ds(3)/dv
    coeffs%ab = ds(3)*ds(3)/dv

    coeffs%ap = - 2.0_dp*(ds(1)*ds(1) + ds(2)*ds(2) + ds(3)*ds(3))/dv

end subroutine

subroutine fs_poisson_v2(imx, jmx, kmx, ds, dv, dx, del_t, coeffs, mi, mj, mk, p, setting, p_ic, bc_types, diverged)
    !!gauss-seidel法により圧力を求める.
    use,intrinsic :: ieee_arithmetic, only : ieee_is_nan
    integer(IP),intent(in) :: imx, jmx, kmx
    real(DP),intent(in) :: ds(3), dv, dx(3)
    real(DP),intent(in) :: del_t
    type(mat_a),intent(in) :: coeffs
    real(DP),intent(in) :: mi(:,2:,2:), mj(2:,:,2:), mk(2:,2:,:)
    real(DP),intent(inout) :: p(:,:,:)
    type(slv_setting),intent(in) :: setting
    real(DP),intent(in) :: p_ic
    type(bc_t),intent(in) :: bc_types(6)
    logical,intent(out) :: diverged

    integer(IP) itr
    real(DP) resid_, den_
    real(DP) b
    real(DP) me, mw, mn, ms, mt, mb
    real(DP),allocatable :: p0(:,:,:)

    integer(IP) :: itr_mx
    real(DP) :: tol
    real(DP) :: alp
    integer(IP) icmx

    integer(IP) i, j, k
    integer(IP) extents_(3)

    extents_(1) = imx
    extents_(2) = jmx
    extents_(3) = kmx

    allocate(p0, source = p)
    itr_mx = setting%itr_max
    tol = setting%tolerance
    alp = setting%relax_coeff

    icmx = (imx - 1)*(jmx - 1)*(kmx - 1)

    p0(:,:,:) = p_ic
    p(:,:,:) = p_ic

    diverged = .false.
    do itr = 1, itr_mx
        do k = 2, kmx
        do j = 2, jmx
        do i = 2, imx
            me = mi(i  ,j  ,k  )*ds(1)
            mw = mi(i-1,j  ,k  )*ds(1)
            mn = mj(i  ,j  ,k  )*ds(2)
            ms = mj(i  ,j-1,k  )*ds(2)
            mt = mk(i  ,j  ,k  )*ds(3)
            mb = mk(i  ,j  ,k-1)*ds(3)
            b = (me - mw + mn - ms + mt - mb)/del_t

            p(i,j,k) = (1.0_dp - alp)*p(i,j,k) + alp*(b - coeffs%aw*p(i-1,j  ,k  ) - coeffs%ae*p(i+1,j  ,k  ) &
                                                        - coeffs%as*p(i  ,j-1,k  ) - coeffs%an*p(i  ,j+1,k  ) &
                                                        - coeffs%ab*p(i  ,j  ,k-1) - coeffs%at*p(i  ,j  ,k+1))/coeffs%ap

        end do
        end do
        end do

        call boundary_condition_pressure(extents_, p, dx, bc_types) !仮想セルの圧力を更新する.

        !収束判定.
        resid_ = 0.0_dp
        den_ = 0.0_dp
        do k = 2, kmx
        do j = 2, jmx
        do i = 2, imx
            resid_ = resid_ + (p(i,j,k) - p0(i,j,k))**2.0_dp
            den_ = den_ + (p(i,j,k))**2.0_dp
        end do
        end do
        end do

        if ( den_ == 0.0_dp ) exit !初期値や条件によってはゼロのままも有り得る. NANを防ぐ.
            
        resid_ = sqrt(resid_/icmx)/sqrt(den_/icmx)

        if ( resid_ <= tol .or. itr >= itr_mx) then
            print "('pressure(',i0,') converged, resid = ',g0)", itr, resid_
            ! print "('exit.')"
            exit
        end if
        
        if ( resid_ > huge(0.0_dp)) then
            print "(A)", "pressure diverged."
            diverged = .true.
            return
        end if

        if ( any(ieee_is_nan(p))) then
            print"(A)", "NAN detected in pressure."
            diverged = .true.
            return
        end if

        if ( (itr > 1000 .and. mod(itr,1000) == 0) ) then
            print "('pressure(',i0,'), resid = ',g0)", itr, resid_
        end if
        
        p0(:,:,:) = p(:,:,:)

    end do


end subroutine

subroutine fs_correction_v2(imx, jmx, kmx, ds, dv, del_t, p, v)
    !!求めた圧力で速度を修正し, n+1段階の速度を求める.
    integer(IP),intent(in) :: imx, jmx, kmx
    real(DP),intent(in) :: ds(3), dv
    real(DP),intent(in) :: del_t
    real(DP),intent(in) :: p(:,:,:)
    real(DP),intent(inout) :: v(:,:,:,:)
        !!中間段階の速度(1st step 終了後の速度).

    integer(IP) i, j, k

    do k = 2, kmx
    do j = 2, jmx
    do i = 2, imx
        !セル中心の速度
        v(1,i,j,k) = v(1,i,j,k) - 0.5_dp*del_t*ds(1)*(p(i+1,j  ,k  ) - p(i-1,j  ,k  ))/dv
        v(2,i,j,k) = v(2,i,j,k) - 0.5_dp*del_t*ds(2)*(p(i  ,j+1,k  ) - p(i  ,j-1,k  ))/dv
        v(3,i,j,k) = v(3,i,j,k) - 0.5_dp*del_t*ds(3)*(p(i  ,j  ,k+1) - p(i  ,j  ,k-1))/dv
    end do
    end do
    end do

end subroutine

subroutine fs_correction_face_flux(imx, jmx, kmx, ds, dv, del_t, p, mi, mj, mk)
    !!求めた圧力で速度を修正し, n+1段階の速度を求める.
    integer(IP),intent(in) :: imx, jmx, kmx
    real(DP),intent(in) :: ds(3), dv
    real(DP),intent(in) :: del_t
    real(DP),intent(in) :: p(:,:,:)
    real(DP),intent(inout) :: mi(:,2:,2:), mj(2:,:,2:), mk(2:,2:,:)

    integer(IP) i, j, k
    !!面の流束を圧力で修正する.

        !面速度 = 面に垂直な流速
    do k = 2, kmx
    do j = 2, jmx
    do i = 1, imx
        mi(i,j,k) = mi(i,j,k) - del_t*ds(1)*(p(i+1,j  ,k  ) - p(i  ,j  ,k  ))/dv
    end do
    end do
    end do

    do k = 2, kmx
    do j = 1, jmx
    do i = 2, imx
        mj(i,j,k) = mj(i,j,k) - del_t*ds(2)*(p(i  ,j+1,k  ) - p(i  ,j  ,k  ))/dv
    end do
    end do
    end do

    do k = 1, kmx
    do j = 2, jmx
    do i = 2, imx
        mk(i,j,k) = mk(i,j,k) - del_t*ds(3)*(p(i  ,j  ,k+1) - p(i,j  ,k  ))/dv
    end do
    end do
    end do
end subroutine
end module fractional_step_m