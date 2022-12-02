module les_solver_m
    use floating_point_parameter_m
    use fractional_step_m, only : calc_convective_and_diffusive_flux
    use fractional_step_implicit_m
    use fluid_field_m
    use mesh_m
    use boundary_condition_m
    use interpolation_m
    use output_m
    use vtk_field_data_m
    implicit none
    type, extends(solver_fs_imp_t), public :: solver_fs_imp_sm_t

    contains
        procedure predict_pseudo_velocity
    end type 
contains

subroutine predict_pseudo_velocity(this, grid, fluid, del_t, bc_types)
    class(solver_fs_imp_sm_t),intent(inout) :: this
    class(equil_mesh_t),intent(in) :: grid
    type(fluid_field_t),intent(inout) :: fluid
    real(DP),intent(in) :: del_t
        !!時間刻み
    type(bc_t),intent(in) :: bc_types(:)
    
    real(dp),allocatable :: conv(:,:,:,:), diff(:,:,:,:), nu_eddy(:,:,:), diff_e(:,:,:,:)
    integer(ip) i, j, k, l, imx, jmx, kmx
    real(dp) nu, re_
    type(attrib_data_holder_t),allocatable :: holders(:)

    call grid%get_extents_sub(imx, jmx, kmx)

    allocate(conv(3,2:imx,2:jmx,2:kmx))
    allocate(diff(3,2:imx,2:jmx,2:kmx), diff_e(3,2:imx,2:jmx,2:kmx))    
    allocate(nu_eddy(1:imx+1,1:jmx+1,1:kmx+1), source=1.0_dp/fluid%reynolds_number())

    ! 渦粘性係数を更新する. nu_eddyはdudr(x,y,z,t)の関数である. 
    ! 陰解法とすると係数行列を逐次更新しなければならないので煩雑となる.
    ! 壁面では渦粘性の効果は小さいので, 陰解法は分子粘性項のみとする. 
    re_ = fluid%reynolds_number()
    call les_0_Smagorinsky(fluid%dudr, re_, grid%dv, nu_eddy, grid%zc)
    ! fluid%kinetic_viscosity(2:imx,2:jmx,2:kmx) = 1.0d0/re_ + nu_eddy(2:imx,2:jmx,2:kmx)

    !粘性係数を仮想セルへ外挿する. 
    !壁面 ... 壁面の粘性係数はnu_eddy = 0なので 1/Reに等しい. 面でこうなるように外挿する.
    !周期境界面 ... 普通の周期境界と同じ.
    !固定流入 ... 面での速度勾配dudrから面の粘性係数を計算. その値を外挿.
    !NOTE :: 現在, サブルーチン開発中. それまで仮想セルは1/Reとする. その場合の誤差など記録しておく.

    call calc_convective_and_diffusive_flux(this, grid, fluid, this%v0, conv, diff)
    call calc_eddy_viscosity_term(this, grid, fluid%dudr, nu_eddy, this%v0, diff_e)

    do k = 2, kmx
    do j = 2, jmx
    do i = 2, imx
        ! nu = fluid%kinetic_viscosity(i,j,k)
        this%rhs_(:,i,j,k) = this%v0(:,i,j,k)*grid%dv(i,j,k) & 
        - del_t*(conv(:,i,j,k) - 0.5_dp*diff(:,i,j,k) - diff_e(:,i,j,k) - fluid%vol_force(:)*grid%dv(i,j,k))
    end do
    end do        
    end do

    call calc_pseudo_velocity_common_core(this, grid%get_extents(), fluid%velocity, bc_types)
    
    !debug
    !渦粘性係数を出力.
    allocate(holders(1))
    holders(1)%name = "Eddy_Viscosity"
    call holders(1)%register_scalar(nu_eddy, [2,2,2], [imx, jmx, kmx])
    call writeout_single_vtk_recti_grid("test_eddy_", 9999, grid%get_extents(), grid%xp, grid%yp, grid%zp, holders=holders)

end subroutine
    
subroutine les_0_Smagorinsky(dudr, re_const, dv, nu, zc)
    !!スマゴリンスキーモデルにより渦粘性係数を計算する. 平行平板間乱流専用.
    !! \nu_{eddy} = (C_SM*\Delta)^{2}*(2 S_{ij}S_{ij})**1.5
    real(dp),intent(in) :: dudr(:,:,:,:,:)
        !!速度勾配テンソル.
    real(dp),intent(in) :: re_const
        !!無次元化された動粘性係数.
    real(dp),intent(in) :: dv(2:,2:,2:)
        !!体積. フィルタ幅に用いる.
    real(dp),intent(inout) :: nu(:,:,:)
        !!更新される粘性係数.
    real(dp),intent(in) :: zc(:)
        !!頂点z座標.
    integer i, j, k
    integer ubnds(3)
    real(dp) delta_, dij(3,3), sij(3,3), reinv_, sij_norm
    real(dp),parameter :: C_SM = 0.173_dp
    
    real(dp) :: damping_
        !!van Driestの減衰関数(1 - exp(-C*z+/A))により壁面漸近挙動を再現する.
        !!壁面方向はz座標とする.
    real(dp) :: A0 = 25.0d0, dz_pls
    integer kmid
    ubnds(:) = ubound(nu) - 1

    reinv_ = 1.0_dp/re_const
    kmid = (2 + ubnds(3))/2
    do k = 2, ubnds(3)
    do j = 2, ubnds(2)
    do i = 2, ubnds(1)
        delta_ = dv(i,j,k)**(1./3.0_dp)
        dij(:,:) = dudr(:,:,i,j,k)
        sij(:,:) = 0.5_dp*(dij + transpose(dij))

        sij_norm = sqrt(sum(2.0_dp*sij*sij))
        !damping function
        !k = 1, k = kmxの壁を想定して減衰関数をかける.
        !z+ = z/l+ = Lz'/(ν/u+) = Lu+/ν*z' = Re*z'
        !減衰関数のz+ は壁(z+ = 0)から中心(z+ = Re*1)までの値の範囲.
        if ( k <= kmid ) then
            dz_pls = re_const*abs(zc(2) - zc(k))
        else if (k >= kmid) then
            dz_pls = re_const*abs(zc(k) - zc(ubnds(3)))
        end if

        damping_ = 1.0d0 - exp(-dz_pls/A0)

        nu(i,j,k) = (damping_*C_SM*delta_)**2.0_dp*sij_norm

    end do
    end do
    end do

end subroutine

subroutine calc_eddy_viscosity_term(this, grid, dudr, nu_eddy, v0, diff)
    !!渦粘性モデル項を計算する.
    use interpolation_m
    class(solver_fs_imp_sm_t),intent(in) :: this
    class(equil_mesh_t),intent(in) :: grid
    real(dp),intent(in) :: dudr(:,:,:,:,:)
    real(dp),intent(in) :: nu_eddy(:,:,:)
    real(DP),intent(in) :: v0(:,:,:,:)
        !!既知段階速度.
    real(dp),intent(out) :: diff(:,2:,2:,2:)
        !!拡散流束. 空間内部のみに存在するので, インデックスは2から始まる.
    real(DP) me, mw, mn, ms, mt, mb
    real(dp),dimension(3) :: dif_e, dif_w, dif_n, dif_s, dif_t, dif_b
    real(DP) nu_f(6)
        !1:w, 2:e, 3:s, 4:n, 5:b, 6:t

    integer(IP) i, j, k, ld, imx, jmx, kmx
    real(dp) ds(3)

    ld = this%settings%diff_type
    
    call grid%get_extents_sub(imx, jmx, kmx)
    
    do k = 2, kmx
    do j = 2, jmx
    do i = 2, imx
        nu_f(1) = harmonic_like_interpolate(nu_eddy(i-1,j,k), nu_eddy(i  ,j,k), grid%dx(i-1), grid%dx(i  ))
        nu_f(2) = harmonic_like_interpolate(nu_eddy(i  ,j,k), nu_eddy(i+1,j,k), grid%dx(i  ), grid%dx(i+1))
        nu_f(3) = harmonic_like_interpolate(nu_eddy(i,j-1,k), nu_eddy(i,j  ,k), grid%dy(j-1), grid%dy(j  ))
        nu_f(4) = harmonic_like_interpolate(nu_eddy(i,j  ,k), nu_eddy(i,j+1,k), grid%dy(j  ), grid%dy(j+1))
        nu_f(5) = harmonic_like_interpolate(nu_eddy(i,j,k-1), nu_eddy(i,j,k  ), grid%dz(k-1), grid%dz(k  ))
        nu_f(6) = harmonic_like_interpolate(nu_eddy(i,j,k  ), nu_eddy(i,j,k+1), grid%dz(k  ), grid%dz(k+1))
        nu_f(:) = 1.0_dp/nu_f(:)

        !----diffusion----
        !::::2nd order central
        if ( ld == 1 ) then
            dif_w(:) = (  v0(:,i-1,j  ,k  ) - v0(:,i  ,j  ,k  ))/(-grid%xc(i-1) + grid%xc(i  ))
            dif_e(:) = (- v0(:,i  ,j  ,k  ) + v0(:,i+1,j  ,k  ))/(-grid%xc(i  ) + grid%xc(i+1))
            dif_s(:) = (  v0(:,i  ,j-1,k  ) - v0(:,i  ,j  ,k  ))/(-grid%yc(j-1) + grid%yc(j  ))
            dif_n(:) = (- v0(:,i  ,j  ,k  ) + v0(:,i  ,j+1,k  ))/(-grid%yc(j  ) + grid%yc(j+1))
            dif_b(:) = (  v0(:,i  ,j  ,k-1) - v0(:,i  ,j  ,k  ))/(-grid%zc(k-1) + grid%zc(k  ))
            dif_t(:) = (- v0(:,i  ,j  ,k  ) + v0(:,i  ,j  ,k+1))/(-grid%zc(k  ) + grid%zc(k+1))
    
        else if ( ld == 2 ) then
            !面勾配*法線を計算する. セル勾配から面へ内挿する. 境界では片側差分となっている.
            dif_w(:) = -(dudr(:,1,i-1,j  ,k  ) + dudr(:,1,i  ,j  ,k  ))*0.5_dp
            dif_e(:) =  (dudr(:,1,i  ,j  ,k  ) + dudr(:,1,i+1,j  ,k  ))*0.5_dp
            dif_s(:) = -(dudr(:,2,i  ,j-1,k  ) + dudr(:,2,i  ,j  ,k  ))*0.5_dp
            dif_n(:) =  (dudr(:,2,i  ,j  ,k  ) + dudr(:,2,i  ,j+1,k  ))*0.5_dp
            dif_b(:) = -(dudr(:,3,i  ,j  ,k-1) + dudr(:,3,i  ,j  ,k  ))*0.5_dp
            dif_t(:) =  (dudr(:,3,i  ,j  ,k  ) + dudr(:,3,i  ,j  ,k+1))*0.5_dp
    
        end if
    
        diff(:,i,j,k) = (nu_f(1)*dif_w(:) + nu_f(2)*dif_e(:))*ds(1) + &
                        (nu_f(3)*dif_s(:) + nu_f(4)*dif_n(:))*ds(2) + &
                        (nu_f(5)*dif_b(:) + nu_f(6)*dif_t(:))*ds(3)
    end do
    end do        
    end do

end subroutine

end module les_solver_m