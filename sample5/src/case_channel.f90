module case_channel_m
    use floating_point_parameter_m
    use case_common_m
    use mesh_m
    use fluid_field_m
    use interpolation_m
    use vtk_field_data_m
    use output_m, only: writeout_single_vtk_recti_grid
    implicit none

    type,extends(case_common_t) :: case_channel_t

        contains
        procedure add_on_pre_process
        procedure phase_post_process
        procedure phase_writeout
    end type
    
contains

subroutine add_on_pre_process(this, grid, fld)
    !!追加で初期状態を管理したい場合に呼び出す.
    class(case_channel_t),intent(inout) :: this
    class(equil_mesh_t),intent(inout) :: grid
    type(fluid_field_t),intent(inout) :: fld

    real(dp) zpls, re_
    real(dp) :: alpha_, beta_
    integer i, j, k, unit, imax, jmax, kmax
    integer :: n_data, n, nj
    real(dp),allocatable :: yj(:,:), uup(:), vvp(:), wwp(:)
    real(dp) um, uvp, rand_v(3)

    alpha_ = (5.0d0*log(30.0d0)+1.0d0)/(2.0d0*log(6.0d0))
    beta_  = 0.20d0*exp(5.0d0/alpha_)
    re_ = fld%reynolds_number()
    call grid%get_extents_sub(imax, jmax, kmax)
    
    do k = 2, kmax
        zpls = re_*(grid%zc(k))
        if ( zpls <= 5.0d0 ) then
            fld%velocity(1,2:imax,2:jmax,k) = zpls
        else if ( 5.0d0   < zpls .and. zpls <= 30.0d0) then
            fld%velocity(1,2:imax,2:jmax,k) = alpha_*log(zpls*beta_)
        else if ( 30.0d0  < zpls .and. zpls <= 395.0d0) then
            fld%velocity(1,2:imax,2:jmax,k) = 2.5d0*log(zpls) + 5.5d0
        else if ( 395.0d0 < zpls .and. zpls <=  760.0d0) then
            fld%velocity(1,2:imax,2:jmax,k) = 2.5d0*log(790.0d0 - zpls) + 5.5d0
        else if ( 760.0d0 < zpls .and. zpls <= 785.0d0 ) then
            fld%velocity(1,2:imax,2:jmax,k) = alpha_*log((790.0d0 - zpls)*beta_)
        else if ( 785.0d0 < zpls ) then
            fld%velocity(1,2:imax,2:jmax,k) = 790.0d0 - zpls
        end if
    end do

    open(newunit=unit, file="ch395_turb_sta.txt", status="old")
        read(unit,*)
        n_data = 96
        allocate(yj(n_data*2,2), uup(n_data*2), vvp(n_data*2), wwp(n_data*2))
        do n = 1, n_data
            read(unit,*) nj, yj(n,1), um, uup(n), wwp(n)
        end do

        read(unit,*)
        read(unit,*)
        do n = 1, n_data
            read(unit,*) nj, yj(n,2), vvp(n), uvp
        end do

        do n = n_data+1, 2*n_data
            yj(n,:) = 2.0d0*395.0d0 - yj(2*n_data-n+1,:)
            uup(n) = uup(2*n_data-n+1)
            vvp(n) = vvp(2*n_data-n+1)
            wwp(n) = wwp(2*n_data-n+1)
        end do

    close(unit)
    !テスト.
    ! open(newunit=unit, file="test.txt")
    !     do n = 1, 2*n_data
    !         write(unit,"(*(g0,1x))") yj(n,1:2), uup(n), vvp(n), wwp(n)
    !     end do
    ! close(unit)

    !擾乱を加える. DNSデータは現在の格子間隔よりも細かい(96点).
    !そのため線形補間により間のデータを計算する.
    !注意) DNSのデータはyが壁方向.
    call random_seed()
    !擾乱の大きさの目安. 
    um = sum(fld%velocity(1,2,2,2:kmax))/real(kmax - 1, dp)
    print"('umean(init) = ', g0)", um
    do k = 2, kmax
    do j = 2, jmax
    do i = 2, imax
        zpls = re_*(grid%zc(k))
        call random_number(rand_v)
        ! rand_v = rand_v - 0.5
        do n = 1, 2*n_data - 1
            !uu+, ww+ に擾乱を加える.
            if ( yj(n,1) <= zpls .and. zpls <= yj(n+1,1) ) then
                fld%velocity(1,i,j,k) = fld%velocity(1,i,j,k) + &
                    sqrt(weighted_average_interpolate(uup(n), uup(n+1), zpls-yj(n,1), yj(n+1,1)-zpls))*rand_v(1)*um*0.05d0
                fld%velocity(2,i,j,k) = fld%velocity(2,i,j,k) + &
                    sqrt(weighted_average_interpolate(wwp(n), wwp(n+1), zpls-yj(n,1), yj(n+1,1)-zpls))*rand_v(2)*um*0.05d0
            end if
            !vv+ 〃
            if ( yj(n,2) <= zpls .and. zpls <= yj(n+1,2) ) then
                fld%velocity(3,i,j,k) = fld%velocity(3,i,j,k) + &
                    sqrt(weighted_average_interpolate(vvp(n), vvp(n+1), zpls-yj(n,2), yj(n+1,2)-zpls))*rand_v(3)*um*0.05d0
            end if
        end do
    end do
    end do
    end do

end subroutine

subroutine phase_post_process(this, grid, fld, step, time)
    !!現時間段階の計算が終わった後に呼び出される. 
    class(case_channel_t),intent(inout) :: this
    class(equil_mesh_t),intent(in) :: grid
    type(fluid_field_t),intent(in) :: fld
    integer(ip),intent(in) :: step
    real(dp),intent(in) :: time

    print "('post process is being called from **case_channel_t** ')"


end subroutine

subroutine phase_writeout(this, grid, fld, nstep)
    !!データの書き出し.
    !!@todo 書き出しサブルーチンの隔離.
    class(case_channel_t),intent(in) :: this
    class(equil_mesh_t),intent(in) :: grid
    type(fluid_field_t),intent(in) :: fld
    integer(ip),intent(in) :: nstep
    character(:),allocatable :: outdir

    outdir = this%output_directory()

    call writeout_(this, grid, fld, outdir//"/"//"output", nstep)


end subroutine

subroutine writeout_(this, grid, fld, basename, current_step)
    !!データを書き出す. フォーマットはvtkファイルとなっている.
    !!@note データ節約のため等間隔メッシュの場合格子座標データを必要としないvtk structured points フォーマットを使用している.
    class(case_channel_t),intent(in) :: this
    class(equil_mesh_t),intent(in) :: grid
    type(fluid_field_t),intent(in) :: fld
    character(*),intent(in) :: basename
    integer(ip),intent(in) :: current_step

    type(attrib_data_holder_t),allocatable :: holders(:)

    allocate(holders(3))
    
    holders(1)%name = "Velocity"
    call holders(1)%register_vector(fld%velocity, [2,2,2], grid%get_extents())
    holders(2)%name = "Pressure"
    call holders(2)%register_scalar(fld%pressure, [2,2,2], grid%get_extents())
    
    call writeout_single_vtk_recti_grid(basename, current_step, grid%get_extents(), &
                                            grid%xp, grid%yp, grid%zp, holders=holders)
    
end subroutine

end module case_channel_m