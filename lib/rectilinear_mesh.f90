module rectilinear_mesh_m
    use floating_point_parameter_m
    use mesh_m, only : equil_mesh_t, calc_geometric_center
    implicit none
    type,extends(equil_mesh_t) :: rectilinear_mesh_t
        !!version : 1.0.0
        !!不等間隔格子クラス.
        contains
        procedure :: init_from_data => construct_from_meshdata
        procedure :: calc_geometry
        procedure,private :: calc_grid_spacing
        procedure,private :: calc_surface_areas_and_volumes
        procedure,private :: check_geometry

    end type
    
contains
subroutine construct_from_meshdata(this, imx, jmx, kmx, r)
    !!クラスを初期化する. 
    !!格子ファイルを読み取り, メンバ変数の構築を行う.
    class(rectilinear_mesh_t),intent(inout) :: this
    real(dp),allocatable,intent(inout) :: r(:,:,:,:)
        !!節点座標.
    integer(ip),intent(in) :: imx, jmx, kmx

    this%imax = imx
    this%jmax = jmx
    this%kmax = kmx

    call move_alloc(r, this%rp)

    !セル中心の計算.
    call calc_geometric_center(this%imax, this%jmax, this%kmax, this%rp, this%rc)

    print"(A)", "Mesh constructed."
    print "(g0, ' <= x <= ', g0)", minval(this%rc(1,2:,2:,2:)), maxval(this%rc(1,2:,2:,2:))
    print "(g0, ' <= y <= ', g0)", minval(this%rc(2,2:,2:,2:)), maxval(this%rc(2,2:,2:,2:))
    print "(g0, ' <= z <= ', g0)", minval(this%rc(3,2:,2:,2:)), maxval(this%rc(3,2:,2:,2:))


    call this%alloc_arrays()

    call this%calc_geometry()

end subroutine 

subroutine calc_geometry(this)
    !!幾何量を計算する.
    !!re-implemented from equil_mesh_t
    class(rectilinear_mesh_t),intent(inout) :: this

    call this%calc_grid_spacing()
    call this%calc_surface_areas_and_volumes()
    call this%check_geometry()

    print "('---- Geometric Properties ----')"
    print "('Sum Of Volume = ', g0.4)", sum(this%dv)
    print "('Max Volume    = ', g0.4, ', at ', *(i0,:,',',1x))", maxval(this%dv), maxloc(this%dv)
    print "('Min Volume    = ', g0.4, ', at ', *(i0,:,',',1x))", minval(this%dv), minloc(this%dv)
    print "('Min Grid Spacing for i-dir. :: ', g0)", minval(this%dx(:))
    print "('Min Grid Spacing for j-dir. :: ', g0)", minval(this%dy(:))
    print "('Min Grid Spacing for k-dir. :: ', g0)", minval(this%dz(:))
end subroutine


subroutine calc_grid_spacing(this)
    !!格子の各方向の格子幅を計算する.
    class(rectilinear_mesh_t),intent(inout) :: this
    
    integer(ip) i, j, k

    do i = 2, this%imax
        this%dx(i) = -this%rp(1,i-1,j,k) + this%rp(1,i,j,k)
    end do

    do j = 2, this%jmax
        this%dy(j) = -this%rp(2,i,j-1,k) + this%rp(2,i,j,k)
    end do

    do k = 2, this%kmax
        this%dz(k) = -this%rp(3,i,j,k-1) + this%rp(3,i,j,k)
    end do
    
end subroutine 


subroutine calc_surface_areas_and_volumes(this)
    !!メッシュの面積と体積を計算する.
    !!@note 事前に格子幅が計算されていなければならない.
    class(rectilinear_mesh_t),intent(inout) :: this

    integer i, j, k
    !i-dir.
    i = 1
    do k = 2, this%kmax
    do j = 2, this%jmax
        !i方向の面積は同じi-格子線上で共通.
        this%dsx(j,k) = this%dy(j)*this%dz(k)
    end do        
    end do

    !j-dir.
    j = 1
    do k = 2, this%kmax
    do i = 2, this%imax
        !j方向の面積は同じj-格子線上で共通.
        this%dsy(i,k) = this%dx(i)*this%dz(k)
    end do
    end do

    !k-dir.
    k = 1
    do j = 2, this%jmax
    do i = 2, this%imax
        !k方向の面積は同じk-格子線上で共通.
        this%dsz(i,j) = this%dx(i)*this%dy(j)
    end do
    end do

    !volume
    do k = 2, this%kmax
    do j = 2, this%jmax
    do i = 2, this%imax
        this%dv(i,j,k) = this%dx(i)*this%dy(j)*this%dz(k)
    end do
    end do
    end do
end subroutine

subroutine check_geometry(this)
    !!計算された幾何量の妥当性をチェックする.
    class(rectilinear_mesh_t),intent(in) :: this

    ! integer i, j, k
    character(*),parameter :: error_header = "check_geometry :: "
    if ( any(this%dv < 0) ) then
        error stop error_header//"negative volume exists."
    end if


end subroutine
    
end module rectilinear_mesh_m