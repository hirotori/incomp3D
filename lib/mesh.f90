module mesh_m
    ! use,intrinsic :: iso_fortran_env, only : IP => int32, DP => real64
    use floating_point_parameter_m
    implicit none
    type equil_mesh_t
        !!version: 1.1.0
        !!等間隔格子クラス.
        !!@note v.1.1.0より, 等間隔の場合でも, 面積, 体積, 格子幅は配列として持たせるようにした.
        integer(IP) imax, jmax, kmax
        real(DP),allocatable :: dsx(:,:)
            !!dsx(2:jmx,2:kmx) : i方向の面の面積.
        real(DP),allocatable :: dsy(:,:)
            !!dsy(2:imx,2:kmx) : j方向の面の面積.
        real(DP),allocatable :: dsz(:,:)
            !!dsz(2:imx,2:jmx) : k方向の面の面積.
        real(DP),allocatable :: dv(:,:,:)
            !!セルの体積.
        real(DP),allocatable :: dx(:)
            !!dx(2:imx) : セルのx方向格子幅.
        real(DP),allocatable :: dy(:)
            !!dy(2:jmx) : セルのy方向格子幅.
        real(DP),allocatable :: dz(:)
            !!dz(2:kmx) : セルのz方向格子幅.
        real(DP),allocatable :: rc(:,:,:,:)
            !!セル中心座標.
        real(DP),allocatable :: rp(:,:,:,:)
            !!節点座標.
        !-----private----
        real(dp),private :: dx_const(3) = 0.0d0
        real(dp),private :: dv_const = 0.0d0
        real(dp),private :: ds_const(3) = 0.0d0
        contains
        procedure,private :: init_equil_spaced_mesh
        procedure,public :: init_from_data
        generic,public :: init => init_equil_spaced_mesh, init_from_data
        procedure,public :: get_extents
        procedure,public :: get_extents_sub
        procedure,public :: get_cell_count
        procedure,public :: alloc_arrays
        procedure,public :: calc_geometry
        procedure,public,non_overridable :: get_equil_dx
    end type
    
contains
pure function get_extents(this) result(extents)
    !!格子点の個数を返す. 
    !!@note 現在のメッシュの取り扱いの上では, 内部セルのインデックス上限と同じである.
    class(equil_mesh_t),intent(in) :: this
    integer(ip),dimension(3) :: extents
        !!各方向の格子点. [imax, jmax, kmax]

    extents(1) = this%imax
    extents(2) = this%jmax
    extents(3) = this%kmax

end function

subroutine get_extents_sub(this, imx, jmx, kmx)
    !!格子点の個数を返す. 
    class(equil_mesh_t),intent(in) :: this
    integer(ip),intent(out) :: imx, jmx, kmx

    imx = this%imax
    jmx = this%jmax
    kmx = this%kmax

end subroutine

integer(ip) function get_cell_count(this)
    class(equil_mesh_t),intent(in) :: this

    get_cell_count = (this%imax - 1)*(this%jmax - 1)*(this%kmax - 1)

end function

subroutine init_equil_spaced_mesh(this, imx, jmx, kmx, lengths)
    !!等間隔直交格子メッシュを生成し初期化する.
    class(equil_mesh_t),intent(inout) :: this
    integer(IP),intent(in) :: imx, jmx, kmx
    real(DP),intent(in) :: lengths(3)

    this%imax = imx
    this%jmax = jmx
    this%kmax = kmx

    allocate(this%rc(3,1:imx,1:jmx,1:kmx))
    allocate(this%rp(3,1:imx,1:jmx,1:kmx))
    
    call create_equil_spaced_mesh(imx, jmx, kmx, this%dx_const, lengths, this%rc, this%rp)
    call calc_equil_mesh_geometry(this)

    call this%alloc_arrays()

    this%dx(:) = this%dx_const(1)
    this%dy(:) = this%dx_const(2)
    this%dz(:) = this%dx_const(3)
    
    this%dsx(2:,2:) = this%ds_const(1)
    this%dsy(2:,2:) = this%ds_const(2)
    this%dsz(2:,2:) = this%ds_const(3)
    
    allocate(this%dv(2:this%imax,2:this%jmax,2:this%kmax), source = this%dv_const)

end subroutine

subroutine init_from_data(this, imx, jmx, kmx, r)
    !!格子をメッシュデータから構築する.
    class(equil_mesh_t),intent(inout) :: this
    integer,intent(in) :: imx, jmx, kmx
    real(dp),allocatable,intent(inout) :: r(:,:,:,:)

end subroutine

function get_equil_dx(this) result(dx_eq)
    class(equil_mesh_t),intent(in) :: this
    real(dp) dx_eq(3)

    dx_eq(:) = this%dx_const(:)

end function

subroutine alloc_arrays(this)
    class(equil_mesh_t),intent(inout) :: this

    allocate(this%dx(2:this%imax), source = 0.0_dp)
    allocate(this%dy(2:this%jmax), source = 0.0_dp)
    allocate(this%dz(2:this%kmax), source = 0.0_dp)

    allocate(this%dsx(2:this%jmax,2:this%kmax), source = 0.0_dp)
    allocate(this%dsy(2:this%imax,2:this%kmax), source = 0.0_dp)
    allocate(this%dsz(2:this%imax,2:this%jmax), source = 0.0_dp)

end subroutine

subroutine create_equil_spaced_mesh(imx, jmx, kmx, dxs, l_xyz, rc, rp)
    !!等間隔直交格子を生成する.
    integer(IP),intent(in) :: imx, jmx, kmx
    real(DP),intent(out) :: dxs(3)
    real(DP),intent(in) :: l_xyz(3)
    real(DP),intent(out) :: rc(:,:,:,:), rp(:,:,:,:)

    integer(IP) i, j, k
    
    ! --------+--------+--------+--------+-- j+1
    !         |        |        |        |
    !         |   +    |   +    |   +    |
    !         |        |(I+1,J+1,K)      |
    ! --------+--------+--------+--------+-- j
    !         |        |        |        |
    !         |   +    |   +    |   +    |
    !         | (I,J,K)|(I+1,J,K)        |
    ! --------+--------+--------+--------+-- j-1        
    !        i-1       i       i+1      i+2         on plane at k=k
            
    dxs(1) = l_xyz(1)/real(imx - 1, DP)
    dxs(2) = l_xyz(2)/real(jmx - 1, DP)
    dxs(3) = l_xyz(3)/real(kmx - 1, DP)

    do k = 1, kmx
    do j = 1, jmx
    do i = 1, imx
        rp(1,i,j,k) = (i - 1)*dxs(1)
        rp(2,i,j,k) = (j - 1)*dxs(2)
        rp(3,i,j,k) = (k - 1)*dxs(3)
    end do
    end do                
    end do

    call calc_geometric_center(imx, jmx, kmx, rp, rc)

    print "('Mesh created.')"
    print "(g0, ' <= x <= ', g0)", minval(rc(1,2:,2:,2:)), maxval(rc(1,2:,2:,2:))
    print "(g0, ' <= y <= ', g0)", minval(rc(2,2:,2:,2:)), maxval(rc(2,2:,2:,2:))
    print "(g0, ' <= z <= ', g0)", minval(rc(3,2:,2:,2:)), maxval(rc(3,2:,2:,2:))

end subroutine

subroutine calc_geometry(this)
    !!幾何量を計算する.
    class(equil_mesh_t),intent(inout) :: this

    call calc_equil_mesh_geometry(this)

end subroutine

subroutine calc_equil_mesh_geometry(this)
    !!幾何量を計算する.
    class(equil_mesh_t),intent(inout) :: this

    !calculate surface area
    this%ds_const(1) = this%dx_const(2)*this%dx_const(3)
    this%ds_const(2) = this%dx_const(1)*this%dx_const(3)
    this%ds_const(3) = this%dx_const(1)*this%dx_const(2)
    !calculate cell volume
    this%dv_const = product(this%dx_const)

end subroutine


subroutine calc_geometric_center(imx, jmx, kmx, rp, rc)
    !!セルの重心を計算する.
    !!@note 等間隔直交格子あるいは不等間隔格子にのみ対応している.
    integer(ip),intent(in) :: imx, jmx, kmx
    real(dp),intent(in) :: rp(:,:,:,:)
    real(dp),intent(inout) :: rc(:,:,:,:)

    integer i, j, k

    do k = 2, kmx
    do j = 2, jmx
    do i = 2, imx
        rc(:,i,j,k) = 0.125_dp*(rp(:,i-1,j-1,k-1) + rp(:,i,j-1,k-1) &
                                + rp(:,i-1,j  ,k-1) + rp(:,i,j  ,k-1) &
                                + rp(:,i-1,j-1,k  ) + rp(:,i,j-1,k  ) &
                                + rp(:,i-1,j  ,k  ) + rp(:,i,j  ,k  ))
    end do
    end do                
    end do

end subroutine
end module mesh_m