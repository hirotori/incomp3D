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
        real(DP),allocatable :: xc(:)
            !!セル中心x座標.
        real(DP),allocatable :: yc(:)
            !!セル中心y座標.
        real(DP),allocatable :: zc(:)
            !!セル中心z座標.
        real(DP),allocatable :: xp(:)
            !!節点x座標.
        real(DP),allocatable :: yp(:)
            !!節点y座標.
        real(DP),allocatable :: zp(:)
            !!節点z座標.
        !-----private----
        real(dp),private :: dx_const(3) = 0.0d0
        real(dp),private :: dv_const = 0.0d0
        real(dp),private :: ds_const(3) = 0.0d0
        ! logical,private :: const_mesh = .true.
        contains
        procedure,private :: init_equil_spaced_mesh
        procedure,public :: init_from_data
        generic,public :: init => init_equil_spaced_mesh, init_from_data
        procedure,public :: get_extents
        procedure,public :: get_extents_sub
        procedure,public :: get_cell_count
        procedure,public :: alloc_arrays
        procedure,public :: calc_geometry
        procedure,public :: calc_ghost_cell_centers
        procedure,public,non_overridable :: get_equil_dx
        procedure,public :: print_self
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

    call this%alloc_arrays()
    
    call create_equil_spaced_mesh(imx, jmx, kmx, this%dx_const, lengths, this%xp, this%yp, this%zp)
    call calc_geometric_center(imx, jmx, kmx, this%xp, this%yp, this%zp, this%xc, this%yc, this%zc)

    call calc_equil_mesh_geometry(this)

    this%dx(:) = this%dx_const(1)
    this%dy(:) = this%dx_const(2)
    this%dz(:) = this%dx_const(3)
    
    this%dsx(2:,2:) = this%ds_const(1)
    this%dsy(2:,2:) = this%ds_const(2)
    this%dsz(2:,2:) = this%ds_const(3)

    call this%calc_ghost_cell_centers()
    
    this%dv(:,:,:) = this%dv_const

    print "('-- Equil-spaced Mesh created. --')"
    call this%print_self()

end subroutine

subroutine init_from_data(this, imx, jmx, kmx, x, y, z)
    !!格子をメッシュデータから構築する.
    class(equil_mesh_t),intent(inout) :: this
    integer,intent(in) :: imx, jmx, kmx
    real(dp),allocatable,intent(inout) :: x(:), y(:), z(:)

end subroutine

function get_equil_dx(this) result(dx_eq)
    class(equil_mesh_t),intent(in) :: this
    real(dp) dx_eq(3)

    dx_eq(:) = this%dx_const(:)

end function

subroutine print_self(this)
    !!コンソールに格子の情報を表示する.
    class(equil_mesh_t),intent(in) :: this

    integer imx, jmx, kmx

    call this%get_extents_sub(imx, jmx, kmx)

    print "('range')"
    print "(1x,g0, ' <= x <= ', g0)", minval(this%xc(2:imx)), maxval(this%xc(2:imx))
    print "(1x,g0, ' <= y <= ', g0)", minval(this%yc(2:jmx)), maxval(this%yc(2:jmx))
    print "(1x,g0, ' <= z <= ', g0)", minval(this%zc(2:kmx)), maxval(this%zc(2:kmx))
    print "('range including ghost cell ::')"
    print "(1x,g0, ' <= x <= ', g0)", minval(this%xc(:)), maxval(this%xc(:))
    print "(1x,g0, ' <= y <= ', g0)", minval(this%yc(:)), maxval(this%yc(:))
    print "(1x,g0, ' <= z <= ', g0)", minval(this%zc(:)), maxval(this%zc(:))
    print "('-- Geometric property --')"
    print "(1x,g0, ' <= ds (x) <= ', g0)", minval(this%dsx), maxval(this%dsx)
    print "(1x,g0, ' <= ds (y) <= ', g0)", minval(this%dsy), maxval(this%dsy)
    print "(1x,g0, ' <= ds (z) <= ', g0)", minval(this%dsz), maxval(this%dsz)
    print "(1x,g0, ' <= dv     <= ', g0)", minval(this%dv), maxval(this%dv)

end subroutine

subroutine alloc_arrays(this)
    class(equil_mesh_t),intent(inout) :: this

    !節点座標の割り付けは,メソッドinit_from_dataを優先する. それ以外はここで割り付ける
    if ( .not. allocated(this%xp) ) allocate(this%xp(this%imax))
    if ( .not. allocated(this%yp) ) allocate(this%yp(this%jmax))
    if ( .not. allocated(this%zp) ) allocate(this%zp(this%kmax)) 
    
    allocate(this%xc(1:this%imax+1))
    allocate(this%yc(1:this%jmax+1))
    allocate(this%zc(1:this%kmax+1))
    
    allocate(this%dx(1:this%imax+1), source = 0.0_dp)
    allocate(this%dy(1:this%jmax+1), source = 0.0_dp)
    allocate(this%dz(1:this%kmax+1), source = 0.0_dp)

    allocate(this%dsx(2:this%jmax,2:this%kmax), source = 0.0_dp)
    allocate(this%dsy(2:this%imax,2:this%kmax), source = 0.0_dp)
    allocate(this%dsz(2:this%imax,2:this%jmax), source = 0.0_dp)

    allocate(this%dv(2:this%imax,2:this%jmax,2:this%kmax), source = 0.0_dp)
end subroutine

subroutine create_equil_spaced_mesh(imx, jmx, kmx, dxs, l_xyz, xp, yp, zp)
    !!等間隔直交格子を生成する.
    integer(IP),intent(in) :: imx, jmx, kmx
    real(DP),intent(out) :: dxs(3)
    real(DP),intent(in) :: l_xyz(3)
    real(dp),intent(inout) :: xp(:), yp(:), zp(:)

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

    !等間隔直交格子. x座標は添え字j,kに依存しない.
    do i = 1, imx
        xp(i) = (i - 1)*dxs(1)
    end do

    do j = 1, jmx
        yp(j) = (j - 1)*dxs(2)
    end do

    do k = 1, kmx
        zp(k) = (k - 1)*dxs(3)
    end do

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


subroutine calc_geometric_center(imx, jmx, kmx, xp, yp, zp, xc, yc, zc)
    !!セルの重心を計算する.
    !!@note 等間隔直交格子あるいは不等間隔格子にのみ対応している.
    integer(ip),intent(in) :: imx, jmx, kmx
    real(dp),intent(in) :: xp(:), yp(:), zp(:)
    real(dp),intent(inout) :: xc(:), yc(:), zc(:)

    integer i, j, k

    do i = 2, imx
        xc(i) = 0.5_dp*(xp(i-1) + xp(i))
    end do

    do j = 2, jmx
        yc(j) = 0.5_dp*(yp(j-1) + yp(j))
    end do

    do k = 2, kmx
        zc(k) = 0.5_dp*(zp(k-1) + zp(k))
    end do

end subroutine

subroutine calc_ghost_cell_centers(this)
    !!仮想セルのセル中心を計算する. 仮想セルの格子幅も計算される.
    !!@note 仮想セルの中心点は, 境界面を挟んで境界セルと面対称な位置に置かれる. 計算には格子幅を用いる.
    class(equil_mesh_t),intent(inout) :: this

    integer i, j, k
    ! real(dp) dx_(3)
    
    !i-left boundary
    i = 1
    this%dx(i) = this%dx(i+1)
    this%xc(i) = this%xc(i+1) - this%dx(i+1)

    !i-right boundary
    i = this%imax
    this%dx(i+1) = this%dx(i)
    this%xc(i+1) = this%xc(i) + this%dx(i)

    !j-left boundary
    !不等間隔直交格子のため, y方向のみ移動する.
    j = 1
    this%dy(j) = this%dy(j+1)
    this%yc(j) = this%yc(j+1) - this%dy(j+1)

    !j-right boundary
    j = this%jmax
    this%dy(j+1) = this%dy(j)
    this%yc(j+1) = this%yc(j) + this%dy(j)

    !k-left boundary
    !不等間隔直交格子のため, z方向のみ移動する.
    k = 1
    this%dz(k) = this%dz(k+1)
    this%zc(k) = this%zc(k+1) - this%dz(k+1)

    !k-right boundary
    k = this%kmax
    this%dz(k+1) = this%dz(k)
    this%zc(k+1) = this%zc(k) + this%dz(k)

end subroutine
end module mesh_m