module mesh_m
    ! use,intrinsic :: iso_fortran_env, only : IP => int32, DP => real64
    use floating_point_parameter_m
    implicit none
    type rectilinear_mesh_t
        integer(IP) imax, jmax, kmax
        real(DP) ds(3)
        real(DP) dv
        real(DP) dx(3)
        real(DP),allocatable :: rc(:,:,:,:)
        real(DP),allocatable :: rp(:,:,:,:)
        contains
        procedure,public :: init => init_mesh
        procedure,public :: get_extents
        procedure,public :: calc_geometry
    end type
    
contains
pure function get_extents(this) result(extents)
    !!格子点の個数を返す. 
    !!@note 現在のメッシュの取り扱いの上では, 内部セルのインデックス上限と同じである.
    class(rectilinear_mesh_t),intent(in) :: this
    integer(ip),dimension(3) :: extents
        !!各方向の格子点. [imax, jmax, kmax]

    extents(1) = this%imax
    extents(2) = this%jmax
    extents(3) = this%kmax

end function

subroutine init_mesh(this, imx, jmx, kmx, lengths)
    !!等間隔直交格子メッシュを生成し初期化する.
    class(rectilinear_mesh_t),intent(inout) :: this
    integer(IP),intent(in) :: imx, jmx, kmx
    real(DP),intent(in) :: lengths(3)

    this%imax = imx
    this%jmax = jmx
    this%kmax = kmx

    allocate(this%rc(3,1:imx,1:jmx,1:kmx))
    allocate(this%rp(3,1:imx,1:jmx,1:kmx))
    
    call create_mesh(imx, jmx, kmx, this%dx, lengths, this%rc, this%rp)
    call this%calc_geometry()

end subroutine

subroutine create_mesh(imx, jmx, kmx, dxs, l_xyz, rc, rp)
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

    print "('Mesh created.')"
    print "(g0, ' <= x <= ', g0)", minval(rc(1,2:,2:,2:)), maxval(rc(1,2:,2:,2:))
    print "(g0, ' <= y <= ', g0)", minval(rc(2,2:,2:,2:)), maxval(rc(2,2:,2:,2:))
    print "(g0, ' <= z <= ', g0)", minval(rc(3,2:,2:,2:)), maxval(rc(3,2:,2:,2:))

end subroutine

subroutine calc_geometry(this)
    !!幾何量を計算する.
    class(rectilinear_mesh_t),intent(inout) :: this

    !calculate surface area
    this%ds(1) = this%dx(2)*this%dx(3)
    this%ds(2) = this%dx(1)*this%dx(3)
    this%ds(3) = this%dx(1)*this%dx(2)
    !calculate cell volume
    this%dv = product(this%dx)

end subroutine
end module mesh_m