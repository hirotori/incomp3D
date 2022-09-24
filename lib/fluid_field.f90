module fluid_field_m
    ! use,intrinsic :: iso_fortran_env, only : IP => int32, DP => real64
    use floating_point_parameter_m
    implicit none
    type fluid_field_t
        real(DP),allocatable :: velocity(:,:,:,:)
        real(DP),allocatable :: mflux_i(:,:,:), mflux_j(:,:,:), mflux_k(:,:,:) 
            !界面の速度(=界面に垂直な速度成分)
        real(DP),allocatable :: pressure(:,:,:)
        ! real(DP) :: reynolds_number = 100.0_dp 
        ! real(DP) :: vol_force(3) = 0.0_dp
        !     !!体積力. 外力として.
        contains
        procedure,public :: init => init_field
    end type
contains
subroutine init_field(this, imx, jmx, kmx, ic_u, ic_p)
    class(fluid_field_t),intent(inout) :: this
    integer(IP) imx, jmx, kmx
    real(DP),intent(in) :: ic_u(3)
    real(DP),intent(in) :: ic_p

    integer(IP) l
    allocate(this%velocity(3,1:imx+1,1:jmx+1,1:kmx+1))
    do l = 1, 3
        this%velocity(l,:,:,:) = ic_u(l)        
    end do
    allocate(this%pressure(1:imx+1,1:jmx+1,1:kmx+1), source = ic_p)

    allocate(this%mflux_i(1:imx, 2:jmx, 2:kmx), source = 0.0_dp)
    allocate(this%mflux_j(2:imx, 1:jmx, 2:kmx), source = 0.0_dp)
    allocate(this%mflux_k(2:imx, 2:jmx, 1:kmx), source = 0.0_dp)

end subroutine

!====================================================================
subroutine calc_gradient_tensor(extents, dv, ds, velocity, dudr)
    !!速度勾配テンソルを計算する. Green-Gaussの勾配を用いる.
    !!@note 速度勾配テンソルは転置されていない. かくんお
    integer(ip),intent(in) :: extents(3)
    real(dp),intent(in) :: dv
    real(dp),intent(in) :: ds(3)
    real(dp),intent(in) :: velocity(:,:,:,:)
    real(dp),intent(inout) :: dudr(:,:,2:,2:,2:)

    integer(ip) i, j, k
    integer(ip) imx, jmx, kmx
    real(dp) uf(3)

    imx = extents(1)
    jmx = extents(2)
    kmx = extents(3)

    dudr(:,:,:,:,:) = 0.0_dp

    !i-dir. (face-loop)
    do k = 2, kmx
    do j = 2, jmx
    do i = 1, imx
        !east and west
        uf(1:3) = 0.5_dp*(velocity(1:3,i,j,k) + velocity(1:3,i+1,j,k))*ds(1:3) ![uf*dy*dz, vf*dz*dx, wf*dx*dy]
        if ( i /= 1   )  dudr(1:3,1,i  ,j,k) = dudr(1:3,1,i  ,j,k) + uf(:) !仮想セルに勾配値を持たせないためi=1は飛ばす.
        if ( i /= imx )  dudr(1:3,1,i+1,j,k) = dudr(1:3,1,i+1,j,k) - uf(:)
        !法線のy,z成分ny, nzはゼロなので, 残りの成分はゼロ.
    end do
    end do        
    end do

    !j-dir. (face-loop)
    do k = 2, kmx
    do j = 1, jmx
    do i = 2, imx
        !north and south
        uf(1:3) = 0.5_dp*(velocity(1:3,i,j,k) + velocity(1:3,i,j+1,k))*ds(1:3) ![uf*dy*dz, vf*dz*dx, wf*dx*dy]
        if ( j /= 1   )  dudr(1:3,2,i,j  ,k) = dudr(1:3,2,i,j  ,k) + uf(:) !仮想セルに勾配値を持たせないためi=1は飛ばす.
        if ( j /= jmx )  dudr(1:3,2,i,j+1,k) = dudr(1:3,2,i,j+1,k) - uf(:)
        !法線のz,x成分nz, nxはゼロなので, 残りの成分はゼロ.
    end do
    end do        
    end do

    !k-dir. (face-loop)
    do k = 1, kmx
    do j = 2, jmx
    do i = 2, imx
        !top and bottom
        uf(1:3) = 0.5_dp*(velocity(1:3,i,j,k) + velocity(1:3,i,j,k+1))*ds(1:3) ![uf*dy*dz, vf*dz*dx, wf*dx*dy]
        if ( k /= 1   )  dudr(1:3,3,i,j,k  ) = dudr(1:3,3,i,j,k  ) + uf(:) !仮想セルに勾配値を持たせないためi=1は飛ばす.
        if ( k /= kmx )  dudr(1:3,3,i,j,k+1) = dudr(1:3,3,i,j,k+1) - uf(:)
        !法線のx,y成分nx, nyはゼロなので, 残りの成分はゼロ.
    end do
    end do        
    end do
    
    dudr(:,:,:,:,:) = dudr(:,:,:,:,:)/dv

end subroutine
end module fluid_field_m