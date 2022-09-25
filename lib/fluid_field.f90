module fluid_field_m
    ! use,intrinsic :: iso_fortran_env, only : IP => int32, DP => real64
    use floating_point_parameter_m
    implicit none
    type fluid_field_t
        real(DP),allocatable :: velocity(:,:,:,:)
            !!流体の速度. 第一引数は成分(1:x,2:y,3:z)
            !!@note 仮想セルには, 面へ内挿したとき境界面の値となるように外挿された値が入っている.
        real(DP),allocatable :: mflux_i(:,:,:), mflux_j(:,:,:), mflux_k(:,:,:) 
            !界面の速度(=界面に垂直な速度成分)
        real(DP),allocatable :: pressure(:,:,:)
        real(dp),allocatable :: dudr(:,:,:,:,:)
            !!速度勾配テンソル.  第一引数は速度成分, 第二引数は微分の方向(x,y,z).
            !!dudr(:,1,i,j,k) = [dudx, dvdx, dwdx]
            !!@note 仮想セルにも勾配は配置されている. 2つのセル勾配から線形内挿したとき境界面の勾配となるように外挿されている.

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
    allocate(this%dudr(1:3,1:3,1:imx+1,1:jmx+1,1:kmx+1), source = 0.0_dp)
end subroutine

!====================================================================
subroutine calc_gradient_tensor(extents, dv, ds, dx, velocity, dudr)
    !!速度勾配テンソルを計算する. Green-Gaussの勾配を用いる.
    !!@note 速度勾配テンソルは転置されていない. かくんお
    integer(ip),intent(in) :: extents(3)
    real(dp),intent(in) :: dv
    real(dp),intent(in) :: ds(3)
    real(dp),intent(in) :: dx(3)
    real(dp),intent(in) :: velocity(:,:,:,:)
    real(dp),intent(inout) :: dudr(:,:,1:,1:,1:)

    integer(ip) i, j, k, l, igc(2), ic(2), inb, jnb, knb
    integer(ip) imx, jmx, kmx
    real(dp) uf(3)

    imx = extents(1)
    jmx = extents(2)
    kmx = extents(3)

    dudr(:,:,:,:,:) = 0.0_dp

    !i-dir. (face-loop)
    ! do k = 2, kmx
    ! do j = 2, jmx
    ! do i = 1, imx
    !     !east and west
    !     uf(1:3) = 0.5_dp*(velocity(1:3,i,j,k) + velocity(1:3,i+1,j,k))*ds(1:3) ![uf*dy*dz, vf*dz*dx, wf*dx*dy]
    !     if ( i /= 1   )  dudr(1:3,1,i  ,j,k) = dudr(1:3,1,i  ,j,k) + uf(:) !仮想セルに勾配値を持たせないためi=1は飛ばす.
    !     if ( i /= imx )  dudr(1:3,1,i+1,j,k) = dudr(1:3,1,i+1,j,k) - uf(:)
    !     !法線のy,z成分ny, nzはゼロなので, 残りの成分はゼロ.
    ! end do
    ! end do        
    ! end do

    ! !j-dir. (face-loop)
    ! do k = 2, kmx
    ! do j = 1, jmx
    ! do i = 2, imx
    !     !north and south
    !     uf(1:3) = 0.5_dp*(velocity(1:3,i,j,k) + velocity(1:3,i,j+1,k))*ds(1:3) ![uf*dy*dz, vf*dz*dx, wf*dx*dy]
    !     if ( j /= 1   )  dudr(1:3,2,i,j  ,k) = dudr(1:3,2,i,j  ,k) + uf(:) !仮想セルに勾配値を持たせないためi=1は飛ばす.
    !     if ( j /= jmx )  dudr(1:3,2,i,j+1,k) = dudr(1:3,2,i,j+1,k) - uf(:)
    !     !法線のz,x成分nz, nxはゼロなので, 残りの成分はゼロ.
    ! end do
    ! end do        
    ! end do

    ! !k-dir. (face-loop)
    ! do k = 1, kmx
    ! do j = 2, jmx
    ! do i = 2, imx
    !     !top and bottom
    !     uf(1:3) = 0.5_dp*(velocity(1:3,i,j,k) + velocity(1:3,i,j,k+1))*ds(1:3) ![uf*dy*dz, vf*dz*dx, wf*dx*dy]
    !     if ( k /= 1   )  dudr(1:3,3,i,j,k  ) = dudr(1:3,3,i,j,k  ) + uf(:) !仮想セルに勾配値を持たせないためi=1は飛ばす.
    !     if ( k /= kmx )  dudr(1:3,3,i,j,k+1) = dudr(1:3,3,i,j,k+1) - uf(:)
    !     !法線のx,y成分nx, nyはゼロなので, 残りの成分はゼロ.
    ! end do
    ! end do        
    ! end do

    ! dudr(:,:,:,:,:) = dudr(:,:,:,:,:)/dv

    !直交格子を仮定して簡略化. 幅2dxの中心差分と等価であるため.
    do k = 2, kmx
    do j = 2, jmx
    do i = 2, imx
        dudr(:,1,i,j,k) = (-velocity(:,i,j,k) + velocity(:,i+1,j  ,k  ))/dx(1)*0.5_dp
        dudr(:,2,i,j,k) = (-velocity(:,i,j,k) + velocity(:,i  ,j+1,k  ))/dx(2)*0.5_dp
        dudr(:,3,i,j,k) = (-velocity(:,i,j,k) + velocity(:,i  ,j  ,k+1))/dx(3)*0.5_dp
    end do
    end do        
    end do

    !dudr for ghost cell
    !粘性フラックスのため. large stencilの場合面の流束を, セル勾配の面への補間として計算する. そのため, 
    !仮想セルにも勾配を持たせておけば場合分けが不要となり都合が良い.
    !なお評価は, 線形補間した面の勾配が(uc - ub)/(dx/2)となるようにする. つまり仮想セルの勾配では無い. 
    !i-direction.
    igc(1) = 1     ;  ic(1) = 2
    igc(2) = imx+1 ;  ic(2) = imx
    do l = 1, 2
        i = igc(l)
        inb = ic(l)
        do k = 2, kmx
        do j = 2, jmx
            uf(:) = (- velocity(:,i,j,k) + velocity(:,inb,j,k))/dx(:)*real(-2*l+3, dp)
            !ufを使いまわしている. ここでは勾配として扱う.
            !velocity(:,i,j,k)が仮想セルug. ug = 2*ub - ucなので(uc-ug)/dx = (uc - 2ub+uc)/dx = (uc - ub)/(dx/2) の片側差分. 
            !最後の部分で, i=1とi=imx+1のときで符号が変わらないように揃える.
            dudr(:,1,i,j,k) = 2.0_dp*uf(:) - dudr(:,1,inb,j,k)
            !線形補間で面に内挿したときにufとなるように外挿する. 
        end do        
        end do    
    end do

    !j-direction.
    igc(1) = 1     ;  ic(1) = 2
    igc(2) = jmx+1 ;  ic(2) = jmx
    do l = 1, 2
        j = igc(l)
        jnb = ic(l)
        do k = 2, kmx
        do i = 2, imx
            uf(:) = (- velocity(:,i,j,k) + velocity(:,i,jnb,k))/dx(:)*real(-2*l+3, dp)
            !ufを使いまわしている. ここでは勾配として扱う.
            !velocity(:,i,j,k)が仮想セルug. ug = 2*ub - ucなので(uc-ug)/dx = (uc - 2ub+uc)/dx = (uc - ub)/(dx/2) の片側差分. 
            dudr(:,2,i,j,k) = 2.0_dp*uf(:) - dudr(:,2,i,jnb,k)
            !線形補間で面に内挿したときにufとなるように外挿する. 
        end do        
        end do    
    end do

    !k-direction.
    igc(1) = 1     ;  ic(1) = 2
    igc(2) = kmx+1 ;  ic(2) = kmx
    do l = 1, 2
        k = igc(l)
        knb = ic(l)
        do j = 2, jmx
        do i = 2, imx
            uf(:) = (- velocity(:,i,j,k) + velocity(:,i,j,knb))/dx(:)*real(-2*l+3, dp)
            !ufを使いまわしている. ここでは勾配として扱う.
            !velocity(:,i,j,k)が仮想セルug. ug = 2*ub - ucなので(uc-ug)/dx = (uc - 2ub+uc)/dx = (uc - ub)/(dx/2) の片側差分. 
            !最後の部分で, i=1とi=imx+1のときで符号が変わらないように揃える.
            dudr(:,3,i,j,k) = 2.0_dp*uf(:) - dudr(:,3,i,j,knb)
            !線形補間で面に内挿したときにufとなるように外挿する. 
        end do        
        end do    
    end do


end subroutine
end module fluid_field_m