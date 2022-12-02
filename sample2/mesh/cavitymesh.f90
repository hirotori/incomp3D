program cvmesh
    use floating_point_parameter_m
    use output_m, only : writeout_as_msh_rectil_format
    implicit none
    integer(ip),parameter :: Nx = 31, Ny = 31, Nz = 1
        !!セルの個数.
    integer(ip),parameter :: imax = Nx + 1, jmax = Ny + 1, kmax = Nz + 1
    real(dp),parameter :: Lx = 1.0d0, Ly = 1.0d0, Lz = 0.1d0
    real(dp),parameter :: PI = -4.0d0*atan(-1.0d0)
    real(dp) dx, dy, dz
    real(dp),dimension(imax) :: xp
    real(dp),dimension(jmax) :: yp
    real(dp),dimension(kmax) :: zp
    integer(ip) i, j, k
    character(32) basename

    print "('Mesh Generator for 2D lid-driven cavity flow')"
    print "('imax = ', i0)", imax
    print "('jmax = ', i0)", jmax

    dx = Lx/real(imax-1, dp)
    dy = Ly/real(jmax-1, dp)
    dz = Lz/real(kmax-1, dp)

    do i = 1, imax
        xp(i) = 0.5d0*Lx*(1.0d0 - cos(PI*real(i-1, dp)/real(imax-1, dp)))
    end do

    do j = 1, jmax
        yp(j) = 0.5d0*Ly*(1.0d0 - cos(PI*real(j-1, dp)/real(jmax-1, dp)))
    end do

    do k = 1, kmax
        zp(k) = real(k-1, dp)*dz
    end do

    write(basename, "('cavity_',i0,'x',i0)") Nx, Ny

    call writeout_as_msh_rectil_format(basename, [imax,jmax,kmax],"Mesh for 2D cavity flow (dimension : 1.0 x 1.0)", &
        xp, yp, zp)

    print "(g0, '<= x <=', g0)", minval(xp), maxval(xp)
    print "(g0, '<= y <=', g0)", minval(yp), maxval(yp)
    print "(g0, '<= z <=', g0)", minval(zp), maxval(zp)

end program cvmesh