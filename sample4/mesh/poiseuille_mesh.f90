program pois_mesh
    use floating_point_parameter_m
    use output_m
    implicit none
    
    integer(ip) icmax, jcmax, kcmax
    integer(ip) imax, jmax, kmax
    real(dp) Lx, Ly, Lz
    real(dp) dx, dz
    real(dp),allocatable :: xp(:), yp(:), zp(:)
    integer(ip) i, j, k
    real(dp) xi_, a_

    !PARAMETER
    icmax = 32
    jcmax = 16
    kcmax = 16
    Lx = 4.0_dp
    Ly = 2.0_dp
    Lz = 2.0_dp


    imax = icmax + 1
    jmax = jcmax + 1
    kmax = kcmax + 1
    dx = Lx/real(imax - 1, dp)
    dz = Lz/real(kmax - 1, dp)
    
    allocate(xp(imax))
    allocate(yp(jmax))
    allocate(zp(kmax))

    do i = 1, imax
        xp(i) = real(i - 1, dp)*dx
    end do
    
    do k = 1, kmax
        zp(k) = real(k - 1, dp)*dz
    end do

    do j = 1, jmax
        xi_ = -1.0_dp + 2.0_dp*real(j - 1, dp)/real(jmax - 1, dp)
        yp(j) = tanh(xi_*atanh(a_))/a_ + 1.0_dp
    end do

    block
        character(32) tmp_
        write(tmp_,"(A,'_', *(i0,:,'_'))") "poiseuille", icmax, jcmax, kcmax

        call writeout_as_msh_format(trim(adjustl(tmp_)), [imax,jmax,kmax], xp, yp, zp)

    end block
end program pois_mesh