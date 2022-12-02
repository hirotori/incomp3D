program meshgen
    use floating_point_parameter_m
    use output_m
    implicit none
    integer(ip),parameter :: icmax = 32, jcmax = 64, kcmax = 64
    integer(ip),parameter :: imax = icmax + 1, jmax = jcmax + 1, kmax = kcmax + 1
    real(dp),allocatable :: xp(:), yp(:), zp(:)
    real(dp) Lx, Ly, Lz, a_, dx, dy, dz, xi_, reynolds_number
    integer i, j, k
    character(128) basename, comment
    
    logical :: swap_ = .false.

    Lx = 6.4_dp
    Ly = 3.2_dp
    Lz = 2.0_dp

    dx = Lx/real(imax - 1, dp)
    dy = Ly/real(jmax - 1, dp)
    dz = Lz/real(kmax - 1, dp)

    reynolds_number = 395.0d0
    
    allocate(xp(imax))
    allocate(yp(jmax))
    allocate(zp(kmax))

    do i = 1, imax
        xp(i) = real(i - 1, dp)*dx
    end do
    
    do j = 1, jmax
        yp(j) = real(j - 1, dp)*dy
    end do

    a_ = 0.95_dp
    do k = 1, kmax
        xi_ = -1.0_dp + 2.0_dp*real(k - 1, dp)/real(kmax - 1, dp)
        zp(k) = tanh(xi_*atanh(a_))/a_ + 1.0_dp
        ! zp(k) = real(k - 1, dp)*dz
    end do
    !zp(1)が不正な値っぽいので(1d-018のオーダー?), 0に変える.
    zp(1) = 0.0_dp
    print "('Mesh Parameter')"
    print "(g0,'<=x<=',g0)", xp(1), xp(imax)
    print "(g0,'<=y<=',g0)", yp(1), yp(jmax)
    print "(g0,'<=z<=',g0)", zp(1), zp(kmax)
    print "('dx+ = ',g0)", dx*reynolds_number
    print "('dy+ = ',g0)", dy*reynolds_number
    print "('dzmin+ = ',g0)", (zp(2) - zp(1))*reynolds_number
    print "('dzmax+ = ',g0)", (zp(kmax/2) - zp(kmax/2-1))*reynolds_number

    if ( swap_ ) then
        write(basename, "(A,'_',i0,'x',i0,'x',i0)") "channel", imax, kmax, jmax
        comment = "3d channel flow"
        call writeout_as_msh_rectil_format(basename, [imax,kmax,jmax], trim(adjustl(comment)), xp, zp, yp)
    else
        write(basename, "(A,'_',i0,'x',i0,'x',i0)") "channel", imax, jmax, kmax
        comment = "3d channel flow"
        call writeout_as_msh_rectil_format(basename, [imax,jmax,kmax], trim(adjustl(comment)), xp, yp, zp)
    end if

end program meshgen