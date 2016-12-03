!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ju@sussex
! 11/2015
! FORTRAN routines to compute low level products
! Main purposes is interfacing with python
! See 1611.01446
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module correlation_functions

contains

    SUBROUTINE  gauleg(ngp, xabsc, weig)
        ! Given a number of gauss point n,
        ! this routines computes the GAUSS-LEGENDRE abscissas (cos) and
        ! weights (sintheta dtheta) for Gaussian Quadrature.
        ! /!\ The routine works for normalized lower and
        ! /!\ upper limits of integration -1.0 & 1.0
        ! /!\ therefore if you work between [0,pi],
        ! /!\ you need to invert the order of the output vectors.
        implicit none
        integer, parameter       :: I4B = 4
        integer, parameter       :: DP = 8
        real(DP), parameter      :: pi = 3.14159265359
        real(DP), parameter      :: EPS = 3.0d-15
        INTEGER  i, j, m
        REAL(DP)  p1, p2, p3, pp, z, z1
        INTEGER, INTENT(IN) :: ngp
        REAL(DP), INTENT(OUT) :: xabsc(ngp), weig(ngp)

        ! Roots are symmetric in the interval - so only need to find half of them
        m = (ngp + 1) / 2

        ! Loop over the desired roots
        do i = 1, m

            z = cos( pi * (i-0.25d0) / (ngp+0.5d0) )
            ! Starting with the above approximation to the ith root,
            ! We enter the main loop of refinement by NEWTON'S method
    100     p1 = 1.0d0
            p2 = 0.0d0
            !*  Loop up the recurrence relation to get the Legendre
            !*  polynomial evaluated at z                 */

            do j = 1, ngp
                p3 = p2
                p2 = p1
                p1 = ((2.0d0*j-1.0d0) * z * p2 - (j-1.0d0)*p3) / j
            enddo

            ! p1 is now the desired Legendre polynomial. We next compute pp,
            ! its derivative, by a standard relation involving also p2, the
            ! polynomial of one lower order.
            pp = ngp*(z*p1-p2)/(z*z-1.0d0)
            z1 = z
            z = z1 - p1/pp ! Newton's Method

            if (dabs(z-z1) .gt. EPS) GOTO  100

            xabsc(i) =  - z                     ! Roots will be bewteen -1.0 & 1.0
            xabsc(ngp+1-i) =  + z               ! and symmetric about the origin
            weig(i) = 2.0d0/((1.0d0-z*z)*pp*pp) ! Compute the weight and its
            weig(ngp+1-i) = weig(i)             ! symmetric counterpart

        end do     ! i loop

    End subroutine gauleg

    subroutine compute_d_ell_1m1(ntheta,lmax,d_ell_1m1)
        ! Compute the Wigner d-matrix (m1=1, m2=-1) for all ell
        ! and all beta, using simple cosine
        ! Rough translation of the C++ code in CLASS
        ! Not used here. Mostly used to perform checks
        ! for the general recursion formula.
        implicit none
        integer, parameter       :: I4B = 4
        integer, parameter       :: DP = 8
        real(DP), parameter      :: pi = 3.14159265359

        integer, intent(in)      :: ntheta, lmax
        real(DP), intent(out)    :: d_ell_1m1(0:ntheta - 1,0:lmax)

        integer(I4B)             :: l, theta
        real(DP)                 :: dl
        real(DP)                 :: cosvec(0:ntheta - 1), sinvec(0:ntheta - 1)
        real(DP)                 :: fac1(0:lmax),fac2(0:lmax),fac3(0:lmax),fac4(0:lmax)
        real(DP)                 :: dlm2, dlm1, dlp1

        do l=2,lmax
            dl = real(l,kind=DP)
            fac1(l) = SQRT((2*dl+3)/(2*dl+1))*(dl+1)*(2*dl+1)/(dl*(dl+2))
            fac2(l) = 1.0/(dl*(dl+1.))
            fac3(l) = SQRT((2*dl+3)/(2*dl-1))*(dl-1)*(dl+1)/(dl*(dl+2))*(dl+1)/dl
            fac4(l) = SQRT(2./(2*dl+3))
        enddo

        CALL generate_cosine_sine_vec(pi,ntheta,cosvec(0:ntheta - 1),sinvec(0:ntheta - 1))

        ! Last loop
        do theta=0,ntheta - 1
            d_ell_1m1(theta,0) = 0.d0
            dlm2 = (1.d0 - cosvec(theta))/2.d0 * SQRT(3.d0/2.d0)
            d_ell_1m1(theta,1) = dlm2 * SQRT(2.d0/3.d0)
            dlm1 = (1.0-cosvec(theta))/2.d0 * (2.d0*cosvec(theta) + 1.d0) * SQRT(5.d0/2.d0)
            d_ell_1m1(theta,2) = dlm1 * SQRT(2.d0/5.d0)

            do l=2,lmax - 1
                dl = real(l,kind=DP)
                dlp1 = fac1(l)*(cosvec(theta) + fac2(l))*dlm1 - fac3(l)*dlm2
                d_ell_1m1(theta,l + 1) = dlp1 * fac4(l)
                dlm2 = dlm1
                dlm1 = dlp1
            enddo
        enddo

    end subroutine

    subroutine generate_cosine_sine_vec(end,npoints,cosvec,sinvec)
        ! Computes vectors of cos(n*theta) and  sin(n*theta) using
        ! recursion formula (Moivre)
        implicit none
        integer, parameter       :: I4B = 4
        integer, parameter       :: DP = 8
        real(DP), parameter      :: pi = 3.14159265359

        integer(I4B), intent(in) :: npoints
        real(DP), intent(in)     :: end
        real(DP), intent(out)    :: cosvec(0:npoints-1), sinvec(0:npoints-1)

        integer(I4B)             :: idx
        real(DP)                 :: dtheta, costheta, sintheta

        ! Step
        dtheta = end / real(npoints,kind=DP)

        ! Initialization
        cosvec(0) = 1.0
        sinvec(0) = 0.0
        costheta = cos(dtheta)
        sintheta = sin(dtheta)

        ! Recursion
        do idx=1, npoints-1
            cosvec(idx) = costheta * cosvec(idx - 1) - sintheta * sinvec(idx - 1)
            sinvec(idx) = sintheta * cosvec(idx - 1) + costheta * sinvec(idx - 1)
        enddo

    end subroutine

    subroutine generate_d_wigner_recursion_pospos(accurate_integration,ntheta,lmax,m1,m2,d_ell_mm)
        ! Fast routine to compute the d-Wigner matrix for all \ell and all angle at fixed m values
        ! The routine is based on 3 recursions.
        ! The first two recursions (hidden) are used to initialize the third recursion
        ! The third recursion comes from spherical harmonics calculation.
        ! To get even more speed-up, cosine and sine are pre-computed and stored in the memory
        ! By choosing accurate_integration, you will use the roots of a Gauss-Legendre quadrature
        ! /!\ The routine computes d_ell_mmp for m>mp, and then apply symmetry if m<mp.
        ! /!\ The routine assumes m - mp = N, l >= N.
        implicit none
        integer, parameter       :: I4B = 4
        integer, parameter       :: DP = 8
        real(DP), parameter      :: pi = 3.14159265359

        integer, intent(in)      :: ntheta, lmax, m1, m2
        real(DP), intent(out)    :: d_ell_mm(0:ntheta - 1,0:lmax)
        logical, intent(in)      :: accurate_integration

        integer(I4B)             :: l, startell, theta, N, im1, im2, FACN1, i, FACN2
        real(DP)                 :: dl, dm1, dm2, dN, sign, di, fac1, halfntheta, fac2
        real(DP)                 :: cosvec(0:ntheta - 1), sinvec(0:ntheta - 1), weig(0:ntheta - 1)
        real(DP)                 :: cosvec_tmp(0:ntheta - 1)
        real(DP)                 :: halfcosvec(0:ntheta - 1), halfsinvec(0:ntheta - 1)
        real(DP)                 :: A(lmax), C(lmax), B(0:ntheta - 1,0:lmax)

        ! Check if m > mp
        ! If not, swap m values, and change sign of the matrix
        if (m2 .lt. 1 .and. m1 .lt. 1) then
            im1 = -m2
            im2 = -m1
        endif

        if (m2 .gt. m1) then
            sign = (-1.d0)**(m1 + m2)
            im1 = m2
            im2 = m1
        else
            sign = 1.d0
            im1 = m1
            im2 = m2
        endif


        dm1 = real(im1, kind=DP)
        dm2 = real(im2, kind=DP)

        ! Initialization
        ! Set accurate_integration to .True. if you are
        ! working with Gauss-legendre quadrature
        if (accurate_integration .eqv. .True.) then
            CALL gauleg(ntheta, cosvec, weig)
            do theta=0,ntheta - 1
                cosvec_tmp(theta) = cosvec(ntheta - 1 - theta)
                sinvec(theta) = SQRT(1.d0 - cosvec_tmp(theta)**2)
                halfcosvec(theta) = SQRT( (1.d0 + cosvec_tmp(theta)) / 2.d0)
                halfsinvec(theta) = SQRT( (1.d0 - cosvec_tmp(theta)) / 2.d0)
            enddo
            cosvec = cosvec_tmp
        else
            CALL generate_cosine_sine_vec(pi,ntheta,cosvec(0:ntheta - 1),sinvec(0:ntheta - 1))
            CALL generate_cosine_sine_vec(pi/2.d0,ntheta,halfcosvec(0:ntheta - 1),halfsinvec(0:ntheta - 1))
        endif


        ! Forbidden values l<m
        startell = MAX(ABS(im1),ABS(im2))
        if (ABS(im1) .gt. 0 .or. ABS(im2) .gt. 0) then
            do l=0,startell-1
                d_ell_mm(:,l)=0.d0
            enddo
        endif

        ! Pre-factors for the last recursion
        do l=MAX(2,MAX(ABS(im1),ABS(im2))),lmax
            dl = real(l, kind=DP)
            A(l) = dl * (2*dl - 1) / SQRT( (dl**2 - dm1**2) * (dl**2 - dm2**2) )
            C(l) = SQRT( ((dl - 1)**2 - dm1**2)*((dl - 1)**2 - dm2**2) ) / ( (2 * dl - 1) * (dl - 1) )
            do theta=0, ntheta - 1
                B(theta,l) = cosvec(theta) - dm1 * dm2 /( dl * (dl - 1) )
            enddo
        enddo

        halfntheta = ntheta/2.d0

        ! Recursion, see Eq. 64 from Blanco et al. 1997
        ! Initialize two first values for the recursion (two internal recursions)
        ! And then major loop
        N = ABS(im1-im2)
        dN = real(N, kind=DP)

        ! Prefactors for the initialization again...
        if (N .eq. 0) then
            FACN1 = 1
            FACN2 = 1
        else
            FACN1 = 1
            FACN2 = 1
            do i=0,N
                FACN1 = FACN1 * (i + 1)
            enddo
            do i=1,N
                FACN2 = FACN2 * i
            enddo
        endif

        do l=startell, lmax ! First loop on l
            dl = real(l, kind=DP)

            ! Prefactors again and again...
            if (N .eq. 0) then
                fac1 = 1.d0
                fac2 = 1.d0
            else
                fac1 = 1.d0
                fac2 = 1.d0
                do i=1,N
                    di = real(i, kind=DP)
                    fac1 = fac1 * (2*dl - di)
                enddo
                do i=0,N - 1
                    di = real(i, kind=DP)
                    fac2 = fac2 * (2*dl - di)
                enddo
            endif

            do theta=0,ntheta - 1 ! Second loop on the angle
                if (l .eq. startell) then ! initial term l - 2
                    if (N .eq. 0) then
                        d_ell_mm(theta,l) = halfcosvec(theta)**(2*l)
                    else
                        d_ell_mm(theta,l) = SQRT(fac2/FACN2) * halfcosvec(theta)**(2*l - N) * (-halfsinvec(theta))**N
                    endif

                elseif (l .eq. startell + 1) then ! initial term l - 1
                    if (N .eq. 0) then
                        d_ell_mm(theta,l) = (dl * cosvec(theta) - dl + 1.d0) * halfcosvec(theta)**(2*l - 2)
                    else
                        d_ell_mm(theta,l) = (dl*(cosvec(theta) - 1.d0) + dN + 1.d0) * &
                         & SQRT(fac1/FACN1) * &
                         & halfcosvec(theta)**(2*l-2-N) * (-halfsinvec(theta))**N
                    endif
                else ! all term l > startell + 2
                    d_ell_mm(theta,l) =  A(l) * (B(theta,l) * d_ell_mm(theta,l - 1) - C(l) * d_ell_mm(theta,l - 2))
                endif
                ! Apply sign convention
                d_ell_mm(theta,l) = sign * d_ell_mm(theta,l)
            enddo
        enddo

    end subroutine

    subroutine compute_correlation_function_temp(accurate_lensing,compute_tgradt,ntheta,clCMB,clpp,xi,lmin,lmax)
        ! Compute the correlation function xi in temperature for all angle
        ! see astro-ph/0601594, Eq. 9.19
        ! It can also compute the power spectrum of T \grad T (polarization not added yet).
        implicit none
        integer, parameter       :: I4B = 4
        integer, parameter       :: DP = 8
        real(DP), parameter      :: pi = 3.14159265359

        integer, intent(in)      :: ntheta, lmax, lmin
        real(DP), intent(in)     :: clCMB(0:lmax), clpp(0:lmax)
        logical, intent(in)      :: accurate_lensing, compute_tgradt
        real(DP), intent(out)    :: xi(0:ntheta - 1)
        real(DP)                 :: d_ell_00(0:ntheta - 1,0:lmax)
        real(DP)                 :: d_ell_11(0:ntheta - 1,0:lmax)
        real(DP)                 :: d_ell_2m2(0:ntheta - 1,0:lmax)
        real(DP)                 :: d_ell_1m1(0:ntheta - 1,0:lmax)
        real(DP)                 :: clg(0:ntheta - 1), clg2(0:ntheta - 1)
        real(DP)                 :: sigma_sq(0:ntheta - 1)
        real(DP)                 :: array_clCMB(0:lmax)
        real(DP)                 :: dl, prefac1, prefac2, dinfinity
        integer(I4B)             :: l, theta, offset, offset2
        dinfinity = 1.e10_dp

        if (compute_tgradt .eqv. .True.) then
            print *,'Computing spectrum of T.gradT'
            offset = 1
        else
            offset = 0

        endif

        CALL generate_d_wigner_recursion_pospos(accurate_lensing,ntheta,lmax,1,1,d_ell_11)
        CALL generate_d_wigner_recursion_pospos(accurate_lensing,ntheta,lmax,1,-1,d_ell_1m1)

        CALL compute_clg(ntheta,lmax,clpp,d_ell_11,clg)
        ! Change the sign of Clg,2 wrt Lewis et al.
        CALL compute_clg(ntheta,lmax,clpp,d_ell_1m1,clg2)

        do theta=0,ntheta - 1
            sigma_sq(theta) = clg(0) - clg(theta)
        enddo

        if (compute_tgradt .eqv. .True.) then
            CALL generate_d_wigner_recursion_pospos(accurate_lensing,ntheta,lmax,1+offset,-1,d_ell_1m1)
        endif

        CALL generate_d_wigner_recursion_pospos(accurate_lensing,ntheta,lmax,0+offset,0,d_ell_00)
        CALL generate_d_wigner_recursion_pospos(accurate_lensing,ntheta,lmax,2+offset,-2,d_ell_2m2)

        do l=0,lmax
            dl = real(l,kind=DP)
            array_clCMB(l) = (2*dl + 1.d0) / (4.d0*pi) * clCMB(l)
        enddo

        xi = 0.d0
        do theta=0,ntheta - 1
            do l=lmin,lmax
                dl = real(l,kind=DP)
                if (compute_tgradt .eqv. .True.) then
                    prefac1 = SQRT(dl * (dl + 1.d0)) * d_ell_00(theta,l) + dl*(dl + 1.d0)/2.d0 * clg2(theta) * &
                    & (SQRT( (dl + 2.d0) * (dl - 1.d0))*d_ell_1m1(theta,l) - SQRT(dl * (dl + 1.d0))*d_ell_00(theta,l))
                    prefac2 = SQRT(dl * (dl + 1.d0))  * (dl*(dl + 1.d0)/(4.d0) )**2 * d_ell_00(theta,l) +  &
                    & (dl + 2.d0)*(dl + 1.d0)*(dl)*(dl - 1.d0) / 16.d0 * &
                    & (SQRT( (dl + 3.d0) * (dl - 2.d0))*d_ell_2m2(theta,l) - SQRT( (dl + 2.d0) * (dl - 1.d0))*d_ell_1m1(theta,l)) ! Second order
                else
                    prefac1 = d_ell_00(theta,l) + dl*(dl + 1.d0)/2.d0 * clg2(theta) * d_ell_1m1(theta,l)
                    prefac2 = (dl*(dl + 1.d0)/(4.d0) )**2 * d_ell_00(theta,l) + &
                    & (dl + 2.d0)*(dl + 1.d0)*(dl)*(dl - 1.d0) / 16.d0 * d_ell_2m2(theta,l) ! Second order
                endif

                xi(theta) = xi(theta) + &
                 & array_clCMB(l) * exp(-dl*(dl + 1.d0)/2.d0*sigma_sq(theta)) * &
                 & (prefac1 + prefac2*clg2(theta)**2)
            enddo
        enddo

    end subroutine

    subroutine compute_correlation_function_pol(accurate_lensing,ntheta,clCMB,clpp,xim,xip,lmin,lmax)
        ! Compute the correlation function xi in polarization (EE and BB) for all angle
        ! see astro-ph/0601594, Eq. 9.20 and 9.21.
        implicit none
        integer, parameter       :: I4B = 4
        integer, parameter       :: DP = 8
        real(DP), parameter      :: pi = 3.14159265359

        integer, intent(in)      :: ntheta, lmax, lmin
        real(DP), intent(in)     :: clCMB(0:lmax), clpp(0:lmax)
        logical, intent(in)      :: accurate_lensing
        real(DP), intent(out)    :: xim(0:ntheta - 1),xip(0:ntheta - 1)
        real(DP)                 :: d_ell_22(0:ntheta - 1,0:lmax)
        real(DP)                 :: d_ell_31(0:ntheta - 1,0:lmax)
        real(DP)                 :: d_ell_3m3(0:ntheta - 1,0:lmax)
        real(DP)                 :: d_ell_11(0:ntheta - 1,0:lmax)
        real(DP)                 :: d_ell_2m2(0:ntheta - 1,0:lmax)
        real(DP)                 :: d_ell_1m1(0:ntheta - 1,0:lmax)
        real(DP)                 :: clg(0:ntheta - 1), clg2(0:ntheta - 1)
        real(DP)                 :: sigma_sq(0:ntheta - 1)
        real(DP)                 :: array_clCMB(0:lmax)
        real(DP)                 :: dl, prefac1m, prefac1p, dinfinity
        integer(I4B)             :: l, theta
        dinfinity = 1.e10_dp

        CALL generate_d_wigner_recursion_pospos(accurate_lensing,ntheta,lmax,2,2,d_ell_22)
        CALL generate_d_wigner_recursion_pospos(accurate_lensing,ntheta,lmax,3,1,d_ell_31)
        CALL generate_d_wigner_recursion_pospos(accurate_lensing,ntheta,lmax,3,-3,d_ell_3m3)
        CALL generate_d_wigner_recursion_pospos(accurate_lensing,ntheta,lmax,1,1,d_ell_11)
        CALL generate_d_wigner_recursion_pospos(accurate_lensing,ntheta,lmax,1,-1,d_ell_1m1)
        CALL generate_d_wigner_recursion_pospos(accurate_lensing,ntheta,lmax,2,-2,d_ell_2m2)

        CALL compute_clg(ntheta,lmax,clpp,d_ell_11,clg)
        ! Change the sign of Clg,2 wrt Lewis et al.
        CALL compute_clg(ntheta,lmax,clpp,d_ell_1m1,clg2)

        do theta=0,ntheta - 1
            sigma_sq(theta) = clg(0) - clg(theta)
        enddo

        do l=0,lmax
            dl = real(l,kind=DP)
            array_clCMB(l) = (2*dl + 1.d0) / (4.d0*pi) * clCMB(l)
        enddo

        xim = 0.d0
        xip = 0.d0
        do theta=0,ntheta - 1
            do l=lmin,lmax
                dl = real(l,kind=DP)
                prefac1p = d_ell_22(theta,l) + dl*(dl + 1.d0)/2.d0 * clg2(theta) * d_ell_31(theta,l)
                prefac1m = d_ell_2m2(theta,l) + dl*(dl + 1.d0)/4.d0 * clg2(theta) * &
                & (d_ell_1m1(theta,l) + d_ell_3m3(theta,l))

                xim(theta) = xim(theta) + &
                 & array_clCMB(l) * exp(-dl*(dl + 1.d0)/2.d0*sigma_sq(theta)) * prefac1m
                xip(theta) = xip(theta) + &
                  & array_clCMB(l) * exp(-dl*(dl + 1.d0)/2.d0*sigma_sq(theta)) * prefac1p
            enddo
        enddo

    end subroutine

    subroutine compute_correlation_function_temppol(accurate_lensing,ntheta,clCMB,clpp,xi,lmin,lmax)
        ! Compute the correlation function xi for TE for all angle
        ! see astro-ph/0601594, Eq. 9.22.
        implicit none
        integer, parameter       :: I4B = 4
        integer, parameter       :: DP = 8
        real(DP), parameter      :: pi = 3.14159265359

        integer, intent(in)      :: ntheta, lmax, lmin
        real(DP), intent(in)     :: clCMB(0:lmax), clpp(0:lmax)
        logical, intent(in)      :: accurate_lensing
        real(DP), intent(out)    :: xi(0:ntheta - 1)
        real(DP)                 :: d_ell_02(0:ntheta - 1,0:lmax)
        real(DP)                 :: d_ell_11(0:ntheta - 1,0:lmax)
        real(DP)                 :: d_ell_3m1(0:ntheta - 1,0:lmax)
        real(DP)                 :: d_ell_1m1(0:ntheta - 1,0:lmax)
        real(DP)                 :: clg(0:ntheta - 1), clg2(0:ntheta - 1)
        real(DP)                 :: sigma_sq(0:ntheta - 1)
        real(DP)                 :: array_clCMB(0:lmax)
        real(DP)                 :: dl, prefac1, dinfinity
        integer(I4B)             :: l, theta
        dinfinity = 1.e10_dp

        CALL generate_d_wigner_recursion_pospos(accurate_lensing,ntheta,lmax,0,2,d_ell_02)
        CALL generate_d_wigner_recursion_pospos(accurate_lensing,ntheta,lmax,1,1,d_ell_11)
        CALL generate_d_wigner_recursion_pospos(accurate_lensing,ntheta,lmax,1,-1,d_ell_1m1)
        CALL generate_d_wigner_recursion_pospos(accurate_lensing,ntheta,lmax,3,-1,d_ell_3m1)

        CALL compute_clg(ntheta,lmax,clpp,d_ell_11,clg)
        ! Change the sign of Clg,2 wrt Lewis et al.
        CALL compute_clg(ntheta,lmax,clpp,d_ell_1m1,clg2)

        do theta=0,ntheta - 1
            sigma_sq(theta) = clg(0) - clg(theta)
        enddo

        do l=0,lmax
            dl = real(l,kind=DP)
            array_clCMB(l) = (2*dl + 1.d0) / (4.d0*pi) * clCMB(l)
        enddo

        xi = 0.d0
        do theta=0,ntheta - 1
            do l=lmin,lmax
                dl = real(l,kind=DP)
                ! Super weird... I had to change the sign in order to get the result right
                ! I suspect the sign of clg2 to be -
                prefac1 = d_ell_02(theta,l) + dl*(dl + 1.d0)/4.d0 * clg2(theta) * &
                & (d_ell_11(theta,l) + d_ell_3m1(theta,l))

                xi(theta) = xi(theta) + &
                 & array_clCMB(l) * exp(-dl*(dl + 1.d0)/2.d0*sigma_sq(theta)) * prefac1
            enddo
        enddo

    end subroutine

    subroutine compute_clg(ntheta,lmax,clpp,d_ell_mm,clg)
        ! Compute C_l,g or C_l,g2 (see astro-ph/0601594)
        implicit none
        integer, parameter       :: I4B = 4
        integer, parameter       :: DP = 8
        real(DP), parameter      :: pi = 3.14159265359

        integer, intent(in)      :: ntheta, lmax
        real(DP), intent(in)     :: clpp(0:lmax)
        real(DP), intent(in)     :: d_ell_mm(0:ntheta - 1,0:lmax)
        real(DP), intent(out)    :: clg(0:ntheta - 1)
        integer(I4B)             :: l
        real(DP)                 :: dl
        real(DP)                 :: array_clpp(0:lmax)

        do l=2,lmax
            dl = real(l,kind=DP)
            array_clpp(l) = (2.d0*dl + 1.d0)*(dl + 1.d0)*dl / (4.d0*pi) * clpp(l)
        enddo
        array_clpp(0)=0.d0
        array_clpp(1)=0.d0

        clg = MATMUL(array_clpp,TRANSPOSE(d_ell_mm))

    end subroutine

    subroutine compute_der_correlation_function_temp(accurate_lensing,ntheta,clCMB,clpp,dxi,l2range,lmin,lmax_tmp,lmax,il2dim)
        ! Compute the derivative wrt phiphi of the correlation function xi in temperature for all angle
        implicit none
        integer, parameter       :: I4B = 4
        integer, parameter       :: DP = 8
        real(DP), parameter      :: pi = 3.14159265359

        integer, intent(in)      :: ntheta,lmin, lmax,lmax_tmp,il2dim
        integer, intent(in)      :: l2range(0:il2dim)
        real(DP), intent(in)     :: clCMB(0:lmax_tmp), clpp(0:lmax_tmp)
        logical, intent(in)      :: accurate_lensing
        real(DP), intent(out)    :: dxi(0:ntheta - 1,0:lmax)
        real(DP)                 :: d_ell_00(0:ntheta - 1,0:lmax_tmp)
        real(DP)                 :: d_ell_11(0:ntheta - 1,0:lmax_tmp)
        ! real(DP)                 :: d_ell_2m2(0:ntheta - 1,0:lmax)
        real(DP)                 :: d_ell_1m1(0:ntheta - 1,0:lmax_tmp)
        real(DP)                 :: clg(0:ntheta - 1), clg2(0:ntheta - 1)
        real(DP)                 :: sigma_sq(0:ntheta - 1),ll(0:lmax_tmp)
        real(DP)                 :: array_clCMB(0:lmax_tmp)
        real(DP)                 :: prefac1(0:ntheta - 1,0:lmax_tmp),prefacl2(0:ntheta - 1,0:lmax_tmp)
        real(DP)                 :: expo(0:ntheta - 1,0:lmax_tmp)
        real(DP)                 :: dl,dl2,dl3, prefac2(0:ntheta - 1),prefac2_scal
        integer(I4B)             :: l,l2,l3, theta,ltmp
        CHARACTER(LEN=13)        :: creturn

        CALL generate_d_wigner_recursion_pospos(accurate_lensing,ntheta,lmax_tmp,0,0,d_ell_00)
        CALL generate_d_wigner_recursion_pospos(accurate_lensing,ntheta,lmax_tmp,1,1,d_ell_11)
        CALL generate_d_wigner_recursion_pospos(accurate_lensing,ntheta,lmax_tmp,1,-1,d_ell_1m1)
        ! CALL generate_d_wigner_recursion_pospos(accurate_lensing,ntheta,lmax,2,-2,d_ell_2m2)

        CALL compute_clg(ntheta,lmax_tmp,clpp,d_ell_11,clg)
        ! Change the sign of Clg,2 wrt Lewis et al.
        CALL compute_clg(ntheta,lmax_tmp,clpp,d_ell_1m1,clg2)

        do theta=0,ntheta - 1
            sigma_sq(theta) = clg(0) - clg(theta)
        enddo

        do l=lmin,lmax_tmp
            dl = real(l,kind=DP)
            array_clCMB(l) = (2*dl + 1.d0)* (dl + 1.d0) * dl / (8.d0*pi) * clCMB(l)
            ll(l) = (2*dl + 1.d0)* (dl + 1.d0) * dl / (4.d0*pi)
            do theta=0,ntheta - 1
                ! prefac1(theta,l) = (d_ell_00(theta,l) + dl*(dl + 1.d0)/2.d0 * clg2(theta) * d_ell_1m1(theta,l))* &
                ! & exp(-dl*(dl + 1.d0)/2.d0*sigma_sq(theta))
                prefac1(theta,l) = (d_ell_00(theta,l) + dl*(dl + 1.d0)/2.d0 * clg2(theta) * d_ell_1m1(theta,l))
                expo(theta,l) = exp(-dl*(dl + 1.d0)/2.d0*sigma_sq(theta))
                prefacl2(theta,l) = -(d_ell_11(0,l) - d_ell_11(theta,l))
            enddo
        enddo

        ! dxi = 0.d0
        ! do theta=0,ntheta - 1
        !     creturn = achar(13)
        !     WRITE( * , 101 , ADVANCE='NO' ) creturn , theta , ntheta
        !     101     FORMAT( a , 'ell : ',i7,' out of a total of ',i7)
        !     do l2=2,lmax
        !         do l3=2,lmax
        !             prefac2_scal = d_ell_1m1(theta,l2)*d_ell_1m1(theta,l3)
        !
        !             dxi(theta,l2) = dxi(theta,l2) + &
        !              & ll(l2) * array_clCMB(l3)  * &
        !              & (prefacl2(theta,l2)*prefac1(theta,l3) + prefac2_scal)
        !         enddo
        !     enddo
        ! enddo

        dxi = 0.d0
        do ltmp=0,il2dim
            l2=l2range(ltmp)
            creturn = achar(13)
            WRITE( * , 101 , ADVANCE='NO' ) creturn , l2 , lmax
            101     FORMAT( a , 'ell : ',i7,' out of a total of ',i7)
            do l3=lmin,lmax_tmp
                prefac2(0:ntheta - 1) = d_ell_1m1(0:ntheta - 1,l2)*d_ell_1m1(0:ntheta - 1,l3)

                dxi(0:ntheta - 1,l2) = dxi(0:ntheta - 1,l2) + &
                 & ll(l2) * array_clCMB(l3) * expo(0:ntheta - 1,l3) * &
                 & (prefacl2(0:ntheta - 1,l2)*prefac1(0:ntheta - 1,l3) + prefac2)
            enddo
        enddo

    end subroutine

    subroutine compute_der_correlation_function_pol(accurate_lensing,ntheta,clCMB,clpp,dxim,&
        & dxip,l2range,lmin,lmax_tmp,lmax,il2dim)
        ! Compute the derivative wrt phiphi of the correlation functions xim and xip in polarization for all angle
        implicit none
        integer, parameter       :: I4B = 4
        integer, parameter       :: DP = 8
        real(DP), parameter      :: pi = 3.14159265359

        integer, intent(in)      :: ntheta,lmin, lmax,lmax_tmp,il2dim
        integer, intent(in)      :: l2range(0:il2dim)
        real(DP), intent(in)     :: clCMB(0:lmax_tmp), clpp(0:lmax_tmp)
        logical, intent(in)      :: accurate_lensing
        real(DP), intent(out)    :: dxip(0:ntheta - 1,0:lmax),dxim(0:ntheta - 1,0:lmax)
        real(DP)                 :: d_ell_11(0:ntheta - 1,0:lmax_tmp)
        real(DP)                 :: d_ell_31(0:ntheta - 1,0:lmax_tmp)
        real(DP)                 :: d_ell_3m3(0:ntheta - 1,0:lmax_tmp)
        real(DP)                 :: d_ell_2m2(0:ntheta - 1,0:lmax_tmp)
        real(DP)                 :: d_ell_22(0:ntheta - 1,0:lmax_tmp)
        real(DP)                 :: d_ell_1m1(0:ntheta - 1,0:lmax_tmp)
        real(DP)                 :: clg(0:ntheta - 1), clg2(0:ntheta - 1)
        real(DP)                 :: sigma_sq(0:ntheta - 1),ll(0:lmax_tmp)
        real(DP)                 :: array_clCMB(0:lmax_tmp)
        real(DP)                 :: expo(0:ntheta - 1,0:lmax_tmp)
        real(DP)                 :: prefac1m(0:ntheta - 1,0:lmax_tmp),prefac2m(0:ntheta - 1)
        real(DP)                 :: prefac1p(0:ntheta - 1,0:lmax_tmp),prefac2p(0:ntheta - 1)
        real(DP)                 :: dl,dl2,dl3, prefacl2(0:ntheta - 1,0:lmax_tmp)
        integer(I4B)             :: l,l2,l3, theta,ltmp
        CHARACTER(LEN=13)        :: creturn

        CALL generate_d_wigner_recursion_pospos(accurate_lensing,ntheta,lmax_tmp,2,2,d_ell_22)
        CALL generate_d_wigner_recursion_pospos(accurate_lensing,ntheta,lmax_tmp,3,1,d_ell_31)
        CALL generate_d_wigner_recursion_pospos(accurate_lensing,ntheta,lmax_tmp,3,-3,d_ell_3m3)
        CALL generate_d_wigner_recursion_pospos(accurate_lensing,ntheta,lmax_tmp,1,1,d_ell_11)
        CALL generate_d_wigner_recursion_pospos(accurate_lensing,ntheta,lmax_tmp,1,-1,d_ell_1m1)
        CALL generate_d_wigner_recursion_pospos(accurate_lensing,ntheta,lmax_tmp,2,-2,d_ell_2m2)

        CALL compute_clg(ntheta,lmax_tmp,clpp,d_ell_11,clg)
        ! Change the sign of Clg,2 wrt Lewis et al.
        CALL compute_clg(ntheta,lmax_tmp,clpp,d_ell_1m1,clg2)

        do theta=0,ntheta - 1
            sigma_sq(theta) = clg(0) - clg(theta)
        enddo

        do l=lmin,lmax_tmp
            dl = real(l,kind=DP)
            array_clCMB(l) = (2*dl + 1.d0)* (dl + 1.d0) * dl / (8.d0*pi) * clCMB(l)
            ll(l) = (2*dl + 1.d0)* (dl + 1.d0) * dl / (4.d0*pi)
            do theta=0,ntheta - 1
                prefac1p(theta,l) = (d_ell_22(theta,l) + dl*(dl + 1.d0)/2.d0 * &
                & clg2(theta) * d_ell_31(theta,l))

                prefac1m(theta,l) = (d_ell_2m2(theta,l) + dl*(dl + 1.d0)/4.d0 * clg2(theta) * &
                & (d_ell_1m1(theta,l) + d_ell_3m3(theta,l)))

                prefacl2(theta,l) = -(d_ell_11(0,l) - d_ell_11(theta,l))

                expo(theta,l) = exp(-dl*(dl + 1.d0)/2.d0*sigma_sq(theta))
            enddo
        enddo

        dxip = 0.d0
        dxim = 0.d0
        do ltmp=0,il2dim
            l2=l2range(ltmp)
            creturn = achar(13)
            WRITE( * , 101 , ADVANCE='NO' ) creturn , l2 , lmax
            101     FORMAT( a , 'ell : ',i7,' out of a total of ',i7)
            do l3=lmin,lmax_tmp
                prefac2p(0:ntheta - 1) = d_ell_1m1(0:ntheta - 1,l2)*d_ell_31(0:ntheta - 1,l3)
                prefac2m(0:ntheta - 1) = d_ell_1m1(0:ntheta - 1,l2) * &
                & (d_ell_1m1(0:ntheta - 1,l3) + d_ell_3m3(0:ntheta - 1,l3))/2.d0

                dxip(0:ntheta - 1,l2) = dxip(0:ntheta - 1,l2) + &
                 & ll(l2) * array_clCMB(l3) * expo(0:ntheta - 1,l3) * &
                 & (prefacl2(0:ntheta - 1,l2)*prefac1p(0:ntheta - 1,l3) + prefac2p)

                dxim(0:ntheta - 1,l2) = dxim(0:ntheta - 1,l2) + &
                 & ll(l2) * array_clCMB(l3) * expo(0:ntheta - 1,l3) * &
                 & (prefacl2(0:ntheta - 1,l2)*prefac1m(0:ntheta - 1,l3) + prefac2m)
            enddo
        enddo

    end subroutine

    subroutine compute_der_correlation_function_temppol(accurate_lensing,ntheta,clCMB,clpp,dxim,l2range,lmin,lmax_tmp,lmax,il2dim)
        ! Compute the derivative wrt phiphi of the correlation functions xim and xip in polarization for all angle
        implicit none
        integer, parameter       :: I4B = 4
        integer, parameter       :: DP = 8
        real(DP), parameter      :: pi = 3.14159265359

        integer, intent(in)      :: ntheta,lmin, lmax,lmax_tmp,il2dim
        integer, intent(in)      :: l2range(0:il2dim)
        real(DP), intent(in)     :: clCMB(0:lmax_tmp), clpp(0:lmax_tmp)
        logical, intent(in)      :: accurate_lensing
        real(DP), intent(out)    :: dxim(0:ntheta - 1,0:lmax)
        real(DP)                 :: d_ell_02(0:ntheta - 1,0:lmax_tmp)
        real(DP)                 :: d_ell_11(0:ntheta - 1,0:lmax_tmp)
        real(DP)                 :: d_ell_3m1(0:ntheta - 1,0:lmax_tmp)
        real(DP)                 :: d_ell_1m1(0:ntheta - 1,0:lmax_tmp)
        real(DP)                 :: clg(0:ntheta - 1), clg2(0:ntheta - 1)
        real(DP)                 :: sigma_sq(0:ntheta - 1),ll(0:lmax_tmp)
        real(DP)                 :: array_clCMB(0:lmax_tmp)
        real(DP)                 :: expo(0:ntheta - 1,0:lmax_tmp)
        real(DP)                 :: prefac1(0:ntheta - 1,0:lmax_tmp),prefac2(0:ntheta - 1)
        real(DP)                 :: dl,dl2,dl3, prefacl2(0:ntheta - 1,0:lmax_tmp)
        integer(I4B)             :: l,l2,l3, theta,ltmp
        CHARACTER(LEN=13)        :: creturn

        CALL generate_d_wigner_recursion_pospos(accurate_lensing,ntheta,lmax_tmp,0,2,d_ell_02)
        CALL generate_d_wigner_recursion_pospos(accurate_lensing,ntheta,lmax_tmp,1,1,d_ell_11)
        CALL generate_d_wigner_recursion_pospos(accurate_lensing,ntheta,lmax_tmp,1,-1,d_ell_1m1)
        CALL generate_d_wigner_recursion_pospos(accurate_lensing,ntheta,lmax_tmp,3,-1,d_ell_3m1)

        CALL compute_clg(ntheta,lmax_tmp,clpp,d_ell_11,clg)
        ! Change the sign of Clg,2 wrt Lewis et al.
        CALL compute_clg(ntheta,lmax_tmp,clpp,d_ell_1m1,clg2)

        do theta=0,ntheta - 1
            sigma_sq(theta) = clg(0) - clg(theta)
        enddo

        do l=lmin,lmax_tmp
            dl = real(l,kind=DP)
            array_clCMB(l) = (2*dl + 1.d0)* (dl + 1.d0) * dl / (8.d0*pi) * clCMB(l)
            ll(l) = (2*dl + 1.d0)* (dl + 1.d0) * dl / (4.d0*pi)
            do theta=0,ntheta - 1

                prefac1(theta,l) = (d_ell_02(theta,l) + dl*(dl + 1.d0)/4.d0 * clg2(theta) * &
                & (d_ell_11(theta,l) + d_ell_3m1(theta,l)))

                prefacl2(theta,l) = -(d_ell_11(0,l) - d_ell_11(theta,l))

                expo(theta,l) = exp(-dl*(dl + 1.d0)/2.d0*sigma_sq(theta))
            enddo
        enddo

        dxim = 0.d0
        do ltmp=0,il2dim
            l2=l2range(ltmp)
            creturn = achar(13)
            WRITE( * , 101 , ADVANCE='NO' ) creturn , l2 , lmax
            101     FORMAT( a , 'ell : ',i7,' out of a total of ',i7)
            do l3=lmin,lmax_tmp
                prefac2(0:ntheta - 1) = d_ell_1m1(0:ntheta - 1,l2) * &
                & (d_ell_11(0:ntheta - 1,l3) + d_ell_3m1(0:ntheta - 1,l3))/2.d0

                dxim(0:ntheta - 1,l2) = dxim(0:ntheta - 1,l2) + &
                 & ll(l2) * array_clCMB(l3) * expo(0:ntheta - 1,l3) * &
                 & (prefacl2(0:ntheta - 1,l2)*prefac1(0:ntheta - 1,l3) + prefac2)
            enddo
        enddo

    end subroutine

    subroutine lensed_spectra_corrfunc_allsky(ntheta,flavor,accurate_lensing,compute_tgradt,clCMB,clpp,clCMB_lensed,lmin,lmax)
        ! Generate lensed spectra using correlation function
        ! Main equations are coming from astro-ph/0601594
        ! The core of the computation lies in a fast computation of Wigner d-matrices using recursion formula.
        implicit none
        integer, parameter       :: I4B = 4
        integer, parameter       :: DP = 8
        real(DP), parameter      :: pi = 3.14159265359

        integer, intent(in)      :: ntheta, lmax, lmin
        real(DP), intent(in)     :: clCMB(0:lmax), clpp(0:lmax)
        logical, intent(in)      :: accurate_lensing, compute_tgradt
        CHARACTER(LEN=13), intent(in) :: flavor
        real(DP), intent(out)    :: clCMB_lensed(0:lmax)

        integer(I4B)             :: l, theta, offset
        real(DP)                 :: dtheta, dl
        real(DP)                 :: fac(0:lmax)
        real(DP)                 :: cosvec(0:ntheta - 1), sinvec(0:ntheta - 1)
        real(DP)                 :: d_ell_00(0:ntheta - 1,0:lmax)
        real(DP)                 :: d_ell_2m2(0:ntheta - 1,0:lmax)
        real(DP)                 :: d_ell_22(0:ntheta - 1,0:lmax)
        real(DP)                 :: d_ell_20(0:ntheta - 1,0:lmax)
        real(DP)                 :: xi(0:ntheta - 1),xim(0:ntheta - 1),xip(0:ntheta - 1)
        real(DP)                 :: xabsc(0:ntheta - 1)
        real(DP)                 :: weig(0:ntheta - 1), weig_tmp(0:ntheta - 1)

        if (accurate_lensing .eqv. .True.) then
            print *,'Accurate lensing'
            CALL gauleg(ntheta, xabsc, weig)
            do theta=0,ntheta - 1
                weig_tmp(theta) = weig(ntheta - 1 - theta)
            enddo
            weig = weig_tmp
        else
            print *,'Simple integration'
            CALL generate_cosine_sine_vec(pi,ntheta,cosvec(0:ntheta - 1),sinvec(0:ntheta - 1))
            dtheta = pi / ntheta
            do theta=0,ntheta - 1
                weig(theta) = sinvec(theta)*dtheta
            enddo
        endif

        if (compute_tgradt .eqv. .True.) then
            print *,'Computing spectrum of T.gradT'
            offset = 1
            do l=0,lmax
                dl = real(l,kind=DP)
                fac(l) = SQRT(dl * (dl + 1.d0))
            enddo
        else
            offset = 0
            do l=0,lmax
                dl = real(l,kind=DP)
                fac(l) = 1.d0
            enddo
        endif

        if (flavor .eq. 'cltt') then
            CALL compute_correlation_function_temp(accurate_lensing,compute_tgradt,ntheta,clCMB,clpp,xi,lmin,lmax)
            CALL generate_d_wigner_recursion_pospos(accurate_lensing,ntheta,lmax,0+offset,0,d_ell_00)
        elseif (flavor .eq. 'clee' .or. flavor .eq. 'clbb') then
            CALL compute_correlation_function_pol(accurate_lensing,ntheta,clCMB,clpp,xim,xip,lmin,lmax)
            CALL generate_d_wigner_recursion_pospos(accurate_lensing,ntheta,lmax,2,-2,d_ell_2m2)
            CALL generate_d_wigner_recursion_pospos(accurate_lensing,ntheta,lmax,2,2,d_ell_22)
        elseif (flavor .eq. 'clte') then
            CALL compute_correlation_function_temppol(accurate_lensing,ntheta,clCMB,clpp,xi,lmin,lmax)
            CALL generate_d_wigner_recursion_pospos(accurate_lensing,ntheta,lmax,2,0,d_ell_20)
        endif

        do l=lmin,lmax
            dl = real(l,kind=DP)
            clCMB_lensed(l) = 0.d0
            do theta=0,ntheta - 1
                if (flavor .eq. 'cltt') then
                    clCMB_lensed(l) = clCMB_lensed(l) + &
                    & (2*pi) * xi(theta) * d_ell_00(theta,l) * weig(theta) / fac(l)
                elseif (flavor .eq. 'clee') then
                    clCMB_lensed(l) = clCMB_lensed(l) + &
                    & (2*pi) * 0.5*(xim(theta)*d_ell_2m2(theta,l) + xip(theta)*d_ell_22(theta,l))&
                    & * weig(theta)
                elseif (flavor .eq. 'clbb') then
                    clCMB_lensed(l) = clCMB_lensed(l) + &
                    & (2*pi) * 0.5*(-xim(theta)*d_ell_2m2(theta,l) + xip(theta)*d_ell_22(theta,l))&
                    & * weig(theta)
                elseif (flavor .eq. 'clte') then
                    clCMB_lensed(l) = clCMB_lensed(l) + &
                    & (2*pi) * xi(theta) * d_ell_20(theta,l) * weig(theta)
                else
                    print *,'Not yet implemented',flavor
                endif
            enddo
        enddo

    end subroutine

    subroutine derivative_dclcmbdclpp_corr_func(ntheta,flavor,accurate_lensing,clCMB, &
        & clpp,dxim,dxip,dellm,dellp,l2range,lmin,il2dim,lmax,lmax_tmp)
        ! Main function to compute derivatives wrt phiphi of the correlation functions (it calls
        ! other functions above)
        implicit none
        integer, parameter       :: I4B = 4
        integer, parameter       :: DP = 8
        real(DP), parameter      :: pi = 3.14159265359

        integer, intent(in)      :: ntheta, lmax, il2dim, lmin, lmax_tmp
        integer, intent(in)      :: l2range(0:il2dim)
        real(DP), intent(in)     :: clCMB(0:lmax_tmp), clpp(0:lmax_tmp)
        logical, intent(in)      :: accurate_lensing
        CHARACTER(LEN=13), intent(in) :: flavor
        real(DP),intent(out)     :: dxim(0:ntheta - 1,0:lmax),dxip(0:ntheta - 1,0:lmax)
        real(DP),intent(out)     :: dellm(0:ntheta - 1,0:lmax),dellp(0:ntheta - 1,0:lmax)

        integer(I4B)             :: l1,l2, theta
        real(DP)                 :: dtheta
        real(DP)                 :: cosvec(0:ntheta - 1), sinvec(0:ntheta - 1)
        real(DP)                 :: d_ell_00(0:ntheta - 1,0:lmax)
        real(DP)                 :: d_ell_2m2(0:ntheta - 1,0:lmax)
        real(DP)                 :: d_ell_22(0:ntheta - 1,0:lmax)
        real(DP)                 :: d_ell_20(0:ntheta - 1,0:lmax)
        real(DP)                 :: dxi(0:ntheta - 1,0:lmax)
        real(DP)                 :: xabsc(0:ntheta - 1)
        real(DP)                 :: weig(0:ntheta - 1), weig_tmp(0:ntheta - 1)
        CHARACTER(LEN=13)        :: creturn

        if (accurate_lensing .eqv. .True.) then
            print *,'Accurate lensing'
            CALL gauleg(ntheta, xabsc, weig)
            do theta=0,ntheta - 1
                weig_tmp(theta) = weig(ntheta - 1 - theta)
            enddo
            weig = weig_tmp
        else
            print *,'Simple integration'
            CALL generate_cosine_sine_vec(pi,ntheta,cosvec(0:ntheta - 1),sinvec(0:ntheta - 1))
            dtheta = pi / ntheta
            do theta=0,ntheta - 1
                weig(theta) = sinvec(theta)*dtheta
            enddo
        endif

        if (flavor .eq. 'cltt') then
            CALL compute_der_correlation_function_temp(accurate_lensing,ntheta,clCMB,clpp,dxim,&
            & l2range,lmin,lmax_tmp,lmax,il2dim)
            CALL generate_d_wigner_recursion_pospos(accurate_lensing,ntheta,lmax,0,0,d_ell_00)
        elseif (flavor .eq. 'clee' .or. flavor .eq. 'clbb') then
            CALL compute_der_correlation_function_pol(accurate_lensing,ntheta,clCMB,clpp,dxim,dxip,&
            & l2range,lmin,lmax_tmp,lmax,il2dim)
            CALL generate_d_wigner_recursion_pospos(accurate_lensing,ntheta,lmax,2,-2,d_ell_2m2)
            CALL generate_d_wigner_recursion_pospos(accurate_lensing,ntheta,lmax,2,2,d_ell_22)
        elseif (flavor .eq. 'clte') then
            CALL compute_der_correlation_function_temppol(accurate_lensing,ntheta,clCMB,clpp,dxim,&
            & l2range,lmin,lmax_tmp,lmax,il2dim)
            CALL generate_d_wigner_recursion_pospos(accurate_lensing,ntheta,lmax,2,0,d_ell_20)
        endif


        if (flavor .eq. 'cltt') then
            do theta=0,ntheta - 1
                dxim(theta,0:lmax) = dxim(theta,0:lmax) * weig(theta)
            enddo
            dellm = d_ell_00
        elseif (flavor .eq. 'clee') then
            do theta=0,ntheta - 1
                dxip(theta,0:lmax) = dxip(theta,0:lmax) * weig(theta)
                dxim(theta,0:lmax) = dxim(theta,0:lmax) * weig(theta)
            enddo
            dellm = d_ell_2m2
            dellp = d_ell_22
        elseif (flavor .eq. 'clbb') then
            do theta=0,ntheta - 1
                dxip(theta,0:lmax) = dxip(theta,0:lmax) * weig(theta)
                dxim(theta,0:lmax) = dxim(theta,0:lmax) * weig(theta)
            enddo
            dellm = d_ell_2m2
            dellp = d_ell_22
        elseif (flavor .eq. 'clte') then
            do theta=0,ntheta - 1
                dxim(theta,0:lmax) = dxim(theta,0:lmax) * weig(theta)
            enddo
            dellm = d_ell_20
        else
            print *,'Not yet implemented ',flavor
        endif

        do l1=0,lmin
            dellm(0:ntheta - 1,l1) = 0.0
            dellp(0:ntheta - 1,l1) = 0.0
        enddo

    end subroutine

end module
