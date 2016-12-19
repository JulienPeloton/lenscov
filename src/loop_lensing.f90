!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright (C) 2016 Peloton
! 10/2015
! FORTRAN routines to compute low level products
! Main purposes is interfacing with python (using f2py)
! See 1611.01446
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module loop_lensing

contains

    subroutine f_functions(name,dclCMB,il1,il2,il3,dW_spin0,dW_spin2,df_l2l1l3,ilmax)
        implicit none

        integer, parameter       :: I4B = 4
        integer, parameter       :: DP = 8
        real(DP), parameter      :: pi = 3.141592

        CHARACTER(LEN=13), intent(in) :: name
        integer(I4B), intent(in) :: il1,il2,il3,ilmax
        real(DP), intent(in)     :: dW_spin0(0:2*ilmax), dW_spin2(0:2*ilmax),dclCMB(0:3,0:2*ilmax)
        real(DP)                 :: dl1,dl2,dl3,dL_l2l1l3,dL_l3l1l2,dS_l1l2l3,dinfinity
        real(DP), intent(out)    :: df_l2l1l3
        integer(I4B)             :: iL

        dl1 = real(il1, kind=DP)
        dl2 = real(il2, kind=DP)
        dl3 = real(il3, kind=DP)

        dinfinity = 1.e20_dp

        ! Compute prefactors
        dL_l2l1l3 = dl1*(dl1+1.d0) - dl2*(dl2+1.d0) + dl3*(dl3+1.d0)
        dL_l3l1l2 = dl1*(dl1+1.d0) + dl2*(dl2+1.d0) - dl3*(dl3+1.d0)
        dS_l1l2l3 = SQRT( ( (2*dl1+1.d0)*(2*dl2+1.d0)*(2*dl3+1.d0) )/(16*pi) )

        iL = il1 + il2 + il3

        ! Compute function f
        if ( name == 'TT' .and. MOD(iL,2) .eq. 0) then
            df_l2l1l3 = (dclCMB(0,il2) * dL_l3l1l2 + dclCMB(0,il3) * dL_l2l1l3) * dS_l1l2l3 * dW_spin0(il1)
        elseif ( name == 'EE' .and. MOD(iL,2) .eq. 0) then
            df_l2l1l3 = (dclCMB(1,il2) * dL_l3l1l2 + dclCMB(1,il3) * dL_l2l1l3) * dS_l1l2l3 * dW_spin2(il1)
        elseif ( name == 'TE' .and. MOD(iL,2) .eq. 0) then
            df_l2l1l3 = (dclCMB(3,il2) * dL_l3l1l2 * dW_spin2(il1) + dclCMB(3,il3) * dL_l2l1l3 * dW_spin0(il1)) * dS_l1l2l3
        elseif ( name == 'BB' ) then
            df_l2l1l3 = 0.d0
        elseif ( name == 'EB' ) then
            df_l2l1l3 = 0.d0
        elseif ( name == 'TB' ) then
            df_l2l1l3 = 0.d0
        else
            df_l2l1l3 = 0.d0
        endif

        ! Check if weird things happened (numerical instabilities or infinities)
        if (df_l2l1l3 .ne. df_l2l1l3 .or. ABS(df_l2l1l3) .gt. dinfinity) then
            df_l2l1l3 = 0.0_dp
        endif

    end subroutine

    subroutine covariance_cmbxcmb_order1_uvupvp(dclCMB,dclpp,dsigmal2l3,il2range,&
        & uup,vvp,uvp,vup,ilmin,ispin0l2,ispin0l3,ispin2l2,ispin2l3,il2dim,ilmax)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! *Description:
        !   This function computes the O(C^{phi phi}) of the covariance of lensed CMB spectra.
    	!   This corresponds to the second term of the RHS in Eq. 52 of 1308.0286
        !   The routine WIG3J compute all allowed l1 given fixed l2,l3,m2,m3.
        !   Symmetries are already taken into account (even l1+l2+l3).
        !   Variables starting with /i/ are int, while variables starting with /d/ are double
        ! *Arguments:
        !   * dclCMB:   IN       unlensed CMB power spectrum (double, (3,1)D array)
        !   * dclpp:    IN       lensing potential power spectrum (double, 1D array)
        !   * il2range: IN       Values of l2 (1D array)
        !   * flavor:   IN       Name of the spectrum (cltt, clee, etc.)
        !   * ilmin:    IN       Minimum ell (int)
        !   * ispin2:   IN       spin value associated to l2 (i.e. m2) (int)
        !   * ispin3:   IN       spin value associated to l3 (i.e. m3) (int)
        !   * il2dim:   IN       Length of the il2range (int)
        !   * ilmax:    IN       Length of dclCMB and dclpp (int)
        !   * dsigmal2l3: OUT    Covariance matrix (double, 2D array)
        ! *Author: Julien Peloton, University of Sussex
        ! *Revision history (mmyy):
        !   * 1015: creation
        !   * 1015: add polarization features
        ! *todo:
        !   * Add explicitely MPI inside the routine (instead of playing with il2range)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        implicit none

        integer, parameter       :: I4B = 4
        integer, parameter       :: DP = 8
        real(DP), parameter      :: pi = 3.141592

        integer(I4B), intent(in) :: ilmin,ilmax,ispin0l2,ispin0l3,ispin2l2,ispin2l3,il2dim
        integer(I4B), intent(in) :: il2range(0:il2dim)
        real(DP), intent(in)     :: dclCMB(0:3,0:2*ilmax), dclpp(0:2*ilmax)
        CHARACTER(LEN=13), intent(in) :: uup,vvp,uvp,vup
        real(DP), intent(out)    :: dsigmal2l3(0:ilmax,0:ilmax)

        integer(I4B)             :: il1, il1min, il1max, il2, il3, iL, indim, il_tmp
        real(DP)                 :: dl1, dl1min, dl1max, dl2, dl3
        real(DP)                 :: dinfinity, dier, dspin0l2, dspin0l3,dspin2l2, dspin2l3
        real(DP)                 :: fuup,fvvp,fuvp,fvup
        CHARACTER(LEN=13)        :: creturn

        real(DP), allocatable    :: dW_l1l2l3_spin0(:), dW_l1l2l3_spin2(:)

        ! Initialize
        dinfinity = 1.e20_dp
        ALLOCATE(dW_l1l2l3_spin0(0:2*ilmax))
        ALLOCATE(dW_l1l2l3_spin2(0:2*ilmax))
        dW_l1l2l3_spin0 = 0.0_dp
        dW_l1l2l3_spin2 = 0.0_dp
        dier = 0.0_dp

        dspin0l2 = real(ispin0l2, kind=DP)
        dspin0l3 = real(ispin0l3, kind=DP)
        dspin2l2 = real(ispin2l2, kind=DP)
        dspin2l3 = real(ispin2l3, kind=DP)

        do il_tmp=0, il2dim
            il2 = il2range(il_tmp)
            creturn = achar(13)
            WRITE( * , 101 , ADVANCE='NO' ) creturn , il2 , ilmax
            101     FORMAT( a , 'ell : ',i7,' out of a total of ',i7)
            dl2 = real(il2, kind=DP)
            do il3=il2, ilmax
                dl3 = real(il3, kind=DP)

                ! Define boundaries for l1 given l2 and l3
                il1min = MAX( abs(il2-il3), abs(ispin0l2+ispin0l3) )
                il1max = il2 + il3
                dl1min = real(il1min, kind=DP)
                dl1max = real(il1max, kind=DP)

                ! Number of allowed l1 given l2 and l3
                indim = il1max - il1min + 1

                ! The fortran routine WIG3J compute all allowed l1 given fixed l2,l3,m2,m3.
                ! Therefore we loop over l2,l3 (m2,m3 are fixed and given by the spins), and at each iteration we grab a
                ! vector dW_l1l2l3 of dimension indim, which contains all 3-j symbols, that we will use to update the final Cl
                ! Be careful for the type of variables to pass. See wig3j_f.f for more infos.
                call WIG3J(dl2, dl3, dspin0l2, dspin0l3, dl1min, dl1max, &
                & dW_l1l2l3_spin0(il1min:il1max), indim, dier)
                call WIG3J(dl2, dl3, dspin2l2, dspin2l3, dl1min, dl1max, &
                & dW_l1l2l3_spin2(il1min:il1max), indim, dier)

                ! Define boundaries for l1 given l2 and l3 AND lmax
                ! Comment for more accuracy, but slower.
                il1max = MIN( il2 + il3 , ilmax)

                do il1=il1min, il1max
                    if (il1 .lt. 2) CYCLE
                    dl1 = real(il1, kind=DP)

                    ! Compute function f
                    call f_functions(uup,dclCMB,il1,il2,il3,dW_l1l2l3_spin0,dW_l1l2l3_spin2,fuup,ilmax)
                    call f_functions(vvp,dclCMB,il1,il2,il3,dW_l1l2l3_spin0,dW_l1l2l3_spin2,fvvp,ilmax)
                    call f_functions(uvp,dclCMB,il1,il2,il3,dW_l1l2l3_spin0,dW_l1l2l3_spin2,fuvp,ilmax)
                    call f_functions(vup,dclCMB,il1,il2,il3,dW_l1l2l3_spin0,dW_l1l2l3_spin2,fvup,ilmax)

                    ! Update the Cl
                    dsigmal2l3(il2,il3) = dsigmal2l3(il2,il3) + &
                    & dclpp(il1) / (2.d0*dl2+1.d0) / (2.d0*dl3+1.d0) * (fuup*fvvp + fuvp*fvup)
                enddo
                dsigmal2l3(il3,il2) = dsigmal2l3(il2,il3)
            enddo
        enddo

        DEALLOCATE(dW_l1l2l3_spin0)
        DEALLOCATE(dW_l1l2l3_spin2)

    end subroutine

    subroutine compute_derivatives_dcttdcpp_mpi(dclCMB,dDcl2TTDcl1pp,il2range,flavor,ilmin &
        & ,ispinl2_x,ispinl3_x,ispinl2_y,ispinl3_y,il2dim,ilmax)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! *Description:
        !   This function computes the derivative of C_\ell^{XX} with respect to C_\ell^{phiphi} for X=T,E.
        !   We use the analytical expression in Eq. 27 of of 1308.0286.
        !   The polarisation version is trivially added.
        !   The routine WIG3J compute all allowed l1 given fixed l2,l3,m2,m3.
        !   Symmetries are already taken into account (even l1+l2+l3).
        !   Variables starting with /i/ are int, while variables starting with /d/ are double
        ! *Arguments:
        !   * dclCMB:   IN       unlensed CMB power spectrum (double, 1D array)
        !   * il2range: IN       Values of l2 (double, 1D array)
        !   * flavor:   IN       Name of the spectrum (character, cltt, clee, etc.)
        !   * ilmin:    IN       Minimum ell (int)
        !   * ispin2:   IN       spin value associated to l2 (i.e. m2) (int)
        !   * ispin3:   IN       spin value associated to l3 (i.e. m3) (int)
        !   * il2dim:   IN       Length of the il2range (int)
        !   * ilmax:    IN       Length of dclCMB and dclpp (int)
        !   * dDcl2TTDcl1pp: OUT Matrix containing the derivatives (double, 2D array)
        ! *Author: Julien Peloton, University of Sussex
        ! *Revision history (mmyy):
        !   * 1015: creation
        !   * 1015: add polarization features
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        implicit none

        integer, parameter       :: I4B = 4
        integer, parameter       :: DP = 8
        real(DP), parameter      :: pi = 3.141592

        integer(I4B), intent(in) :: ilmin,ispinl2_x,ispinl3_x,ispinl2_y,ispinl3_y,ilmax,il2dim
        integer(I4B), intent(in) :: il2range(0:il2dim)
        real(DP), intent(in)     :: dclCMB(0:2*ilmax)
        CHARACTER(LEN=13), intent(in) :: flavor
        real(DP), intent(out)    :: dDcl2TTDcl1pp(0:ilmax,0:ilmax)

        integer(I4B)             :: il1, il1min, il1max, il2, il3, iL, indim, il_tmp
        real(DP)                 :: dl1, dl1min, dl1max, dl2, dl3, dZl2l3
        real(DP)                 :: dL_l2l3l1, dS_l1l2l3, dF_l2l3l1, dinfinity, dier
        real(DP)                 :: dspinl2_x, dspinl3_x, dspinl2_y, dspinl3_y
        CHARACTER(LEN=13)        :: creturn

        real(DP), allocatable    :: dW_l1l2l3_x(:), dW_l1l2l3_y(:)

        ! Initialize
        dinfinity = 1.e20_dp
        ALLOCATE(dW_l1l2l3_x(0:2*ilmax))
        ALLOCATE(dW_l1l2l3_y(0:2*ilmax))
        dW_l1l2l3_x = 0.0_dp
        dW_l1l2l3_y = 0.0_dp
        dier = 0.0_dp

        dspinl2_x = real(ispinl2_x, kind=DP)
        dspinl3_x = real(ispinl3_x, kind=DP)
        dspinl2_y = real(ispinl2_y, kind=DP)
        dspinl3_y = real(ispinl3_y, kind=DP)

        !!!! PRECOMPUTE DERIVATIVE !!!!
        do il_tmp=0, il2dim
            il2 = il2range(il_tmp)
            creturn = achar(13)
            WRITE( * , 101 , ADVANCE='NO' ) creturn , il2 , ilmax
            101     FORMAT( a , 'ell : ',i7,' out of a total of ',i7)
            dl2 = real(il2, kind=DP)
            do il3=ilmin, ilmax
                dl3 = real(il3, kind=DP)

                ! Define boundaries for l1 given l2 and l3
                il1min = MAX( abs(il2-il3), abs(ispinl2_x+ispinl3_x) )
                il1max = il2 + il3
                dl1min = real(il1min, kind=DP)
                dl1max = real(il1max, kind=DP)

                ! Number of allowed l1 given l2 and l3
                indim = il1max - il1min + 1

                ! The fortran routine WIG3J compute all allowed l1 given fixed l2,l3,m2,m3.
                ! Therefore we loop over l2,l3 (m2,m3 are fixed and given by the spins), and at each iteration we grab a
                ! vector dW_l1l2l3 of dimension indim, which contains all 3-j symbols, that we will use to update the final Cl
                ! Be careful for the type of variables to pass. See wig3j_f.f for more infos.
                call WIG3J(dl2, dl3, dspinl2_x, dspinl3_x, dl1min, dl1max, &
                & dW_l1l2l3_x(il1min:il1max), indim, dier)
                if (flavor == 'clte') then
                    call WIG3J(dl2, dl3, dspinl2_y, dspinl3_y, dl1min, dl1max, &
                    & dW_l1l2l3_y(il1min:il1max), indim, dier)
                endif

                ! Define boundaries for l1 given l2 and l3 AND lmax
                ! il1max = MIN( il2 + il3 , ilmax )

                do il1=il1min, il1max
                    if (il1 .lt. 2) CYCLE
                    dl1 = real(il1, kind=DP)

                    ! For the scalar case, the wigner are non-zero only if iL is even
                    iL = il1 + il2 + il3
                    if (MOD(iL,2) .ne. 0 .and. flavor == 'cltt') CYCLE
                    ! For the EE, the weights are non-zero only if iL is even
                    if (MOD(iL,2) .ne. 0 .and. flavor == 'clee') CYCLE
                    ! For the TE, the weights are non-zero only if iL is even
                    if (MOD(iL,2) .ne. 0 .and. flavor == 'clte') CYCLE
                    ! For the BB, the weights are non-zero only if iL is even
                    if (MOD(iL,2) .ne. 1 .and. flavor == 'clbb') CYCLE

                    ! Compute prefactors in Eq 52 or 73 (astro-ph/0001303)
                    dL_l2l3l1 = dl1*(dl1+1.d0) - dl2*(dl2+1.d0) + dl3*(dl3+1.d0)
                    dS_l1l2l3 = SQRT( ( (2*dl1+1.d0)*(2*dl2+1.d0)*(2*dl3+1.d0) )/(16*pi) )

                    ! Compute function F in Eq 52 or 73 (astro-ph/0001303)
                    if (flavor == 'clte') then
                        dF_l2l3l1 = (dL_l2l3l1 * dS_l1l2l3)**2 * dW_l1l2l3_x(il1)*dW_l1l2l3_y(il1)
                    else
                        dF_l2l3l1 = (dL_l2l3l1 * dS_l1l2l3 * dW_l1l2l3_x(il1))**2
                    endif


                    ! Check if weird things happened (numerical instabilities or infinities)
                    if (dF_l2l3l1 .ne. dF_l2l3l1 .or. ABS(dF_l2l3l1) .gt. dinfinity) then
                        dF_l2l3l1 = 0.0_dp
                    endif

                    ! Update the Cl
                    dDcl2TTDcl1pp(il2,il3) = dDcl2TTDcl1pp(il2,il3) + &
                    & dclCMB(il1) * dF_l2l3l1 / (2*dl2 + 1.d0)

                enddo
                if (flavor == 'cltt') then
                    dZl2l3 = dl3 * (dl3 + 1.d0) * (2*dl3 + 1.d0) * dl2 * (dl2 + 1.d0) * dclCMB(il2) / (8 * pi)
                else if (flavor == 'clee') then
                    dZl2l3 = dl3 * (dl3 + 1.d0) * (2*dl3 + 1.d0) * (dl2**2 + dl2 - 4.d0) * dclCMB(il2) / (8 * pi)
                else if (flavor == 'clte') then
                    dZl2l3 = dl3 * (dl3 + 1.d0) * (2*dl3 + 1.d0) * (dl2**2 + dl2 - 2.d0) * dclCMB(il2) / (8 * pi)
                else if (flavor == 'clbb') then
                    dZl2l3 = 0.d0
                endif
                dDcl2TTDcl1pp(il2,il3) = dDcl2TTDcl1pp(il2,il3) - dZl2l3
            enddo
        enddo

        DEALLOCATE(dW_l1l2l3_x)
        DEALLOCATE(dW_l1l2l3_y)

    end subroutine

    subroutine compute_bias_zero_xyxy_mpi(dclxx,dclyy,dclxy,N0bias,il2range,flavor,ilmin, &
        & ispinl2_x,ispinl3_x,ispinl2_y,ispinl3_y, &
        & il2dim,ilmax)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! *Description:
        !   This function returns sum(g_ell*f_ell) = ( 2*ell + 1 )/N_0 AND NOT N_0 directly.
        !   This is done this way because the l2 values are splitted among processors,
        !   and it allows easy external parallelization (and therefore sum up terms with reduce with python)
        !   (Recall that N_0 = (2ell +1)/(sum gell fell))
        !   The routine WIG3J compute all allowed l1 given fixed l2,l3,m2,m3.
        !   Symmetries are already taken into account (even l1+l2+l3).
        !   Variables starting with /i/ are int, while variables starting with /d/ are double
        ! *Arguments:
        !   * dclCMB:   IN       unlensed CMB power spectrum (double, 1D array)
        !   * dclCMBexp:IN       Lensed CMB power spectrum with noise (double, 1D array)
        !   * il2range: IN       Values of l2 (1D array)
        !   * flavor:   IN       Name of the spectrum (character, cltt, clee, etc.)
        !   * ilmin:    IN       Minimum ell (int)
        !   * ispin2:   IN       spin value associated to l2 (i.e. m2) (int)
        !   * ispin3:   IN       spin value associated to l3 (i.e. m3) (int)
        !   * il2dim:   IN       Length of the il2range (int)
        !   * ilmax:    IN       Length of dclCMB and dclpp (int)
        !   * N0bias:   OUT      ( 2*ell + 1 )/N_0 (double, 1D array)
        ! *Author: Julien Peloton, University of Sussex
        ! *Revision history (mmyy):
        !   * 1015: creation
        !   * 1015: add polarization features
        !   * 1115: switch from function to subroutine
        ! *todo:
        !   * Add explicitely MPI inside the routine (instead of playing with il2range)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        implicit none

        integer, parameter       :: I4B = 4
        integer, parameter       :: DP = 8
        real(DP), parameter      :: pi = 3.141592

        integer(I4B), intent(in) :: ilmin,ilmax,il2dim
        integer(I4B), intent(in) :: ispinl2_x,ispinl3_x,ispinl2_y,ispinl3_y
        integer(I4B), intent(in) :: il2range(0:il2dim)
        real(DP), intent(in)     :: dclxx(0:1,0:2*ilmax), dclyy(0:1,0:2*ilmax), dclxy(0:1,0:2*ilmax)
        CHARACTER(LEN=13), intent(in) :: flavor
        real(DP), intent(out)    :: N0bias(0:ilmax)

        integer(I4B)             :: il1, il1min, il1max, il2, il3, iL, indim,il_tmp, parity
        real(DP)                 :: dl1, dl1min, dl1max, dl2, dl3
        real(DP)                 :: dL_l3l2l1, dL_l1l2l3, dS_l1l2l3, dF_l1l2l3
        real(DP)                 :: dinfinity, dier
        real(DP)                 :: dspinl2_x, dspinl3_x, dspinl2_y, dspinl3_y
        CHARACTER(LEN=13)        :: creturn

        real(DP), allocatable    :: dW_l1l2l3_x(:), dW_l1l2l3_y(:)

        ! Initialize
        dinfinity = 1.e20_dp
        ALLOCATE(dW_l1l2l3_x(0:2*ilmax))
        ALLOCATE(dW_l1l2l3_y(0:2*ilmax))
        dW_l1l2l3_x = 0.0_dp
        dW_l1l2l3_y = 0.0_dp
        dier = 0.0_dp

        dspinl2_x = real(ispinl2_x, kind=DP)
        dspinl3_x = real(ispinl3_x, kind=DP)
        dspinl2_y = real(ispinl2_y, kind=DP)
        dspinl3_y = real(ispinl3_y, kind=DP)

        if (flavor == 'cltttt' .or. flavor == 'cleeee' .or. &
            & flavor == 'clbbbb' .or. flavor == 'clttee' .or. &
            & flavor == 'clttte' .or. flavor == 'cleete') then
            parity = 0
        else if (flavor == 'cltete') then
            parity = 0
        else if (flavor == 'cltbtb') then
            parity = 1
        else if (flavor == 'clebeb') then
            parity = 1
        else if (flavor == 'clebtb') then
            parity = 1
        endif

        do il_tmp=0, il2dim
            il2 = il2range(il_tmp)
            creturn = achar(13)
            WRITE( * , 101 , ADVANCE='NO' ) creturn , il2 , ilmax
            101     FORMAT( a , 'ell : ',i7,' out of a total of ',i7)

            dl2 = real(il2, kind=DP)
            do il3=ilmin, ilmax
                dl3 = real(il3, kind=DP)

                ! Define boundaries for l1 given l2 and l3
                il1min = MAX( abs(il2-il3), abs(ispinl2_x+ispinl3_x) )
                il1max = il2 + il3
                dl1min = real(il1min, kind=DP)
                dl1max = real(il1max, kind=DP)

                ! Number of allowed l1 given l2 and l3
                indim = il1max - il1min + 1

                ! The fortran routine WIG3J compute all allowed l1 given fixed l2,l3,m2,m3.
                ! Therefore we loop over l2,l3 (m2,m3 are fixed and given by the spins), and at each iteration we grab a
                ! vector dW_l1l2l3 of dimension indim, which contains all 3-j symbols, that we will use to update the final Cl
                ! Be careful for the type of variables to pass. See wig3j_f.f for more infos.
                if (flavor == 'cltete' .or. flavor == 'clttee' .or. flavor == 'clttte' .or. flavor == 'cleete') then
                    call WIG3J(dl2, dl3, dspinl2_x, dspinl3_x, dl1min, dl1max, &
                    & dW_l1l2l3_x(il1min:il1max), indim, dier)
                    call WIG3J(dl2, dl3, dspinl2_y, dspinl3_y, dl1min, dl1max, &
                    & dW_l1l2l3_y(il1min:il1max), indim, dier)
                else
                    ! Includes TT, EE, BB, EB.
                    ! This includes TB as well which get only spin2 wigner.
                    ! use spin_y to include TB as well (spin_x = spin_y for all but TE)
                    call WIG3J(dl2, dl3, dspinl2_y, dspinl3_y, dl1min, dl1max, &
                    & dW_l1l2l3_x(il1min:il1max), indim, dier)
                endif

                il1max = MIN( il2 + il3 , ilmax)

                do il1=il1min, il1max
                    if (il1 .lt. 2) CYCLE
                    dl1 = real(il1, kind=DP)

                    ! For TT, EE, BB and TE cases, the wigner are non-zero only if iL is even
                    iL = il1 + il2 + il3
                    if (MOD(iL,2) .ne. 0 .and. parity .eq. 0) CYCLE
                    ! For the TB and EB cases, the weights are non-zero only if iL is odd
                    if (MOD(iL,2) .ne. 1 .and. parity .eq. 1) CYCLE

                    if (flavor == 'cltttt' .or. flavor == 'cleeee' .or. flavor == 'clbbbb') then
                        ! Compute prefactors in Eq 52 or 73 (astro-ph/0001303)
                        dL_l3l2l1 =  dl1*(dl1+1.d0) + dl2*(dl2+1.d0) - dl3*(dl3+1.d0)
                        dL_l1l2l3 = -dl1*(dl1+1.d0) + dl2*(dl2+1.d0) + dl3*(dl3+1.d0)
                        dS_l1l2l3 = SQRT( ( (2*dl1+1.d0)*(2*dl2+1.d0)*(2*dl3+1.d0) )/(16*pi) )

                        ! Compute function F in Eq 52 or 73 (astro-ph/0001303)
                        dF_l1l2l3 = (dclxx(0,il1) * dL_l3l2l1 + dclyy(0,il3) * dL_l1l2l3) * &
                        & dS_l1l2l3 * dW_l1l2l3_x(il1)

                        ! Check if weird things happened (numerical instabilities or infinities)
                        if (dF_l1l2l3 .ne. dF_l1l2l3 .or. ABS(dF_l1l2l3) .gt. dinfinity) then
                            dF_l1l2l3 = 0.0_dp
                        endif

                        ! Update the Cl
                        N0bias(il2) = N0bias(il2) + dF_l1l2l3**2 / (2.d0*dclxx(1,il3)*dclyy(1,il1))

                    else if (flavor == 'clttee') then
                        ! Compute prefactors in Eq 52 or 73 (astro-ph/0001303)
                        dL_l3l2l1 =  dl1*(dl1+1.d0) + dl2*(dl2+1.d0) - dl3*(dl3+1.d0)
                        dL_l1l2l3 = -dl1*(dl1+1.d0) + dl2*(dl2+1.d0) + dl3*(dl3+1.d0)
                        dS_l1l2l3 = SQRT( ( (2*dl1+1.d0)*(2*dl2+1.d0)*(2*dl3+1.d0) )/(16*pi) )

                        ! Compute function F in Eq 52 or 73 (astro-ph/0001303)
                        dF_l1l2l3 = (dclxx(0,il1) * dL_l3l2l1 + dclxx(0,il3) * dL_l1l2l3) * &
                        & dS_l1l2l3 * dW_l1l2l3_x(il1) * &
                        (dclyy(0,il1) * dL_l3l2l1 + dclyy(0,il3) * dL_l1l2l3) * &
                        & dS_l1l2l3 * dW_l1l2l3_y(il1)

                        ! Check if weird things happened (numerical instabilities or infinities)
                        if (dF_l1l2l3 .ne. dF_l1l2l3 .or. ABS(dF_l1l2l3) .gt. dinfinity) then
                            dF_l1l2l3 = 0.0_dp
                        endif

                        ! Update the Cl
                        N0bias(il2) = N0bias(il2) + dF_l1l2l3 * dclxy(1,il3)*dclxy(1,il1) / &
                        (2.d0*dclxx(1,il3)*dclxx(1,il1)*dclyy(1,il3)*dclyy(1,il1))

                    else if (flavor == 'cltete') then
                        ! Compute prefactors in Eq 52 or 73 (astro-ph/0001303)
                        dL_l3l2l1 =  dl1*(dl1+1.d0) + dl2*(dl2+1.d0) - dl3*(dl3+1.d0)
                        dL_l1l2l3 = -dl1*(dl1+1.d0) + dl2*(dl2+1.d0) + dl3*(dl3+1.d0)
                        dS_l1l2l3 = SQRT( ( (2*dl1+1.d0)*(2*dl2+1.d0)*(2*dl3+1.d0) )/(16*pi) )

                        ! Compute function F in Eq 52 or 73 (astro-ph/0001303)
                        dF_l1l2l3 = (dclxy(0,il1) * dL_l3l2l1 * dW_l1l2l3_x(il1) + &
                        &  dclxy(0,il3) * dL_l1l2l3 * dW_l1l2l3_y(il1)) * dS_l1l2l3

                        ! Check if weird things happened (numerical instabilities or infinities)
                        if (dF_l1l2l3 .ne. dF_l1l2l3 .or. ABS(dF_l1l2l3) .gt. dinfinity) then
                            dF_l1l2l3 = 0.0_dp
                        endif

                        ! Update the Cl
                        N0bias(il2) = N0bias(il2) + dF_l1l2l3**2 * &
                        & (dclxx(1,il1) * dclyy(1,il3) - dclxy(1,il3) * dclxy(1,il1)) / &
                        & (dclxx(1,il3) * dclxx(1,il1) * dclyy(1,il3) * dclyy(1,il1) - (dclxy(1,il3) * dclxy(1,il1))**2)

                    else if (flavor == 'clttte') then
                        ! Compute prefactors in Eq 52 or 73 (astro-ph/0001303)
                        dL_l3l2l1 =  dl1*(dl1+1.d0) + dl2*(dl2+1.d0) - dl3*(dl3+1.d0)
                        dL_l1l2l3 = -dl1*(dl1+1.d0) + dl2*(dl2+1.d0) + dl3*(dl3+1.d0)
                        dS_l1l2l3 = SQRT( ( (2*dl1+1.d0)*(2*dl2+1.d0)*(2*dl3+1.d0) )/(16*pi) )

                        ! Compute function F in Eq 52 or 73 (astro-ph/0001303)
                        dF_l1l2l3 = (dclxx(0,il1) * dL_l3l2l1 + dclxx(0,il3) * dL_l1l2l3) * dS_l1l2l3 * dW_l1l2l3_x(il1) * &
                        & (dclxy(0,il1) * dL_l3l2l1 * dW_l1l2l3_x(il1) + dclxy(0,il3) * dL_l1l2l3 * dW_l1l2l3_y(il1)) * dS_l1l2l3

                        ! Check if weird things happened (numerical instabilities or infinities)
                        if (dF_l1l2l3 .ne. dF_l1l2l3 .or. ABS(dF_l1l2l3) .gt. dinfinity) then
                            dF_l1l2l3 = 0.0_dp
                        endif

                        ! Update the Cl
                        N0bias(il2) = N0bias(il2) + dF_l1l2l3 * &
                        & ( (dclxx(1,il1)*dclyy(1,il3) - dclxy(1,il3)*dclxy(1,il1)) * dclxx(1,il3) * dclxy(1,il1) + &
                        &   (dclyy(1,il1)*dclxx(1,il3) - dclxy(1,il3)*dclxy(1,il1)) * dclxy(1,il3) * dclxx(1,il1) ) / &
                        & (2.d0 * dclxx(1,il3) * dclxx(1,il1) * &
                        & ( dclxx(1,il3)*dclxx(1,il1)*dclyy(1,il3)*dclyy(1,il1) - (dclxy(1,il3)*dclxy(1,il1))**2) )

                    else if (flavor == 'cleete') then
                        ! Compute prefactors in Eq 52 or 73 (astro-ph/0001303)
                        dL_l3l2l1 =  dl1*(dl1+1.d0) + dl2*(dl2+1.d0) - dl3*(dl3+1.d0)
                        dL_l1l2l3 = -dl1*(dl1+1.d0) + dl2*(dl2+1.d0) + dl3*(dl3+1.d0)
                        dS_l1l2l3 = SQRT( ( (2*dl1+1.d0)*(2*dl2+1.d0)*(2*dl3+1.d0) )/(16*pi) )

                        ! Compute function F in Eq 52 or 73 (astro-ph/0001303)
                        dF_l1l2l3 = (dclyy(0,il1) * dL_l3l2l1 + dclyy(0,il3) * dL_l1l2l3) * dS_l1l2l3 * dW_l1l2l3_y(il1) * &
                        & (dclxy(0,il1) * dL_l3l2l1 * dW_l1l2l3_x(il1) + dclxy(0,il3) * dL_l1l2l3 * dW_l1l2l3_y(il1)) * dS_l1l2l3

                        ! Check if weird things happened (numerical instabilities or infinities)
                        if (dF_l1l2l3 .ne. dF_l1l2l3 .or. ABS(dF_l1l2l3) .gt. dinfinity) then
                            dF_l1l2l3 = 0.0_dp
                        endif

                        ! Update the Cl
                        N0bias(il2) = N0bias(il2) + dF_l1l2l3 * &
                        & ( (dclxx(1,il1)*dclyy(1,il3) - dclxy(1,il3)*dclxy(1,il1)) * dclxy(1,il3) * dclyy(1,il1) + &
                        &   (dclyy(1,il1)*dclxx(1,il3) - dclxy(1,il3)*dclxy(1,il1)) * dclyy(1,il3) * dclxy(1,il1) ) / &
                        & (2.d0 * dclyy(1,il3) * dclyy(1,il1) * &
                        & ( dclxx(1,il3)*dclxx(1,il1)*dclyy(1,il3)*dclyy(1,il1) - (dclxy(1,il3)*dclxy(1,il1))**2) )

                    else if (flavor == 'clebtb') then
                        ! Compute prefactors in Eq 52 or 73 (astro-ph/0001303)
                        dL_l3l2l1 =  dl1*(dl1+1.d0) + dl2*(dl2+1.d0) - dl3*(dl3+1.d0)
                        dL_l1l2l3 = -dl1*(dl1+1.d0) + dl2*(dl2+1.d0) + dl3*(dl3+1.d0)
                        dS_l1l2l3 = SQRT( ( (2*dl1+1.d0)*(2*dl2+1.d0)*(2*dl3+1.d0) )/(16*pi) )

                        ! Compute function F in Eq 52 or 73 (astro-ph/0001303)
                        dF_l1l2l3 = (-dclyy(0,il1) * dL_l3l2l1 + dclxx(0,il3) * dL_l1l2l3) * dS_l1l2l3 * dW_l1l2l3_x(il1) * &
                        & dL_l1l2l3 * dS_l1l2l3 * dW_l1l2l3_x(il1)

                        ! Check if weird things happened (numerical instabilities or infinities)
                        if (dF_l1l2l3 .ne. dF_l1l2l3 .or. ABS(dF_l1l2l3) .gt. dinfinity) then
                            dF_l1l2l3 = 0.0_dp
                        endif

                        ! Update the Cl
                        N0bias(il2) = N0bias(il2) + dF_l1l2l3 * dclxy(1,il3) / dclyy(1,il1)

                    else if (flavor == 'clebeb') then
                        ! Compute prefactors in Eq 52 or 73 (astro-ph/0001303)
                        dL_l3l2l1 =  dl1*(dl1+1.d0) + dl2*(dl2+1.d0) - dl3*(dl3+1.d0)
                        dL_l1l2l3 = -dl1*(dl1+1.d0) + dl2*(dl2+1.d0) + dl3*(dl3+1.d0)
                        dS_l1l2l3 = SQRT( ( (2*dl1+1.d0)*(2*dl2+1.d0)*(2*dl3+1.d0) )/(16*pi) )

                        ! Compute function F in Eq 52 or 73 (astro-ph/0001303)
                        dF_l1l2l3 = (dclyy(0,il1) * dL_l3l2l1 - dclxx(0,il3) * dL_l1l2l3) * dS_l1l2l3 * dW_l1l2l3_x(il1)


                        ! Check if weird things happened (numerical instabilities or infinities)
                        if (dF_l1l2l3 .ne. dF_l1l2l3 .or. ABS(dF_l1l2l3) .gt. dinfinity) then
                            dF_l1l2l3 = 0.0_dp
                        endif

                        ! Update the Cl
                        N0bias(il2) = N0bias(il2) + dF_l1l2l3**2 / (dclxx(1,il3)*dclyy(1,il1))

                    else if (flavor == 'cltbtb') then
                        ! Compute prefactors in Eq 52 or 73 (astro-ph/0001303)
                        dL_l1l2l3 = -dl1*(dl1+1.d0) + dl2*(dl2+1.d0) + dl3*(dl3+1.d0)
                        dS_l1l2l3 = SQRT( ( (2*dl1+1.d0)*(2*dl2+1.d0)*(2*dl3+1.d0) )/(16*pi) )

                        ! For TB, cl_xy is TE (see Tab1 in Okamoto and Hu)
                        dF_l1l2l3 = dclxy(0,il3) * dL_l1l2l3 * dS_l1l2l3 * dW_l1l2l3_x(il1)

                        ! Check if weird things happened (numerical instabilities or infinities)
                        if (dF_l1l2l3 .ne. dF_l1l2l3 .or. ABS(dF_l1l2l3) .gt. dinfinity) then
                            dF_l1l2l3 = 0.0_dp
                        endif

                        ! Update the Cl
                        N0bias(il2) = N0bias(il2) + dF_l1l2l3**2 / (dclxx(1,il3)*dclyy(1,il1))
                    endif

                enddo
            enddo
        enddo

        DEALLOCATE(dW_l1l2l3_x)
        DEALLOCATE(dW_l1l2l3_y)

    end subroutine

    subroutine compute_lensed_spectra_xy_mpi(dclCMB,dclpp,clanalytic,il2range,flavor,ilmin,ispinl2_x,ispinl3_x, &
        & ispinl2_y,ispinl3_y,il2dim,ilmax)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! *Description:
        !   Given an unlensed spectrum (in the form XY, where X= T, E), this function computes
        !   the second term of RHS of Eq. 25 from 1308.0286 (perturbative computation for lensed spectrum).
        !   The routine WIG3J compute all allowed l1 given fixed l2,l3,m2,m3.
        !   Symmetries are already taken into account (even l1+l2+l3).
        !   Variables starting with /i/ are int, while variables starting with /d/ are double
        ! *Arguments:
        !   * dclCMB:   IN       unlensed CMB power spectrum (double, 1D array)
        !   * dclpp:    IN       lensing potential power spectrum (double, 1D array)
        !   * il2range: IN       Values of l2 (1D array)
        !   * flavor:   IN       Name of the spectrum (character, cltt, clee, etc.)
        !   * ilmin:    IN       Minimum ell (int)
        !   * ispin2_x:   IN       spin value associated to l2 (i.e. m2) (int) (spin0)
        !   * ispin3_x:   IN       spin value associated to l3 (i.e. m3) (int) (spin0)
        !   * ispin2_y:   IN       spin value associated to l2 (i.e. m2) (int) (spin2)
        !   * ispin3_y:   IN       spin value associated to l3 (i.e. m3) (int) (spin2)
        !   * il2dim:   IN       Length of the il2range (int)
        !   * ilmax:    IN       Length of dclCMB and dclpp (int)
        !   * clanalytic: OUT    Part of the lensed spectrum (double, 1D array)
        ! *Author: Julien Peloton, University of Sussex
        ! *Revision history (mmyy):
        !   * 1015: creation
        !   * 1015: add polarization features
        !   * 1115: switch from function to subroutine
        ! *todo:
        !   * Add explicitely MPI inside the routine (instead of playing with il2range)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        implicit none

        integer, parameter       :: I4B = 4
        integer, parameter       :: DP = 8
        real(DP), parameter      :: pi = 3.141592

        integer(I4B), intent(in) :: ilmin,ilmax,il2dim
        integer(I4B), intent(in) :: ispinl2_x,ispinl3_x,ispinl2_y,ispinl3_y
        integer(I4B), intent(in) :: il2range(0:il2dim)
        real(DP), intent(in)     :: dclCMB(0:1,0:2*ilmax), dclpp(0:2*ilmax)
        CHARACTER(LEN=13), intent(in) :: flavor
        real(DP), intent(out)    :: clanalytic(0:ilmax)

        integer(I4B)             :: il1, il1min, il1max, il2, il3, iL, indim, il_tmp
        real(DP)                 :: dl1, dl1min, dl1max, dl2, dl3
        real(DP)                 :: dL_l2l3l1, dS_l1l2l3, dF_l2l3l1, dinfinity, dier
        real(DP)                 :: dspinl2_x, dspinl3_x, dspinl2_y, dspinl3_y

        real(DP), allocatable    :: dW_l1l2l3_0(:), dW_l1l2l3_2(:)
        CHARACTER(LEN=13)        :: creturn


        !Initialize cl array
        do il1=0, ilmax
            clanalytic(il1)=0.0_dp
        enddo

        ! Initialize
        dinfinity = 1.e20_dp
        ALLOCATE(dW_l1l2l3_0(0:2*ilmax))
        ALLOCATE(dW_l1l2l3_2(0:2*ilmax))
        dW_l1l2l3_0 = 0.0_dp
        dW_l1l2l3_2 = 0.0_dp
        dier = 0.0_dp

        dspinl2_x = real(ispinl2_x, kind=DP)
        dspinl3_x = real(ispinl3_x, kind=DP)
        dspinl2_y = real(ispinl2_y, kind=DP)
        dspinl3_y = real(ispinl3_y, kind=DP)

        do il_tmp=0, il2dim
            il2 = il2range(il_tmp)
            creturn = achar(13)
            WRITE( * , 101 , ADVANCE='NO' ) creturn , il2 , ilmax
            101     FORMAT( a , 'ell : ',i7,' out of a total of ',i7)
            dl2 = real(il2, kind=DP)
            do il3=ilmin, ilmax
                dl3 = real(il3, kind=DP)

                ! Define boundaries for l1 given l2 and l3
                il1min = MAX( abs(il2-il3), abs(ispinl2_x+ispinl3_x) )
                il1max = il2 + il3
                dl1min = real(il1min, kind=DP)
                dl1max = real(il1max, kind=DP)

                ! Number of allowed l1 given l2 and l3
                indim = il1max - il1min + 1

                ! The fortran routine WIG3J compute all allowed l1 given fixed l2,l3,m2,m3.
                ! Therefore we loop over l2,l3 (m2,m3 are fixed and given by the spins), and at each iteration we grab a
                ! vector dW_l1l2l3 of dimension indim, which contains all 3-j symbols, that we will use to update the final Cl
                ! Be careful for the type of variables to pass. See wig3j_f.f for more infos.
                call WIG3J(dl2, dl3, dspinl2_x, dspinl3_x, dl1min, dl1max, &
                & dW_l1l2l3_0(il1min:il1max), indim, dier)
                call WIG3J(dl2, dl3, dspinl2_y, dspinl3_y, dl1min, dl1max, &
                & dW_l1l2l3_2(il1min:il1max), indim, dier)

                ! Define boundaries for l1 given l2 and l3 AND lmax
                ! il1max = MIN( il2 + il3 , ilmax )

                do il1=il1min, il1max
                    if (il1 .lt. 2) CYCLE
                    dl1 = real(il1, kind=DP)

                    ! For the scalar case, the wigner are non-zero only if iL is even
                    iL = il1 + il2 + il3
                    if (MOD(iL,2) .ne. 0 .and. flavor == 'cltt') CYCLE
                    ! For the EE, the weights are non-zero only if iL is even
                    if (MOD(iL,2) .ne. 0 .and. flavor == 'clee') CYCLE
                    ! For the TE, the weights are non-zero only if iL is even
                    if (MOD(iL,2) .ne. 0 .and. flavor == 'clte') CYCLE
                    ! For the EE, the weights are non-zero only if iL is even
                    if (MOD(iL,2) .ne. 1 .and. flavor == 'clbb') CYCLE

                    ! Compute prefactors in Eq 52 or 73 (astro-ph/0001303)
                    dL_l2l3l1 = dl1*(dl1+1.d0) - dl2*(dl2+1.d0) + dl3*(dl3+1.d0)
                    dS_l1l2l3 = SQRT( ( (2*dl1+1.d0)*(2*dl2+1.d0)*(2*dl3+1.d0) )/(16*pi) )

                    ! Compute function F in Eq 52 or 73 (astro-ph/0001303)
                    dF_l2l3l1 = dL_l2l3l1**2 * dS_l1l2l3**2 * dW_l1l2l3_0(il1) * dW_l1l2l3_2(il1)

                    ! Check if weird things happened (numerical instabilities or infinities)
                    if (dF_l2l3l1 .ne. dF_l2l3l1 .or. ABS(dF_l2l3l1) .gt. dinfinity) then
                        dF_l2l3l1 = 0.0_dp
                    endif

                    ! Update the Cl
                    clanalytic(il2) = clanalytic(il2) + ( dclCMB(0,il1) + (-1)**iL * dclCMB(1,il1) ) *&
                    & dclpp(il3) * dF_l2l3l1 / (2*dl2 + 1.d0)

                enddo
            enddo
        enddo

        DEALLOCATE(dW_l1l2l3_0)
        DEALLOCATE(dW_l1l2l3_2)

    end subroutine

    subroutine XX_phiphi_disconnected_mpi(dclxx,dclyy,dclxy,il2range, &
        & Bl2l3, &
        & flavor, &
        & ilmin,ispinl2_x,ispinl3_x,ispinl2_y,ispinl3_y, &
        & il2dim,ilmax )
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! *Description:
        !   This subroutine computes part of the derivative dN^{(0)}/dC^{CMB}
        !   The subroutine WIG3J compute all allowed l1 given fixed l2,l3,m2,m3.
        !   Symmetries are already taken into account (even l1+l2+l3).
        !   Variables starting with /i/ are int, while variables starting with /d/ are double
        ! *Arguments:
        !   * dclCMB:   IN       unlensed CMB power spectrum (double, 1D array)
        !   * dclCMBexp:IN       Lensed CMB power spectrum with noise (double, 1D array)
        !   * il2range: IN       Values of l2 (1D array)
        !   * flavor:   IN       Name of the spectrum (character, cltt, clee, etc.)
        !   * ilmin:    IN       Minimum ell (int)
        !   * ispin2:   IN       spin value associated to l2 (i.e. m2) (int)
        !   * ispin3:   IN       spin value associated to l3 (i.e. m3) (int)
        !   * il2dim:   IN       Length of the il2range (int)
        !   * ilmax:    IN       Length of dclCMB and dclpp (int)
        !   * Bl2l3:    OUT      Derivatives of N0 wrt C_{\ell,exp}^{XX} (double, 2D array)
        ! *Author: Julien Peloton, University of Sussex
        ! *Revision history (mmyy):
        !   * 1015: creation
        !   * 1015: add polarization features
        ! *todo:
        !   * Add explicitely MPI inside the routine (instead of playing with il2range)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        implicit none

        integer, parameter       :: I4B = 4
        integer, parameter       :: DP = 8
        real(DP), parameter      :: pi = 3.141592

        integer(I4B), intent(in) :: ilmin,ilmax,il2dim
        integer(I4B), intent(in) :: ispinl2_x,ispinl3_x,ispinl2_y,ispinl3_y
        integer(I4B), intent(in) :: il2range(0:il2dim)
        real(DP), intent(in)     :: dclxx(0:1,0:2*ilmax), dclyy(0:1,0:2*ilmax), dclxy(0:1,0:2*ilmax)
        CHARACTER(LEN=13), intent(in) :: flavor
        real(DP), intent(out)  :: Bl2l3(0:ilmax,0:ilmax)

        integer(I4B)             :: il1, il1min, il1max, il2, il3, iL, indim,il_tmp, parity
        real(DP)                 :: dl1, dl1min, dl1max, dl2, dl3
        real(DP)                 :: dL_l3l2l1, dL_l1l2l3, dS_l1l2l3, dF_l1l2l3
        real(DP)                 :: dinfinity, dier
        real(DP)                 :: dspinl2_x, dspinl3_x, dspinl2_y, dspinl3_y
        CHARACTER(LEN=13)        :: creturn

        real(DP), allocatable    :: dW_l1l2l3_x(:), dW_l1l2l3_y(:)

        ! Initialize
        dinfinity = 1.e20_dp
        ALLOCATE(dW_l1l2l3_x(0:2*ilmax))
        ALLOCATE(dW_l1l2l3_y(0:2*ilmax))
        dW_l1l2l3_x = 0.0_dp
        dW_l1l2l3_y = 0.0_dp
        dier = 0.0_dp

        dspinl2_x = real(ispinl2_x, kind=DP)
        dspinl3_x = real(ispinl3_x, kind=DP)
        dspinl2_y = real(ispinl2_y, kind=DP)
        dspinl3_y = real(ispinl3_y, kind=DP)

        if (flavor == 'cltt' .or. flavor == 'clee' .or. flavor == 'clbb') then
            parity = 0
        else if (flavor == 'clte') then
            parity = 0
        else if (flavor == 'cltete_te' .or. flavor == 'cltete_tt' .or. flavor == 'cltete_ee') then
            parity = 0
        else if (flavor == 'clttee_te') then
            parity = 0
        else if (flavor == 'clttte_tt' .or. flavor == 'clttte_te' .or. flavor == 'cleete_ee' .or. flavor == 'cleete_te') then
            parity = 0
        else if (flavor == 'cltb') then
            parity = 1
        else if (flavor == 'clebeb_ee' .or. flavor == 'clebeb_bb') then
            parity = 1
        endif

        do il_tmp=0, il2dim
            il2 = il2range(il_tmp)
            creturn = achar(13)
            WRITE( * , 101 , ADVANCE='NO' ) creturn , il2 , ilmax
            101     FORMAT( a , 'ell : ',i7,' out of a total of ',i7)

            dl2 = real(il2, kind=DP)
            do il3=ilmin, ilmax
                dl3 = real(il3, kind=DP)

                ! Define boundaries for l1 given l2 and l3
                ! il1min = MAX( abs(il2-il3), abs(ispinl2+ispinl3) )
                il1min = MAX( abs(il2-il3), abs(ispinl2_x+ispinl3_x) )
                il1max = il2 + il3
                dl1min = real(il1min, kind=DP)
                dl1max = real(il1max, kind=DP)

                ! Number of allowed l1 given l2 and l3
                indim = il1max - il1min + 1

                ! The fortran routine WIG3J compute all allowed l1 given fixed l2,l3,m2,m3.
                ! Therefore we loop over l2,l3 (m2,m3 are fixed and given by the spins), and at each iteration we grab a
                ! vector dW_l1l2l3 of dimension indim, which contains all 3-j symbols, that we will use to update the final Cl
                ! Be careful for the type of variables to pass. See wig3j_f.f for more infos.
                if (flavor == 'clte' .or. &
                & flavor == 'cltete_te' .or. flavor == 'cltete_tt' .or. flavor == 'cltete_ee' .or. &
                & flavor == 'clttee_te' .or. flavor == 'clttte_tt' .or. flavor == 'clttte_te' .or. &
                & flavor == 'cleete_ee' .or. flavor == 'cleete_te') then
                    call WIG3J(dl2, dl3, dspinl2_x, dspinl3_x, dl1min, dl1max, &
                    & dW_l1l2l3_x(il1min:il1max), indim, dier)
                    call WIG3J(dl2, dl3, dspinl2_y, dspinl3_y, dl1min, dl1max, &
                    & dW_l1l2l3_y(il1min:il1max), indim, dier)
                else
                    ! Includes TT, EE, BB, EB.
                    ! This includes TB as well which get only spin2 wigner.
                    ! use spin_y to include TB as well (spin_x = spin_y for all but TE)
                    call WIG3J(dl2, dl3, dspinl2_y, dspinl3_y, dl1min, dl1max, &
                    & dW_l1l2l3_x(il1min:il1max), indim, dier)
                endif

                ! Define boundaries for l1 given l2 and l3 AND lmax
                il1max = MIN( il2 + il3 , ilmax)

                do il1=il1min, il1max
                    ! need to remove first two multipoles from the sum
                    if (il1 .lt. 2) CYCLE
                    dl1 = real(il1, kind=DP)

                    ! For TT, EE, BB and TE cases, the wigner are non-zero only if iL is even
                    iL = il1 + il2 + il3
                    if (MOD(iL,2) .ne. 0 .and. parity .eq. 0) CYCLE
                    ! For the TB and EB cases, the weights are non-zero only if iL is odd
                    if (MOD(iL,2) .ne. 1 .and. parity .eq. 1) CYCLE

                    if (flavor == 'cltt' .or. flavor == 'clee' .or. flavor == 'clbb') then
                        ! Compute prefactors in Eq 52 or 73 (astro-ph/0001303)
                        dL_l3l2l1 =  dl1*(dl1+1.d0) + dl2*(dl2+1.d0) - dl3*(dl3+1.d0)
                        dL_l1l2l3 = -dl1*(dl1+1.d0) + dl2*(dl2+1.d0) + dl3*(dl3+1.d0)
                        dS_l1l2l3 = SQRT( ( (2*dl1+1.d0)*(2*dl2+1.d0)*(2*dl3+1.d0) )/(16*pi) )

                        ! Compute function F in Eq 52 or 73 (astro-ph/0001303)
                        dF_l1l2l3 = (dclxx(0,il1) * dL_l3l2l1 + dclyy(0,il3) * dL_l1l2l3) * &
                        & dS_l1l2l3 * dW_l1l2l3_x(il1) / (2*dclxx(1,il3)*dclyy(1,il1))

                        ! Check if weird things happened (numerical instabilities or infinities)
                        if (dF_l1l2l3 .ne. dF_l1l2l3 .or. ABS(dF_l1l2l3) .gt. dinfinity) then
                            dF_l1l2l3 = 0.0_dp
                        endif

                        ! Update the Cl
                        Bl2l3(il2,il3) = Bl2l3(il2,il3) + dF_l1l2l3**2 * dclxx(1,il1)

                    else if (flavor == 'clttee_te') then
                        ! Compute prefactors in Eq 52 or 73 (astro-ph/0001303)
                        dL_l3l2l1 =  dl1*(dl1+1.d0) + dl2*(dl2+1.d0) - dl3*(dl3+1.d0)
                        dL_l1l2l3 = -dl1*(dl1+1.d0) + dl2*(dl2+1.d0) + dl3*(dl3+1.d0)
                        dS_l1l2l3 = SQRT( ( (2*dl1+1.d0)*(2*dl2+1.d0)*(2*dl3+1.d0) )/(16*pi) )

                        ! Compute function F in Eq 52 or 73 (astro-ph/0001303)
                        dF_l1l2l3 = (dclxx(0,il1) * dL_l3l2l1 + dclxx(0,il3) * dL_l1l2l3) * &
                        & dS_l1l2l3 * dW_l1l2l3_x(il1) * &
                        (dclyy(0,il1) * dL_l3l2l1 + dclyy(0,il3) * dL_l1l2l3) * &
                        & dS_l1l2l3 * dW_l1l2l3_y(il1)

                        ! Check if weird things happened (numerical instabilities or infinities)
                        if (dF_l1l2l3 .ne. dF_l1l2l3 .or. ABS(dF_l1l2l3) .gt. dinfinity) then
                            dF_l1l2l3 = 0.0_dp
                        endif

                        ! Update the Cl
                        Bl2l3(il2,il3) = Bl2l3(il2,il3) + dF_l1l2l3 * dclxy(1,il1) / &
                        (4.d0*dclxx(1,il3)*dclxx(1,il1)*dclyy(1,il3)*dclyy(1,il1))

                    else if (flavor == 'clttte_tt') then
                        ! Compute prefactors in Eq 52 or 73 (astro-ph/0001303)
                        dL_l3l2l1 =  dl1*(dl1+1.d0) + dl2*(dl2+1.d0) - dl3*(dl3+1.d0)
                        dL_l1l2l3 = -dl1*(dl1+1.d0) + dl2*(dl2+1.d0) + dl3*(dl3+1.d0)
                        dS_l1l2l3 = SQRT( ( (2*dl1+1.d0)*(2*dl2+1.d0)*(2*dl3+1.d0) )/(16*pi) )

                        ! Compute function F in Eq 52 or 73 (astro-ph/0001303)
                        dF_l1l2l3 = (dclxx(0,il1) * dL_l3l2l1 + dclxx(0,il3) * dL_l1l2l3) * dS_l1l2l3 * dW_l1l2l3_x(il1) * &
                        & (dclxy(0,il1) * dL_l3l2l1 * dW_l1l2l3_x(il1) + dclxy(0,il3) * dL_l1l2l3 * dW_l1l2l3_y(il1)) * dS_l1l2l3

                        ! Check if weird things happened (numerical instabilities or infinities)
                        if (dF_l1l2l3 .ne. dF_l1l2l3 .or. ABS(dF_l1l2l3) .gt. dinfinity) then
                            dF_l1l2l3 = 0.0_dp
                        endif

                        ! Update the Cl
                        Bl2l3(il2,il3) = Bl2l3(il2,il3) + dF_l1l2l3 * &
                        & (dclxx(1,il1)*dclyy(1,il3) - dclxy(1,il3)*dclxy(1,il1)) * dclxy(1,il1) / &
                        & (2.d0 * dclxx(1,il3) * dclxx(1,il1) * &
                        & ( dclxx(1,il3)*dclxx(1,il1)*dclyy(1,il3)*dclyy(1,il1) - (dclxy(1,il3)*dclxy(1,il1))**2) )

                    else if (flavor == 'clttte_te') then
                        ! Compute prefactors in Eq 52 or 73 (astro-ph/0001303)
                        dL_l3l2l1 =  dl1*(dl1+1.d0) + dl2*(dl2+1.d0) - dl3*(dl3+1.d0)
                        dL_l1l2l3 = -dl1*(dl1+1.d0) + dl2*(dl2+1.d0) + dl3*(dl3+1.d0)
                        dS_l1l2l3 = SQRT( ( (2*dl1+1.d0)*(2*dl2+1.d0)*(2*dl3+1.d0) )/(16*pi) )

                        ! Compute function F in Eq 52 or 73 (astro-ph/0001303)
                        dF_l1l2l3 = (dclxx(0,il1) * dL_l3l2l1 + dclxx(0,il3) * dL_l1l2l3) * dS_l1l2l3 * dW_l1l2l3_x(il1) * &
                        & (dclxy(0,il1) * dL_l3l2l1 * dW_l1l2l3_x(il1) + dclxy(0,il3) * dL_l1l2l3 * dW_l1l2l3_y(il1)) * dS_l1l2l3

                        ! Check if weird things happened (numerical instabilities or infinities)
                        if (dF_l1l2l3 .ne. dF_l1l2l3 .or. ABS(dF_l1l2l3) .gt. dinfinity) then
                            dF_l1l2l3 = 0.0_dp
                        endif

                        ! Update the Cl
                        Bl2l3(il2,il3) = Bl2l3(il2,il3) + dF_l1l2l3 * &
                        & (dclyy(1,il1)*dclxx(1,il3) - dclxy(1,il3)*dclxy(1,il1)) * dclxx(1,il1) / &
                        & (2.d0 * dclxx(1,il3) * dclxx(1,il1) * &
                        & ( dclxx(1,il3)*dclxx(1,il1)*dclyy(1,il3)*dclyy(1,il1) - (dclxy(1,il3)*dclxy(1,il1))**2) )

                    else if (flavor == 'cleete_ee') then
                        ! Compute prefactors in Eq 52 or 73 (astro-ph/0001303)
                        dL_l3l2l1 =  dl1*(dl1+1.d0) + dl2*(dl2+1.d0) - dl3*(dl3+1.d0)
                        dL_l1l2l3 = -dl1*(dl1+1.d0) + dl2*(dl2+1.d0) + dl3*(dl3+1.d0)
                        dS_l1l2l3 = SQRT( ( (2*dl1+1.d0)*(2*dl2+1.d0)*(2*dl3+1.d0) )/(16*pi) )

                        ! Compute function F in Eq 52 or 73 (astro-ph/0001303)
                        dF_l1l2l3 = (dclyy(0,il1) * dL_l3l2l1 + dclyy(0,il3) * dL_l1l2l3) * dS_l1l2l3 * dW_l1l2l3_y(il1) * &
                        & (dclxy(0,il1) * dL_l3l2l1 * dW_l1l2l3_x(il1) + dclxy(0,il3) * dL_l1l2l3 * dW_l1l2l3_y(il1)) * dS_l1l2l3

                        ! Check if weird things happened (numerical instabilities or infinities)
                        if (dF_l1l2l3 .ne. dF_l1l2l3 .or. ABS(dF_l1l2l3) .gt. dinfinity) then
                            dF_l1l2l3 = 0.0_dp
                        endif

                        ! Update the Cl
                        Bl2l3(il2,il3) = Bl2l3(il2,il3) + dF_l1l2l3 * &
                        & (dclyy(1,il1)*dclxx(1,il3) - dclxy(1,il3)*dclxy(1,il1)) * dclxy(1,il1) / &
                        & (2.d0 * dclyy(1,il3) * dclyy(1,il1) * &
                        & ( dclxx(1,il3)*dclxx(1,il1)*dclyy(1,il3)*dclyy(1,il1) - (dclxy(1,il3)*dclxy(1,il1))**2) )

                    else if (flavor == 'cleete_te') then
                        ! Compute prefactors in Eq 52 or 73 (astro-ph/0001303)
                        dL_l3l2l1 =  dl1*(dl1+1.d0) + dl2*(dl2+1.d0) - dl3*(dl3+1.d0)
                        dL_l1l2l3 = -dl1*(dl1+1.d0) + dl2*(dl2+1.d0) + dl3*(dl3+1.d0)
                        dS_l1l2l3 = SQRT( ( (2*dl1+1.d0)*(2*dl2+1.d0)*(2*dl3+1.d0) )/(16*pi) )

                        ! Compute function F in Eq 52 or 73 (astro-ph/0001303)
                        dF_l1l2l3 = (dclyy(0,il1) * dL_l3l2l1 + dclyy(0,il3) * dL_l1l2l3) * dS_l1l2l3 * dW_l1l2l3_y(il1) * &
                        & (dclxy(0,il1) * dL_l3l2l1 * dW_l1l2l3_x(il1) + dclxy(0,il3) * dL_l1l2l3 * dW_l1l2l3_y(il1)) * dS_l1l2l3

                        ! Check if weird things happened (numerical instabilities or infinities)
                        if (dF_l1l2l3 .ne. dF_l1l2l3 .or. ABS(dF_l1l2l3) .gt. dinfinity) then
                            dF_l1l2l3 = 0.0_dp
                        endif

                        ! Update the Cl
                        Bl2l3(il2,il3) = Bl2l3(il2,il3) + dF_l1l2l3 * &
                        & (dclxx(1,il1)*dclyy(1,il3) - dclxy(1,il3)*dclxy(1,il1)) * dclyy(1,il1) / &
                        & (2.d0 * dclyy(1,il3) * dclyy(1,il1) * &
                        & ( dclxx(1,il3)*dclxx(1,il1)*dclyy(1,il3)*dclyy(1,il1) - (dclxy(1,il3)*dclxy(1,il1))**2) )

                    else if (flavor == 'cltete_te') then
                        ! Compute prefactors in Eq 52 or 73 (astro-ph/0001303)
                        dL_l3l2l1 =  dl1*(dl1+1.d0) + dl2*(dl2+1.d0) - dl3*(dl3+1.d0)
                        dL_l1l2l3 = -dl1*(dl1+1.d0) + dl2*(dl2+1.d0) + dl3*(dl3+1.d0)
                        dS_l1l2l3 = SQRT( ( (2*dl1+1.d0)*(2*dl2+1.d0)*(2*dl3+1.d0) )/(16*pi) )

                        ! Compute function F in Eq 52 or 73 (astro-ph/0001303)
                        dF_l1l2l3 = (dclxy(0,il1) * dL_l3l2l1 * dW_l1l2l3_x(il1) + &
                        &  dclxy(0,il3) * dL_l1l2l3 * dW_l1l2l3_y(il1)) * dS_l1l2l3

                        ! Check if weird things happened (numerical instabilities or infinities)
                        if (dF_l1l2l3 .ne. dF_l1l2l3 .or. ABS(dF_l1l2l3) .gt. dinfinity) then
                            dF_l1l2l3 = 0.0_dp
                        endif

                        ! Update the Cl
                        Bl2l3(il2,il3) = Bl2l3(il2,il3) + dF_l1l2l3**2 * &
                        & (dclxx(1,il1)*dclyy(1,il3) - dclxy(1,il3)*dclxy(1,il1)) * &
                        & (dclyy(1,il1)*dclxx(1,il3) - dclxy(1,il3)*dclxy(1,il1)) / &
                        & (dclxx(1,il3) * dclxx(1,il1) * dclyy(1,il3) * dclyy(1,il1) - &
                        & (dclxy(1,il3) * dclxy(1,il1))**2)**2 * &
                        & dclxy(1,il1)

                    else if (flavor == 'cltete_tt') then
                        ! Compute prefactors in Eq 52 or 73 (astro-ph/0001303)
                        dL_l3l2l1 =  dl1*(dl1+1.d0) + dl2*(dl2+1.d0) - dl3*(dl3+1.d0)
                        dL_l1l2l3 = -dl1*(dl1+1.d0) + dl2*(dl2+1.d0) + dl3*(dl3+1.d0)
                        dS_l1l2l3 = SQRT( ( (2*dl1+1.d0)*(2*dl2+1.d0)*(2*dl3+1.d0) )/(16*pi) )

                        ! Compute function F in Eq 52 or 73 (astro-ph/0001303)
                        dF_l1l2l3 = (dclxy(0,il1) * dL_l3l2l1 * dW_l1l2l3_x(il1) + &
                        &  dclxy(0,il3) * dL_l1l2l3 * dW_l1l2l3_y(il1)) * dS_l1l2l3

                        ! Check if weird things happened (numerical instabilities or infinities)
                        if (dF_l1l2l3 .ne. dF_l1l2l3 .or. ABS(dF_l1l2l3) .gt. dinfinity) then
                            dF_l1l2l3 = 0.0_dp
                        endif

                        ! Update the Cl
                        Bl2l3(il2,il3) = Bl2l3(il2,il3) + dF_l1l2l3**2 * &
                        & (dclxx(1,il1)*dclyy(1,il3) - dclxy(1,il3)*dclxy(1,il1)) * &
                        & (dclxx(1,il1)*dclyy(1,il3) - dclxy(1,il3)*dclxy(1,il1)) / &
                        & (dclxx(1,il3) * dclxx(1,il1) * dclyy(1,il3) * dclyy(1,il1) - &
                        & (dclxy(1,il3) * dclxy(1,il1))**2)**2 * &
                        & dclyy(1,il1)

                    else if (flavor == 'cltete_ee') then
                        ! Compute prefactors in Eq 52 or 73 (astro-ph/0001303)
                        dL_l3l2l1 =  dl1*(dl1+1.d0) + dl2*(dl2+1.d0) - dl3*(dl3+1.d0)
                        dL_l1l2l3 = -dl1*(dl1+1.d0) + dl2*(dl2+1.d0) + dl3*(dl3+1.d0)
                        dS_l1l2l3 = SQRT( ( (2*dl1+1.d0)*(2*dl2+1.d0)*(2*dl3+1.d0) )/(16*pi) )

                        ! Compute function F in Eq 52 or 73 (astro-ph/0001303)
                        dF_l1l2l3 = (dclxy(0,il1) * dL_l3l2l1 * dW_l1l2l3_x(il1) + &
                        &  dclxy(0,il3) * dL_l1l2l3 * dW_l1l2l3_y(il1)) * dS_l1l2l3

                        ! Check if weird things happened (numerical instabilities or infinities)
                        if (dF_l1l2l3 .ne. dF_l1l2l3 .or. ABS(dF_l1l2l3) .gt. dinfinity) then
                            dF_l1l2l3 = 0.0_dp
                        endif

                        ! Update the Cl
                        Bl2l3(il2,il3) = Bl2l3(il2,il3) + dF_l1l2l3**2 * &
                        & (dclyy(1,il1)*dclxx(1,il3) - dclxy(1,il3)*dclxy(1,il1)) * &
                        & (dclyy(1,il1)*dclxx(1,il3) - dclxy(1,il3)*dclxy(1,il1)) / &
                        & (dclxx(1,il3) * dclxx(1,il1) * dclyy(1,il3) * dclyy(1,il1) - &
                        & (dclxy(1,il3) * dclxy(1,il1))**2)**2 * &
                        & dclxx(1,il1)

                    else if (flavor == 'clebeb_ee') then
                        ! Compute prefactors in Eq 52 or 73 (astro-ph/0001303)
                        dL_l3l2l1 =  dl1*(dl1+1.d0) + dl2*(dl2+1.d0) - dl3*(dl3+1.d0)
                        dL_l1l2l3 = -dl1*(dl1+1.d0) + dl2*(dl2+1.d0) + dl3*(dl3+1.d0)
                        dS_l1l2l3 = SQRT( ( (2*dl1+1.d0)*(2*dl2+1.d0)*(2*dl3+1.d0) )/(16*pi) )

                        ! Compute function F in Eq 52 or 73 (astro-ph/0001303)
                        dF_l1l2l3 = (dclyy(0,il1) * dL_l3l2l1 - dclxx(0,il3) * dL_l1l2l3) * &
                        & dS_l1l2l3 * dW_l1l2l3_x(il1) / (dclxx(1,il3)*dclyy(1,il1))

                        ! Check if weird things happened (numerical instabilities or infinities)
                        if (dF_l1l2l3 .ne. dF_l1l2l3 .or. ABS(dF_l1l2l3) .gt. dinfinity) then
                            dF_l1l2l3 = 0.0_dp
                        endif

                        ! Update the Cl
                        Bl2l3(il2,il3) = Bl2l3(il2,il3) + dF_l1l2l3**2 * dclyy(1,il1)

                    else if (flavor == 'clebeb_bb') then
                        ! Compute prefactors in Eq 52 or 73 (astro-ph/0001303)
                        dL_l3l2l1 =  dl1*(dl1+1.d0) + dl2*(dl2+1.d0) - dl3*(dl3+1.d0)
                        dL_l1l2l3 = -dl1*(dl1+1.d0) + dl2*(dl2+1.d0) + dl3*(dl3+1.d0)
                        dS_l1l2l3 = SQRT( ( (2*dl1+1.d0)*(2*dl2+1.d0)*(2*dl3+1.d0) )/(16*pi) )

                        ! Compute function F in Eq 52 or 73 (astro-ph/0001303)
                        dF_l1l2l3 = (dclyy(0,il1) * dL_l3l2l1 - dclxx(0,il3) * dL_l1l2l3) * &
                        & dS_l1l2l3 * dW_l1l2l3_x(il1) / (dclxx(1,il3)*dclyy(1,il1))

                        ! Check if weird things happened (numerical instabilities or infinities)
                        if (dF_l1l2l3 .ne. dF_l1l2l3 .or. ABS(dF_l1l2l3) .gt. dinfinity) then
                            dF_l1l2l3 = 0.0_dp
                        endif

                        ! Update the Cl
                        Bl2l3(il2,il3) = Bl2l3(il2,il3) + dF_l1l2l3**2 * dclxx(1,il1)

                    else if (flavor == 'cltb') then
                        ! Not done yet
                        print *, 'Not available yet.'
                    endif

                enddo
            enddo
        enddo

        DEALLOCATE(dW_l1l2l3_x)
        DEALLOCATE(dW_l1l2l3_y)

    end subroutine

    subroutine precompute_trispb_fortran(dclxx1,dclyy1,dclxy1, &
        & dclxx2,dclyy2,dclxy2, &
        & il2range, &
        & Bl2l3, &
        & flavor1,flavor2,ilmin, &
        & ispinl2_x1,ispinl3_x1,ispinl2_y1,ispinl3_y1, &
        & ispinl2_x2,ispinl3_x2,ispinl2_y2,ispinl3_y2, &
        & il2dim,ilmax )
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! *Description:
        !   This subroutine computes parts used for the Type B trispectrum.
        !   Variables starting with /i/ are int, while variables starting with /d/ are double
        ! *Arguments:
        ! *Author: Julien Peloton, University of Sussex
        ! *Revision history (mmyy):
        !   * 1015: creation
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        implicit none

        integer, parameter       :: I4B = 4
        integer, parameter       :: DP = 8
        real(DP), parameter      :: pi = 3.141592

        integer(I4B), intent(in) :: ilmin,ilmax,il2dim
        integer(I4B), intent(in) :: ispinl2_x1,ispinl3_x1,ispinl2_y1,ispinl3_y1
        integer(I4B), intent(in) :: ispinl2_x2,ispinl3_x2,ispinl2_y2,ispinl3_y2
        integer(I4B), intent(in) :: il2range(0:il2dim)
        real(DP), intent(in)     :: dclxx1(0:1,0:2*ilmax), dclyy1(0:1,0:2*ilmax), dclxy1(0:1,0:2*ilmax)
        real(DP), intent(in)     :: dclxx2(0:1,0:2*ilmax), dclyy2(0:1,0:2*ilmax), dclxy2(0:1,0:2*ilmax)
        CHARACTER(LEN=13), intent(in) :: flavor1, flavor2
        real(DP), intent(out)  :: Bl2l3(0:ilmax,0:ilmax)

        integer(I4B)             :: il1, il1min, il1max, il2, il3, iL, indim,il_tmp, parity
        real(DP)                 :: dl1, dl1min, dl1max, dl2, dl3
        real(DP)                 :: dL_l3l2l1, dL_l1l2l3, dS_l1l2l3, dF_l1l2l3, dG_l1l2l3
        real(DP)                 :: dinfinity, dier
        real(DP)                 :: dspinl2_x1, dspinl3_x1, dspinl2_y1, dspinl3_y1
        real(DP)                 :: dspinl2_x2, dspinl3_x2, dspinl2_y2, dspinl3_y2
        CHARACTER(LEN=13)        :: creturn

        real(DP), allocatable    :: dW_l1l2l3_x1(:), dW_l1l2l3_y1(:)
        real(DP), allocatable    :: dW_l1l2l3_x2(:), dW_l1l2l3_y2(:)

        ! Initialize
        dinfinity = 1.e20_dp
        ALLOCATE(dW_l1l2l3_x1(0:2*ilmax))
        ALLOCATE(dW_l1l2l3_y1(0:2*ilmax))
        ALLOCATE(dW_l1l2l3_x2(0:2*ilmax))
        ALLOCATE(dW_l1l2l3_y2(0:2*ilmax))
        dW_l1l2l3_x1 = 0.0_dp
        dW_l1l2l3_y1 = 0.0_dp
        dW_l1l2l3_x2 = 0.0_dp
        dW_l1l2l3_y2 = 0.0_dp
        dier = 0.0_dp

        dspinl2_x1 = real(ispinl2_x1, kind=DP)
        dspinl3_x1 = real(ispinl3_x1, kind=DP)
        dspinl2_y1 = real(ispinl2_y1, kind=DP)
        dspinl3_y1 = real(ispinl3_y1, kind=DP)

        dspinl2_x2 = real(ispinl2_x2, kind=DP)
        dspinl3_x2 = real(ispinl3_x2, kind=DP)
        dspinl2_y2 = real(ispinl2_y2, kind=DP)
        dspinl3_y2 = real(ispinl3_y2, kind=DP)

        !! Need to check only one flavor (odd*even = 0)
        if (flavor1 == 'tt' .or. flavor1 == 'ee' .or. flavor1 == 'bb') then
            parity = 0
        else if (flavor1 == 'te' .or. flavor1 == 'et') then
            parity = 0
        else if (flavor1 == 'eb' .or. flavor1 == 'be') then
            parity = 1
        else if (flavor1 == 'tb' .or. flavor1 == 'bt') then
            parity = 1
        endif

        do il_tmp=0, il2dim
            il2 = il2range(il_tmp)
            creturn = achar(13)
            WRITE( * , 101 , ADVANCE='NO' ) creturn , il2 , ilmax
            101     FORMAT( a , 'ell : ',i7,' out of a total of ',i7)

            dl2 = real(il2, kind=DP)
            do il3=ilmin, ilmax
                dl3 = real(il3, kind=DP)

                ! Define boundaries for l1 given l2 and l3
                ! il1min = MAX( abs(il2-il3), abs(ispinl2+ispinl3) )
                il1min = MAX( abs(il2-il3), 0 )
                il1max = il2 + il3
                dl1min = real(il1min, kind=DP)
                dl1max = real(il1max, kind=DP)

                ! Number of allowed l1 given l2 and l3
                indim = il1max - il1min + 1

                ! The fortran routine WIG3J compute all allowed l1 given fixed l2,l3,m2,m3.
                ! Therefore we loop over l2,l3 (m2,m3 are fixed and given by the spins), and at each iteration we grab a
                ! vector dW_l1l2l3 of dimension indim, which contains all 3-j symbols, that we will use to update the final Cl
                ! Be careful for the type of variables to pass. See wig3j_f.f for more infos.
                if (flavor1 == 'te' .or. flavor1 == 'et' .or. flavor1 == 'tb' .or. flavor1 == 'bt') then
                    call WIG3J(dl2, dl3, dspinl2_x1, dspinl3_x1, dl1min, dl1max, &
                    & dW_l1l2l3_x1(il1min:il1max), indim, dier)
                    call WIG3J(dl2, dl3, dspinl2_y1, dspinl3_y1, dl1min, dl1max, &
                    & dW_l1l2l3_y1(il1min:il1max), indim, dier)
                else if (flavor1 == 'tt' .or. flavor1 == 'ee' .or. flavor1 == 'bb' &
                    & .or. flavor1 == 'eb' .or. flavor1 == 'be') then
                    ! Includes TT, EE, BB, EB, BE.
                    ! This includes TB as well which get only spin2 wigner.
                    ! use spin_y to include TB as well (spin_x = spin_y for all but TE)
                    call WIG3J(dl2, dl3, dspinl2_y1, dspinl3_y1, dl1min, dl1max, &
                    & dW_l1l2l3_x1(il1min:il1max), indim, dier)
                endif

                if (flavor2 == 'te' .or. flavor2 == 'et' .or. flavor2 == 'tb' .or. flavor2 == 'bt') then
                    call WIG3J(dl2, dl3, dspinl2_x2, dspinl3_x2, dl1min, dl1max, &
                    & dW_l1l2l3_x2(il1min:il1max), indim, dier)
                    call WIG3J(dl2, dl3, dspinl2_y2, dspinl3_y2, dl1min, dl1max, &
                    & dW_l1l2l3_y2(il1min:il1max), indim, dier)
                else if (flavor2 == 'tt' .or. flavor2 == 'ee' .or. flavor2 == 'bb' &
                    & .or. flavor2 == 'eb' .or. flavor2 == 'be') then
                    ! Includes TT, EE, BB, EB, BE.
                    ! This includes TB as well which get only spin2 wigner.
                    ! use spin_y to include TB as well (spin_x = spin_y for all but TE)
                    call WIG3J(dl2, dl3, dspinl2_y2, dspinl3_y2, dl1min, dl1max, &
                    & dW_l1l2l3_x2(il1min:il1max), indim, dier)
                endif

                il1max = MIN( il2 + il3 , ilmax)

                do il1=il1min, il1max
                    ! need to remove first two multipoles from the sum
                    if (il1 .lt. 2) CYCLE
                    dl1 = real(il1, kind=DP)

                    ! For TT, EE, BB and TE cases, the wigner are non-zero only if iL is even
                    iL = il1 + il2 + il3
                    if (MOD(iL,2) .ne. 0 .and. parity .eq. 0) CYCLE
                    ! For the TB and EB cases, the weights are non-zero only if iL is odd
                    if (MOD(iL,2) .ne. 1 .and. parity .eq. 1) CYCLE

                    !! Prefactors
                    dL_l3l2l1 =  dl1*(dl1+1.d0) + dl2*(dl2+1.d0) - dl3*(dl3+1.d0)
                    dL_l1l2l3 = -dl1*(dl1+1.d0) + dl2*(dl2+1.d0) + dl3*(dl3+1.d0)
                    dS_l1l2l3 = SQRT( ( (2*dl1+1.d0)*(2*dl2+1.d0)*(2*dl3+1.d0) )/(16*pi) )

                    !! g factor for simple XX
                    if (flavor1 == 'tt' .or. flavor1 == 'ee' .or. flavor1 == 'bb') then
                        ! Compute function F in Eq 52 or 73 (astro-ph/0001303)
                        dG_l1l2l3 = (dclxx1(0,il1) * dL_l3l2l1 + dclyy1(0,il3) * dL_l1l2l3) * &
                        & dS_l1l2l3 * dW_l1l2l3_x1(il1) / (2*dclxx1(1,il3)*dclyy1(1,il1))

                    else if (flavor1 == 'eb') then
                        ! Compute function F in Eq 52 or 73 (astro-ph/0001303)
                        dG_l1l2l3 = (dclxx1(0,il3) * dL_l1l2l3 - dclyy1(0,il1) * dL_l3l2l1) * &
                        & dS_l1l2l3 * dW_l1l2l3_x1(il1) / (dclxx1(1,il3)*dclyy1(1,il1))

                    else if (flavor1 == 'be') then
                        ! Compute function F in Eq 52 or 73 (astro-ph/0001303)
                        dG_l1l2l3 = (dclyy1(0,il3) * dL_l1l2l3 - dclxx1(0,il1) * dL_l3l2l1) * &
                        & dS_l1l2l3 * dW_l1l2l3_x1(il1) / (dclyy1(1,il3)*dclxx1(1,il1))

                    else if (flavor1 == 'te') then
                        dG_l1l2l3 = (dclxy1(0,il1) * dL_l3l2l1 * dW_l1l2l3_x1(il1) + &
                        &  dclxy1(0,il3) * dL_l1l2l3 * dW_l1l2l3_y1(il1)) * dS_l1l2l3 / &
                        & (dclxx1(1,il3) * dclyy1(1,il1))

                    else if (flavor1 == 'et') then
                        dG_l1l2l3 = (dclxy1(0,il1) * dL_l3l2l1 * dW_l1l2l3_y1(il1) + &
                        &  dclxy1(0,il3) * dL_l1l2l3 * dW_l1l2l3_x1(il1)) * dS_l1l2l3 / &
                        & (dclyy1(1,il3) * dclxx1(1,il1))

                    else if (flavor1 == 'tb') then
                        dG_l1l2l3 = dclxy1(0,il3) * dL_l1l2l3 * dW_l1l2l3_y1(il1) * dS_l1l2l3 / &
                        & (dclxx1(1,il3) * dclyy1(1,il1))

                    else if (flavor1 == 'bt') then
                        dG_l1l2l3 = dclxy1(0,il3) * dL_l1l2l3 * dW_l1l2l3_y1(il1) * dS_l1l2l3 / &
                        & (dclyy1(1,il3) * dclxx1(1,il1))

                    endif

                    !! g factor for simple XX
                    if (flavor2 == 'tt' .or. flavor2 == 'ee' .or. flavor2 == 'bb') then
                        ! Compute function F in Eq 52 or 73 (astro-ph/0001303)
                        dF_l1l2l3 = (dclxx2(0,il1) * dL_l3l2l1 + dclyy2(0,il3) * dL_l1l2l3) * dS_l1l2l3 * dW_l1l2l3_x2(il1)

                    else if (flavor2 == 'eb') then
                        dF_l1l2l3 = (dclxx2(0,il3) * dL_l1l2l3 - dclyy2(0,il1) * dL_l3l2l1) * dS_l1l2l3 * dW_l1l2l3_x2(il1)

                    else if (flavor2 == 'be') then
                        dF_l1l2l3 = (dclyy2(0,il3) * dL_l1l2l3 - dclxx2(0,il1) * dL_l3l2l1) * dS_l1l2l3 * dW_l1l2l3_x2(il1)

                    else if (flavor2 == 'te') then
                        ! Compute function F in Eq 52 or 73 (astro-ph/0001303)
                        dF_l1l2l3 = (dclxy2(0,il1) * dL_l3l2l1 * dW_l1l2l3_x2(il1) + &
                        &            dclxy2(0,il3) * dL_l1l2l3 * dW_l1l2l3_y2(il1)) * dS_l1l2l3

                    else if (flavor2 == 'et') then
                        ! Compute function F in Eq 52 or 73 (astro-ph/0001303)
                        dF_l1l2l3 = (dclxy2(0,il1) * dL_l3l2l1 * dW_l1l2l3_y2(il1) + &
                        &            dclxy2(0,il3) * dL_l1l2l3 * dW_l1l2l3_x2(il1)) * dS_l1l2l3

                    else if (flavor2 == 'tb' .or. flavor2 == 'bt') then
                        dF_l1l2l3 = dclxy2(0,il3) * dL_l1l2l3 * dW_l1l2l3_y2(il1) * dS_l1l2l3

                    endif

                    Bl2l3(il2,il3) = Bl2l3(il2,il3) + dG_l1l2l3 * dF_l1l2l3
                enddo
            enddo
        enddo

        DEALLOCATE(dW_l1l2l3_x1)
        DEALLOCATE(dW_l1l2l3_y1)
        DEALLOCATE(dW_l1l2l3_x2)
        DEALLOCATE(dW_l1l2l3_y2)

    end subroutine

end module
