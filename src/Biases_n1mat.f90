    module LensingBiases
    implicit none
    ! Order 1 2 3 = T E B
    !estimator order TT, EE, EB, TE, TB, BB

    integer, parameter :: n_est = 6
    integer, parameter :: i_TT=1,i_EE=2,i_EB=3,i_TE=4,i_TB=5, i_BB=6
    character(2) :: estnames(n_est) = ['TT','EE','EB','TE','TB','BB']
    integer  :: lumped_indices(2,n_est)=[[1,1],[2,2],[2,3],[1,2],[1,3],[3,3]]
    integer, parameter :: lumpsFor(3,3) = [[i_TT, i_TE, i_TB],[i_TE,i_EE,i_EB],[i_TB, i_EB, i_BB]]

    integer :: renorm_deltaL = 20
    integer, parameter :: dp = kind(1.d0)
    real(dp), parameter :: pi =  3.1415927, twopi=2*pi
    integer, parameter :: lmaxmax = 8000
    real(dp) Norms(lmaxmax,n_est) , NormsCurl(lmaxmax,n_est)
    real(dp) :: CX(lmaxmax), CE(lmaxmax),CB(lmaxmax), CT(lmaxmax),CPhi(lmaxmax)
    real(dp) :: CXf(lmaxmax), CEf(lmaxmax),CBf(lmaxmax), CTf(lmaxmax)
    !real(dp) :: ddCt(lmaxmax),ddCphi(lmaxmax),splineCtphi(lmaxmax),ddCtphi(lmaxmax),ddCX(lmaxmax), ddCE(lmaxmax),ddCB(lmaxmax)
    real(dp) :: NT(lmaxmax), NP(lmaxmax)
    real(dp) :: CEobs(lmaxmax), CTobs(lmaxmax), CBobs(lmaxmax)


    integer Phi_Sample(lmaxmax), nPhiSample
    real(dp) dPhi_Sample(lmaxmax)


    real(dp) :: noise_fwhm_deg != 1.d0/60
    real(dp) :: muKArcmin != 1.0
    real(dp) ::  NoiseVar, NoiseVarP  !=  ( muKArcmin * pi/ 180 / 60.) ** 2, NoiseVarP=NoiseVar*2

    integer :: lmaxout! = 3150 !2150
    integer :: lmax! = 6000 !2050*2
    integer ::lmax_TT! = 6000 !the L_Max for TT, even if higher for pol
    integer :: LMin
    integer :: lmin_filter
    real :: acc=1
    logical :: doCurl = .false.

    logical :: UsePhiSampling = .true. !for N0, otherwise use L step

    integer, parameter :: dL=20, dPhiL=20
    integer, parameter :: Lstep = 20

    ! character(LEN=*), parameter :: dir = './' !'C:\Work\F90\LensingBiases\'
    character(LEN=50) dir
    character(LEN=50) vartag
    character(LEN=:), allocatable :: root
    contains

    subroutine SetPhiSampling(lmx, sampling)
    integer lmx
    integer i,ix, dL,Lix, L
    logical sampling

    print *, 'sampling', sampling

    if (.not. sampling) then
        Lix=0
        do L=LMin, lmx, Lstep
            LIx=Lix+1
            Phi_Sample(Lix)=L
        end do
        nPhiSample = Lix
    else
        ix=0

        do i=2, 110, nint(10/acc)
            ix=ix+1
            Phi_Sample(ix)=i
        end do

        dL =nint(30/acc)
        do i=Phi_Sample(ix)+dL, 580, dL
            ix=ix+1
            Phi_Sample(ix)=i
        end do
        dL =nint(100/acc)
        do i=Phi_Sample(ix)+dL, lmx/2, dL
            ix=ix+1
            Phi_Sample(ix)=i
        end do
        dL =nint(300/acc)
        do i=Phi_Sample(ix)+dL, lmx, dL
            ix=ix+1
            Phi_Sample(ix)=i
        end do

        nPhiSample =  ix
    end if

    dPhi_Sample(1) = (Phi_Sample(2)-Phi_Sample(1))/2.
    do i=2, nPhiSample-1
        dPhi_Sample(i) = (Phi_Sample(i+1)-Phi_Sample(i-1))/2.
    end do
    dPhi_Sample(nPhiSample) = (Phi_Sample(nPhiSample)-Phi_Sample(nPhiSample-1))

    end subroutine SetPhiSampling

    subroutine NoiseInit(AN,ANP)
    real(dp) xlc, sigma2, AN, ANP
    integer l

    xlc= 180*sqrt(8.*log(2.))/pi
    sigma2 = (noise_fwhm_deg/xlc)**2
    do l=2, lmax
        NT(L) = AN*exp(l*(l+1)*sigma2)
        if (l>lmax_TT) NT(L) = NT(L) + ( 0.000001*pi/180/60.)**2 *exp(l*(l+1)*(15./60./xlc)**2)
        NP(L) = ANP*exp(l*(l+1)*sigma2)
    end do

    end subroutine NoiseInit

    subroutine ReadPhiPhi(Filename)
    character(LEN=*) filename
    integer file_id
    character(LEN=1024) InLine
    integer L, status
    real(dp) T, E, B, TE, phi

    open(file=Filename, newunit = file_id, form='formatted', status='old', iostat=status)
    if (status/=0) stop 'error opening Cl'
    CPhi=0
    do
        read(file_id, '(a)', iostat=status) InLine
        if (InLine=='') cycle
        if (InLIne(1:1)=='#') cycle
        if (status/=0) exit
        read(InLine,*, iostat=status)  l, T, E, B , TE, phi
        if (status/=0) then
            read(InLine,*, iostat=status)  l, phi
        end if
        if (L<1) cycle
        if (L> Lmax) exit
        CPhi(L) = phi * twopi/real(L*(L+1),dp)**2 !*(L/100.)**(-0.2) !!!!!
    end do
    close(file_id)

    end subroutine ReadPhiPhi

    subroutine ReadPower(Filename)
    character(LEN=*) filename
    integer file_id
    character(LEN=1024) InLine
    integer L, status
    real(dp) T, E, B, TE, Phi
    logical :: newform = .false.

    open(file=Filename, newunit = file_id, form='formatted', status='old', iostat=status)
    if (status/=0) stop 'error opening Cl'
    CT=0
    CE=0
    CB=0
    CX=0
    do
        read(file_id, '(a)', iostat=status) InLine
        if (status/=0) exit
        if (InLine=='') cycle
        if (InLIne(1:1)=='#') then
            !!  newform = .true.
            cycle
        end if
        if (newform) then
            read(InLine,*, iostat=status)  l, T, E, TE, Phi
            B=0
            if (status/=0) stop 'error reading power'
        else
            read(InLine,*, iostat=status)  l, T, E, B , TE
            if (status/=0) then
                B=1
                read(InLine,*, iostat=status)  l, T, E, TE
            end if
        end if
        if (L> Lmax) exit
        if (newform .and. L>=1) then
            CPhi(L) = phi * twopi/real(L*(L+1),dp)**2
        end if
        if (L<2) cycle
        CT(L) = T * twopi/(l*(l+1))
        CE(L) = E * twopi/(l*(l+1))
        CB(L) = B * twopi/(l*(l+1))
        CX(L) = TE * twopi/(l*(l+1))
    end do
    CTf=CT
    CEf=CE
    CBf=CB
    CXf=CX
    do L=2,lmax
        !        CT(L) = CT(L) * (l/1500.d0)**0.01
    end do
    close(file_id)

    end subroutine ReadPower

    subroutine ReadFilter(Filename)
    character(LEN=*) filename
    integer file_id
    character(LEN=1024) InLine
    integer L, status
    real(dp) T

    open(file=Filename, newunit = file_id, form='formatted', status='old', iostat=status)
    if (status/=0) stop 'error opening Cl'
    CTobs=0
    do
        read(file_id, '(a)', iostat=status) InLine
        if (status/=0) exit
        if (InLine(1:1)=='#') cycle
        read(InLine,*, iostat=status) L, T
        if (L<2) cycle
        if (T==0) then
            lmin_filter = L+1
        else
            CTobs(L) = 1/T
        end if
        if (L >= Lmax) exit
    end do
    print *,'Lmin,Lmax filter = ', lmin_filter, L
    lmax = min(lmax,L)
    close(file_id)

    end subroutine ReadFilter


    subroutine ReadRuthFile(FileName)
    character(LEN=*) filename
    integer file_id
    character(LEN=1024) InLine
    integer L, status
    real(dp) E, B, Etot, Btot,  phi

    open(file=Filename, newunit = file_id, form='formatted', status='old', iostat=status)
    if (status/=0) stop 'error opening Ruth Cl'
    CT=1
    CE=0
    CB=0
    CEobs=0
    CBobs=0
    CTobs=1
    CX=1
    Phi=0
    do
        read(file_id, '(a)', iostat=status) InLine
        if (status/=0) exit
        read(InLine,*)  l, E, B , Etot, Btot, phi
        if (l<2) cycle
        if (L> Lmax) exit
        CE(L) = E
        CB(L) = B
        CEobs(L)=Etot
        CBobs(L)=Btot
        CPhi(L) = phi
        !        if (L<=40) CPhi(L)=0
    end do
    close(file_id)

    end subroutine ReadRuthFile

    subroutine WriteMatrixLine(file_id,N0)
    integer i,j, file_id
    real(dp) N0(:,:)

    do i=1, n_est
        do j=1,n_est
            write (file_id,'(1E16.6)',Advance='NO') N0(i,j)
        end do
    end do
    write(file_id,'(a)') ''

    end subroutine WriteMatrixLine

    subroutine getNorm(WantIntONly,sampling)
    integer L, Lix, l1, nphi, phiIx, L2int
    real(dp) dphi
    real(dp) L1Vec(2), L2vec(2), LVec(2)
    real(dp) phi, cos2L1L2, sin2
    real(dP) norm, L2
    real(dp) cosfac,f12(n_est),f21(n_est),Win21(n_est),Win12(n_est), fac
    integer, parameter :: dL1=1
    integer file_id, i,j, icurl
    logical isCurl
    logical WantIntONly
    logical sampling
    real(dp) N0(n_est,n_est), N0_L(n_est,n_est)

    call SetPhiSampling(lmaxout, sampling)

    !Order TT, EE, EB, TE, TB, BB
    do icurl = 0,1
        isCurl = icurl==1
        if (IsCurl) then
            if (.not. doCurl) cycle
            open(file=trim(dir)//'N0'//trim(vartag)//'_Curl.dat', newunit = file_id, form='formatted', status='replace')
        else
            open(file=trim(dir)//'N0'//trim(vartag)//'.dat', newunit = file_id, form='formatted', status='replace')
        end if


        do Lix =1, nPhiSample
            L = Phi_Sample(Lix)
            Lvec(1) = L
            LVec(2)= 0
            N0_L=0
            !$OMP PARALLEL DO default(shared), private(L1,nphi,dphi,N0,PhiIx,phi,L1vec,L2vec,L2,L2int,cos2L1L2,sin2,cosfac,f12,f21,Win12,Win21,i,j), &
            !$OMP reduction(+:N0_L)
            do L1=lmin_filter, lmax, dL1

            nphi=(2*L1+1)
            dphi=(2*Pi/nphi)
            N0=0

            do phiIx=0,(nphi-1)/2
                !                  do phiIx=-(nphi-1)/2, (nphi-1)/2
                phi= dphi*PhiIx
                L1vec(1)=L1*cos(phi)
                L1vec(2)=L1*sin(phi)
                L2vec = Lvec-L1vec
                L2=(sqrt(L2vec(1)**2+L2vec(2)**2))
                if (L2<lmin_filter .or. L2>lmax) cycle
                L2int=nint(L2)

                if (isCurl) then
                    call getResponseFull(L1vec(2)*L,L2vec(2)*L, L1vec,real(L1,dp),L1, L2vec,L2, L2int,f12,f21, CT, CE, CX, CB)
                    call getWins(L1vec(2)*L,L2vec(2)*L, L1vec,real(L1,dp),L1, L2vec,L2, L2int,  Win12, Win21)
                else
                    call getResponseFull(L1vec(1)*L,L2vec(1)*L, L1vec,real(L1,dp),L1, L2vec,L2, L2int,f12,f21, CT, CE, CX, CB)
                    call getWins(L1vec(1)*L,L2vec(1)*L, L1vec,real(L1,dp),L1, L2vec,L2, L2int,  Win12, Win21)
                end if

                do i=1,n_est
                    N0(i,i) = N0(i,i) + f12(i)*Win12(i)
                end do


                !Important to use symmetric form here if only doing half PhiIx integral
                N0(i_TT,i_TE) = N0(i_TT,i_TE) + Win12(i_TT)*(Win12(i_TE)*CTobs(L1)*CXf(L2int) + Win21(i_TE)*CXf(L1)*CTobs(L2int))

                N0(i_TT,i_EE) = N0(i_TT,i_EE) + 2*Win12(i_TT)*Win12(i_EE)*CXf(L1)*CXf(L2int)

                N0(i_EE,i_TE) = N0(i_EE,i_TE) + Win12(i_EE)*(Win12(i_TE)*CXf(L1)*CEobs(L2int) + Win21(i_TE)*CEobs(L1)*CXf(L2int))

                N0(i_EB,i_TB) = N0(i_EB,i_TB) + (Win12(i_EB)*(Win12(i_TB)*CXf(L1)*CBobs(L2int)) + &
                    Win21(i_EB)*(Win21(i_TB)*CXf(L2int)*CBobs(L1)))/2

                if (PhiIx==0) N0 = N0/2
            end do
            fac = dphi* L1*dL1 *2
            N0_L = N0_L + N0 * fac

        end do
        !$OMP END PARALLEL DO
        N0_L = N0_L/(twopi**2)

        N0=0
        do i=1,n_est
            do j=i,n_est
                !N0_TT_TE(L)= Norm_T(L)*Norm_X(L)*N0_L(i_TT, i_TE)
                if (WantIntONly) then
                    N0(i,j) = N0_L(i,j)
                else
                    N0(i,j) = N0_L(i,j)/N0_L(i,i)/N0_L(j,j)
                end if
                N0(j,i) = N0(i,j)
            end do
        end do

        do i=1,n_est
            Norms(L,i) = N0(i,i)
        end do

        norm = real(L*(L+1),dp)**2/twopi
        write (file_id,'(1I5, 1E16.6)',Advance='NO') L, CPhi(L)*norm

        call WriteMatrixLine(file_id,N0)


        !write (file_id,'(1I5, 6E16.6)') L, norm_T(L),norm_E(L), norm_EB(L), norm_X(L),N0_TT_TE(L), CPhi(L)*norm
        !        print *, L,  norm*Norm_T(L),norm*Norm_E(L),norm*Norm_EB(L),norm*Norm_X(L)
        print *, L, CPhi(L)*norm, N0(i_TT,i_TT)*norm
    end do
    close(file_id)
    end do

    end subroutine getNorm

    subroutine loadNorm()
    integer ell, file_id, L,i
    real(dp) N0(n_est,n_est), dum


    open(file=trim(dir)//'N0'//trim(vartag)//'.dat', newunit = file_id, form='formatted', status='old')
    do L=LMin, lmaxout, Lstep
        read(file_id,*) ell, dum, N0
        if (L/=ell) stop 'wrong N0 file'
        do i=1,n_est
            Norms(L,i) = N0(i,i)
        end do
    end do
    close(file_id)
    end subroutine loadNorm

    subroutine loadNormCurl()
    integer ell, file_id, L,i
    real(dp) N0(n_est,n_est), dum


    open(file=trim(dir)//'N0'//trim(vartag)//'_Curl.dat', newunit = file_id, form='formatted', status='old')
    do L=LMin, lmaxout, Lstep
        read(file_id,*) ell, dum, N0
        if (L/=ell) stop 'wrong N0 file'
        do i=1,n_est
            NormsCurl(L,i) = N0(i,i)
        end do
    end do
    close(file_id)
    end subroutine loadNormCurl

    subroutine WriteRanges(name)
    integer PhiLix, L, file_id
    character(LEN=*), intent(in) :: name

    open(file=trim(dir)//trim(name)//trim(vartag)//'_Lin.dat', newunit = file_id, form='formatted', status='replace')
    do PhiLix=1, nPhiSample
        write (file_id, '(1I6,1E16.6)') Phi_Sample(PhiLix), dPhi_Sample(PhiLix)
    end do
    close(file_id)

    open(file=trim(dir)//trim(name)//trim(vartag)//'_Lout.dat', newunit = file_id, form='formatted', status='replace')
    do L=LMin, lmaxout, Lstep
        write (file_id, *) L
    end do
    close(file_id)

    end subroutine WriteRanges

    subroutine WriteMatrix(name, matrix)
    integer Lix, L, file_id, PhiLix, L1
    character(LEN=*), intent(in) :: name
    real(dp), intent(in) :: matrix(:,:)

    Lix=0
    open(file=trim(dir)//trim(name)//trim(vartag)//'_matrix.dat', newunit = file_id, form='formatted', status='replace')

    write (file_id, '(1I6)', advance='NO') 0
    do PhiLix=1, nPhiSample
        write (file_id, '(1I16)', advance='NO') Phi_Sample(PhiLIx)
    end do
    write(file_id,'(a)') ''
    do L=LMin, lmaxout, Lstep
        Lix=Lix+1
        write (file_id, '(1I6)', advance='NO') L
        do PhiLix=1, nPhiSample
            L1 = Phi_Sample(PhiLIx)
            write (file_id, '(1E16.6)', advance='NO') matrix(Lix,PhiLix)/((L1*(L1+1.d0))**2/twopi) * (L*(L+1.d0))**2/twopi
        end do
        write(file_id,'(a)') ''
    end do
    close(file_id)

    end subroutine WriteMatrix


    subroutine GetN1Temp()
    integer L, l1, nphi, phiIx, PhiL_nphi, PhiL_phi_ix, L2int, L3int, L4int
    integer PhiL
    real(dp) dphi, PhiL_phi_dphi
    real(dp) L1Vec(2), L2vec(2), LVec(2), L3Vec(2),L4Vec(2), phiLVec(2)
    real(dp) N1, phi, PhiL_phi, N1_L1, N1_PhiL
    real(dP) norm, L2, L4, L3
    real(dp) g12, g34, tmp, fact, dPh
    real(dp), allocatable :: Matrix(:,:), MatrixL1(:)
    integer file_id,  Lix, PhiLix
    character(LEN=12) TTtag

    call loadNorm()
    TTtag='N1_TT'
    call WriteRanges(TTtag)
    open(file=trim(dir)//trim(TTtag)//trim(vartag)//'.dat', newunit = file_id, form='formatted', status='replace')

    allocate(matrix((lmaxout-Lmin)/Lstep+1,nPhiSample))
    allocate(matrixL1(nPhiSample))

    matrix=0
    Lix=0
    do L=LMin, lmaxout, Lstep
        Lix=Lix+1
        Lvec(1) = L
        LVec(2)= 0
        N1=0
        do L1=max(lmin_filter,dL/2), lmax, dL
            N1_L1 = 0
            matrixL1=0
            nphi=(2*L1+1)
            if (L1>3*dL) nphi=2*nint(L1/real(2*dL))+1
            dphi=(2*Pi/nphi)
            !$OMP PARALLEL DO default(shared), private(PhiIx,phi, L1vec,L2vec, L2, L2int, g12, tmp, L3, L3int, PhiL_nphi, PhiL_phi_dphi, PhiL_phi_ix, PhiL_phi, L3vec, L4vec, L4), &
            !$OMP private(L4int,  g34, PhiL, PhiLix, PhiLvec, fact,N1_PhiL, dPh), schedule(STATIC), reduction(+:N1_L1)
            do phiIx=0,(nphi-1)/2 ! -(nphi-1)/2, (nphi-1)/2
                phi= dphi*PhiIx
                L1vec(1)=L1*cos(phi)
                L1vec(2)=L1*sin(phi)
                L2vec = Lvec-L1vec
                L2=(sqrt(L2vec(1)**2+L2vec(2)**2))
                if (L2<lmin_filter .or. L2>lmax) cycle
                L2int=nint(L2)
                g12= (L1vec(1)*L*CTf(L1) + L2vec(1)*L*CTf(L2int))/(2*CTobs(L1)*CTobs(L2int))
                N1_PhiL=0
                do PhiLIx = 1, nPhiSample
                    PhiL = Phi_Sample(PhiLIx)
                    dPh = dPhi_Sample(PhiLIx)
                    PhiL_nphi=(2*PhiL+1)
                    if (phiL>20) PhiL_nphi=2*nint(real(PhiL_nphi)/dPh/2)+1
                    PhiL_phi_dphi=(2*Pi/PhiL_nphi)
                    tmp=0
                    do PhiL_phi_ix=-(PhiL_nphi-1)/2, (PhiL_nphi-1)/2
                        PhiL_phi= PhiL_phi_dphi*PhiL_phi_ix
                        PhiLvec(1)=PhiL*cos(PhiL_phi)
                        PhiLvec(2)=PhiL*sin(PhiL_phi)

                        L3vec= PhiLvec - L1vec
                        L3 = sqrt(L3vec(1)**2+L3vec(2)**2)
                        if (L3>=lmin_filter .and. L3<=lmax) then
                            L3int = nint(L3)
                            L4vec = -Lvec-L3vec
                            L4 = sqrt(L4vec(1)**2+L4vec(2)**2)
                            if (L4>=lmin_filter .and. L4<=lmax) then
                                L4int=nint(L4)
                                g34= (L3vec(1)*CTf(L3int) + L4vec(1)*CTf(L4int))/(2*CTobs(L3int)*CTobs(L4int))
                                tmp = tmp -g34*(dot_product(L1vec,PhiLvec)*CT(L1)+dot_product(L3vec,PhiLvec)*CT(L3int))&
                                    *(dot_product(L2Vec,PhiLvec)*CT(L2int)+dot_product(L4Vec,PhiLvec)*CT(L4int))
                            end if
                        end if

                    end do
                    if (phiIx/=0) tmp=tmp*2 !integrate 0-Pi for phi_L1
                    fact=tmp* PhiL_phi_dphi* L*PhiL
                    !$OMP ATOMIC
                    matrixL1(phiLix)=matrixL1(phiLix)+ fact*g12*dPh
                    N1_PhiL= N1_PhiL + fact*Cphi(PhiL)*dPh
                end do
                N1_L1 =N1_L1 + N1_PhiL*g12
            end do
            !$OMP END PARALLEL DO
            !            print *, L1, N1_L1* dphi* L1
            matrix(Lix,:)=matrix(Lix,:) + matrixL1(:)*dphi*L1*dL
            N1= N1 + N1_L1 * dphi* L1*dL
        end do
        fact = -2*Norms(L,i_TT)**2  / (twopi**4)
        matrix(Lix,:)=matrix(Lix,:)*fact
        N1 = N1*fact

        !Get normalized and unnormalized
        write (file_id,'(1I5, 2E15.5)') L, N1,  N1/Norms(L,i_TT)**2   !dot_product(matrix(Lix,:), CPhi(Phi_Sample(1:nPhiSample))*dPhi_Sample(1:nPhiSample))
        norm = real(L*(L+1),dp)**2/twopi
        print *, '***', L,  norm*N1, N1
    end do
    close(file_id)

    call WriteMatrix(TTtag, matrix)

    end subroutine GetN1Temp

    subroutine getResponseFull(L_dot_L1,L_dot_L2, L1vec,L1,L1int, L2vec,L2, L2int,f12,f21, CT, CE, CX, CB)
    integer, intent(in) :: L1int, L2int
    real(dp), intent(in) :: L1vec(2),L2vec(2), L1,L2, L_dot_L1,L_dot_L2
    real(dp) cosfac, sin2, cos2L1L2
    real(dp) f12(n_est), f21(n_est)
    real(dp) CT(:), CE(:), CX(:), CB(:)

    f12(i_TT)= (L_dot_L1*CT(L1) + L_dot_L2*CT(L2int))

    cosfac= dot_product(L1vec,L2vec)/real(L1*L2,dp)
    cos2L1L2 =2*cosfac**2-1
    f12(i_EE)= (L_dot_L1*CE(L1int) + L_dot_L2*CE(L2int))*cos2L1L2

    if (n_est>=i_BB) then
        f12(i_BB) = (L_dot_L1*CB(L1) + L_dot_L2*CB(L2int))*cos2L1L2
    end if
    f21=f12

    sin2=  2*cosfac*(L1vec(2)*L2vec(1)-L1vec(1)*L2vec(2))/(L2*L1)
    f12(i_EB) = (L_dot_L1*CE(L1int) - L_dot_L2*CB(L2int))*sin2
    f21(i_EB) = -(L_dot_L2*CE(L2int) - L_dot_L1*CB(L1int))*sin2

    f12(i_TE)=  L_dot_L1*CX(L1int)*cos2L1L2 + L_dot_L2*CX(L2int)
    f21(i_TE)=  L_dot_L2*CX(L2int)*cos2L1L2 + L_dot_L1*CX(L1int)

    f12(i_TB) = L_dot_L1*CX(L1int)*sin2
    f21(i_TB) = -L_dot_L2*CX(L2int)*sin2

    end subroutine getResponseFull

    subroutine getResponsefid(L_dot_L1,L_dot_L2, L1vec,L1,L1int, L2vec,L2, L2int,f12,f21)
    integer, intent(in) :: L1int, L2int
    real(dp), intent(in) :: L1vec(2),L2vec(2), L1,L2, L_dot_L1,L_dot_L2
    real(dp) f12(n_est), f21(n_est)

    call getResponseFull(L_dot_L1,L_dot_L2, L1vec,L1,L1int, L2vec,L2, L2int,f12,f21, CTf, CEf, CXf, CBf)

    end subroutine getResponsefid

    subroutine getResponse(L_dot_L1,L_dot_L2, L1vec,L1,L1int, L2vec,L2, L2int,f12,f21)
    integer, intent(in) :: L1int, L2int
    real(dp), intent(in) :: L1vec(2),L2vec(2), L1,L2, L_dot_L1,L_dot_L2
    real(dp) f12(n_est), f21(n_est)

    call getResponseFull(L_dot_L1,L_dot_L2, L1vec,L1,L1int, L2vec,L2, L2int,f12,f21, CT, CE, CX, CB)

    end subroutine getResponse


    subroutine getWins(L_dot_L1,L_dot_L2, L1vec,L1,L1int, L2vec,L2, L2int, Win12, Win21)
    integer, intent(in) :: L1int, L2int
    real(dp), intent(in) :: L1vec(2),L2vec(2), L1,L2, L_dot_L1,L_dot_L2
    real(dp) f12(n_est), f21(n_est)
    real(dp), intent(out) :: Win12(n_est)
    real(dp), intent(out), optional :: Win21(n_est)

    call getResponsefid(L_dot_L1,L_dot_L2, L1vec,L1,L1int, L2vec,L2, L2int, f12, f21)

    Win12(i_TT) = f12(i_TT)/(2*CTobs(L1int)*CTobs(L2int))

    Win12(i_EE) = f12(i_EE)/(2*CEobs(L1int)*CEobs(L2int))

    Win12(i_EB) = f12(i_EB)/(CEobs(L1int)*CBobs(L2int))

    Win12(i_TE) = (f12(i_TE)* CEobs(L1int)*CTobs(L2int) - f21(i_TE)*CX(L1int)*CX(L2int))&
        /(CTobs(L1int)*CEobs(L2int)*CTobs(L2int)*CEobs(L1int) - (CX(L1int)*CX(L2int))**2)

    Win12(i_TB) = f12(i_TB)/(CTobs(L1int)*CBobs(L2int))

    if (n_est>=i_BB) then
        Win12(i_BB) = f12(i_BB)/(2*CBobs(L1int)*CBobs(L2int))
    end if

    if (present(Win21)) then
        Win21=Win12

        Win21(i_TE) = (f21(i_TE)* CTobs(L1int)*CEobs(L2int) - f12(i_TE)*CX(L1int)*CX(L2int))&
            /(CTobs(L1int)*CEobs(L2int)*CTobs(L2int)*CEobs(L1int) - (CX(L1int)*CX(L2int))**2)

        Win21(i_EB) = f21(i_EB)/(CEobs(L2int)*CBobs(L1int))
        Win21(i_TB) = f21(i_TB)/(CTobs(L2int)*CBobs(L1int))
    end if

    end subroutine getWins

    function responseFor(i,j, f12,f21)
    integer, intent(in) :: i,j
    real(dp) f12(n_est), f21(n_est)
    integer ix
    real(dp) responseFor

    if (j>=i) then
        ix= lumpsFor(i,j)
        if (ix<=n_est) then
            responseFor = f12(ix)
        else
            responseFor=0
        end if
    else
        ix= lumpsFor(j,i)
        if (ix<=n_est) then
            responseFor = f21(ix)
        else
            responseFor=0
        end if
    end if

    end function responseFor

    subroutine GetN1General()
    !general form
    integer L, l1, nphi, phiIx, PhiL_nphi, PhiL_phi_ix, L2int, L3int, L4int
    integer PhiL
    real(dp) dphi, PhiL_phi_dphi
    real(dp) L1Vec(2), L2vec(2), LVec(2), L3Vec(2),L4Vec(2), phiLVec(2)
    real(dp) phi, PhiL_phi
    real(dP) L2, L4, L3
    real(dp) dPh
    real(dp) phiL_dot_L2, phiL_dot_L3, phiL_dot_L1, phiL_dot_L4
    real(dp) fact(n_est,n_est),tmp(n_est,n_est), N1(n_est,n_est), N1_L1(n_est,n_est),N1_PhiL(n_est,n_est)
    real(dp) Win12(n_est), Win34(n_est), Win43(n_est)
    real(dp) WinCurl12(n_est), WinCurl34(n_est), WinCurl43(n_est), tmpCurl(n_est,n_est), &
        factCurl(n_est,n_est),N1_PhiL_Curl(n_est,n_est), N1_L1_Curl(n_est,n_est),  N1_Curl(n_est,n_est)
    real(dp) f24(n_est), f13(n_est),f31(n_est), f42(n_est)
    integer file_id, file_id_Curl, PhiLix, Lix
    integer ij(2),pq(2), est1, est2
    real(dp) tmpPS, tmpPSCurl, N1_PhiL_PS, N1_PhiL_PS_Curl, N1_L1_PS_Curl, N1_L1_PS, N1_PS, N1_PS_Curl
    integer file_id_PS
    real(dp) this13, this24
    character(LEN=10) outtag

    outtag = 'N1_All'
    call loadNorm()
    call loadNormCurl()

    call WriteRanges(outtag)
    open(file=trim(dir)//trim(outtag)//trim(vartag)//'.dat', newunit = file_id, form='formatted', status='replace')
    open(file=trim(dir)//trim(outtag)//trim(vartag)//'_Curl.dat', newunit = file_id_Curl, form='formatted', status='replace')
    open(file=trim(dir)//trim(outtag)//trim(vartag)//'_PS.dat', newunit = file_id_PS, form='formatted', status='replace')

    Lix=0
    do L=LMin, lmaxout, Lstep
        Lix=Lix+1
        Lvec(1) = L
        LVec(2)= 0
        N1=0
        N1_Curl = 0
        N1_PS=0
        N1_PS_Curl=0
        do L1=max(lmin_filter,dL/2), lmax, dL
            N1_L1 = 0
            N1_L1_Curl = 0
            N1_L1_PS = 0
            N1_L1_PS_Curl = 0


            nphi=(2*L1+1)
            if (L1>3*dL) nphi=2*nint(L1/real(2*dL))+1
            dphi=(2*Pi/nphi)

            !$OMP PARALLEL DO default(shared), private(PhiIx,phi,PhiL_nphi, PhiL_phi_dphi, PhiL_phi_ix, PhiL_phi,PhiLix, dPh), &
            !$OMP private(L1vec,L2,L2vec, L2int,  L3, L3vec, L3int, L4, L4vec, L4int),&
            !$OMP private(tmp, Win12, Win34, Win43, fact,phiL_dot_L1, phiL_dot_L2, phiL_dot_L3, phiL_dot_L4), &
            !$OMP private(tmpCurl, WinCurl12, WinCurl34, WinCurl43, factCurl, N1_PhiL_Curl), &
            !$OMP private(f13, f31, f42, f24, ij, pq, est1, est2,this13,this24), &
            !$OMP private(tmpPS, tmpPSCurl, N1_PhiL_PS, N1_PhiL_PS_Curl), &

            !$OMP private(PhiL, PhiLVec, N1_PhiL), schedule(STATIC), reduction(+:N1_L1), reduction(+:N1_L1_Curl), &
            !$OMP reduction(+:N1_L1_PS), reduction(+:N1_L1_PS_Curl)
            !do phiIx= -(nphi-1)/2, (nphi-1)/2
            do phiIx=0,(nphi-1)/2 !
            phi= dphi*PhiIx
            L1vec(1)=L1*cos(phi)
            L1vec(2)=L1*sin(phi)
            L2vec = Lvec-L1vec
            L2=(sqrt(L2vec(1)**2+L2vec(2)**2))
            if (L2<lmin_filter .or. L2>lmax) cycle
            L2int=nint(L2)

            call getWins(L*L1vec(1),L*L2vec(1), L1vec,real(L1,dp),L1, L2vec,L2, L2int,  Win12)
            call getWins(L*L1vec(2),L*L2vec(2), L1vec,real(L1,dp),L1, L2vec,L2, L2int,  WinCurl12)

            N1_PhiL=0
            N1_PhiL_Curl=0
            N1_PhiL_PS=0
            N1_PhiL_PS_Curl=0
            do PhiLIx = 1, nPhiSample
                PhiL = Phi_Sample(PhiLIx)
                dPh = dPhi_Sample(PhiLIx)
                PhiL_nphi=(2*PhiL+1)
                if (phiL>20) PhiL_nphi=2*nint(real(PhiL_nphi)/dPh/2)+1
                PhiL_phi_dphi=(2*Pi/PhiL_nphi)
                tmp=0
                tmpCurl=0
                tmpPS = 0
                tmpPSCurl = 0
                do PhiL_phi_ix=-(PhiL_nphi-1)/2, (PhiL_nphi-1)/2
                    PhiL_phi= PhiL_phi_dphi*PhiL_phi_ix
                    PhiLvec(1)=PhiL*cos(PhiL_phi)
                    PhiLvec(2)=PhiL*sin(PhiL_phi)
                    L3vec= PhiLvec - L1vec
                    L3 = sqrt(L3vec(1)**2+L3vec(2)**2)
                    if (L3>=lmin_filter .and. L3<=lmax) then
                        L3int = nint(L3)

                        L4vec = -Lvec-L3vec
                        L4 = sqrt(L4vec(1)**2+L4vec(2)**2)
                        L4int=nint(L4)
                        if (L4>=lmin_filter .and. L4<=lmax) then
                            call getWins(-L*L3vec(1),-L*L4vec(1), L3vec,L3,L3int, L4vec,L4, L4int,  Win34, Win43)
                            call getWins(-L*L3vec(2),-L*L4vec(2), L3vec,L3,L3int, L4vec,L4, L4int,  WinCurl34, WinCurl43)

                            phiL_dot_L1=dot_product(PhiLVec,L1vec)
                            phiL_dot_L2=-dot_product(PhiLVec,L2vec)
                            phiL_dot_L3=dot_product(PhiLVec,L3vec)
                            phiL_dot_L4=-dot_product(PhiLVec,L4vec)

                            call getResponse(phiL_dot_L1,phiL_dot_L3, L1vec,real(L1,dp),L1, L3vec,L3, L3int,  f13, f31)
                            call getResponse(phiL_dot_L2,phiL_dot_L4, L2vec,L2,L2int, L4vec,L4, L4int,  f24, f42)

                            do est1=1,n_est
                                ij=lumped_indices(:,est1)
                                do est2=est1,n_est
                                    pq=lumped_indices(:,est2)
                                    this13 = responseFor(ij(1),pq(1),f13,f31)
                                    this24 = responseFor(ij(2),pq(2),f24,f42)
                                    tmp(est1,est2)=tmp(est1,est2)+this13*this24*Win34(est2)
                                    tmpCurl(est1,est2)=tmpCurl(est1,est2)+this13*this24*WinCurl34(est2)

                                    this13 = responseFor(ij(1),pq(2),f13,f31)
                                    this24 = responseFor(ij(2),pq(1),f24,f42)
                                    tmp(est1,est2)=tmp(est1,est2)+this13*this24*Win43(est2)
                                    tmpCurl(est1,est2)=tmpCurl(est1,est2)+this13*this24*WinCurl43(est2)
                                end do
                            end do
                            tmpPS = tmpPS + Win43(1) + Win34(1)
                            tmpPSCurl = tmpPSCurl + WinCurl43(1) + WinCurl34(1)
                        end if
                    end if
                end do
                if (phiIx/=0) tmp=tmp*2 !integrate 0-Pi for phi_L1
                if (phiIx/=0) tmpCurl=tmpCurl*2 !integrate 0-Pi for phi_L1
                if (phiIx/=0) tmpPS=tmpPS*2 !integrate 0-Pi for phi_L1
                if (phiIx/=0) tmpPSCurl=tmpPSCurl*2 !integrate 0-Pi for phi_L1
                fact = tmp* PhiL_phi_dphi* PhiL
                factCurl = tmpCurl* PhiL_phi_dphi* PhiL
                N1_PhiL= N1_PhiL + fact * Cphi(PhiL)*dPh
                N1_PhiL_Curl= N1_PhiL_Curl + factCurl * Cphi(PhiL)*dPh

                N1_PhiL_PS = N1_PhiL_PS + tmpPS* PhiL_phi_dphi* PhiL*dPh  !/ PhiL**2
                N1_PhiL_PS_Curl = N1_PhiL_PS_CUrl + tmpPSCurl* PhiL_phi_dphi* PhiL*dPh  !/ PhiL**2

            end do
            do est1=1,n_est
                N1_PhiL(est1,:)=N1_PhiL(est1,:)*Win12(est1)
                N1_PhiL_Curl(est1,:)=N1_PhiL_Curl(est1,:)*WinCurl12(est1)
            end do
            N1_L1 = N1_L1+N1_PhiL
            N1_L1_Curl = N1_L1_Curl+N1_PhiL_Curl
            N1_L1_PS = N1_L1_PS + N1_PhiL_PS * Win12(1)
            N1_L1_PS_Curl = N1_L1_PS_Curl + N1_PhiL_PS_Curl * WinCurl12(1)
        end do
        !$OMP END PARALLEL DO
        N1= N1 + N1_L1 * dphi* L1*dL
        N1_Curl= N1_Curl + N1_L1_Curl * dphi* L1*dL
        N1_PS = N1_PS + N1_L1_PS *dphi*L1*dL
        N1_PS_Curl = N1_PS_Curl + N1_L1_PS_Curl *dphi*L1*dL

    end do

    do est1=1,n_est
        do est2=est1,n_est
            N1(est1,est2) = norms(L,est1)*norms(L,est2)*N1(est1,est2) / (twopi**4)
            N1(est2,est1) = N1(est1,est2)
            N1_Curl(est1,est2) = normsCurl(L,est1)*normsCurl(L,est2)*N1_Curl(est1,est2) / (twopi**4)
            N1_Curl(est2,est1) = N1_Curl(est1,est2)
        end do
    end do
    N1_PS = norms(L,1)*norms(L,1)*N1_PS / (twopi**4)
    N1_PS_Curl = normsCurl(L,1)*normsCurl(L,1)*N1_PS_Curl / (twopi**4)

    write(file_id,'(1I5)',advance='NO') L
    call WriteMatrixLine(file_id, N1)
    write(file_id_Curl,'(1I5)',advance='NO') L
    call WriteMatrixLine(file_id_Curl, N1_Curl)

    write(file_id_PS,*) L, N1_PS, N1_PS_Curl


    print *, 'Phi',L, N1(i_TT,i_TT), N1(i_tb,i_eb), N1(i_eb, i_ee)
    print *, 'Psi',L, N1_Curl(i_TT,i_TT), N1_Curl(i_tb,i_eb), N1_Curl(i_eb, i_ee)
    print *, 'PS', L, N1_PS, N1_PS_Curl

    !norm = real(L*(L+1),dp)**2/twopi
    !       print *, '***', L,  norm*N1_EEEB, N1_EEEB !, dot_product(matrix(Lix,:), CPhi(Phi_Sample(1:nPhiSample))*dPhi_Sample(1:nPhiSample))
    end do
    close(file_id)
    close(file_id_Curl)
    close(file_id_PS)

    end subroutine GetN1General


    subroutine GetN1MatrixGeneral()
    !general form
    integer L, l1, nphi, phiIx, PhiL_nphi, PhiL_phi_ix, L2int, L3int, L4int
    integer PhiL
    real(dp) dphi, PhiL_phi_dphi
    real(dp) L1Vec(2), L2vec(2), LVec(2), L3Vec(2),L4Vec(2), phiLVec(2)
    real(dp) phi, PhiL_phi
    real(dP) L2, L4, L3
    real(dp) dPh
    real(dp) phiL_dot_L2, phiL_dot_L3, phiL_dot_L1, phiL_dot_L4
    real(dp) matrixfact(n_est,n_est), fact(n_est,n_est),tmp(n_est,n_est), N1(n_est,n_est), N1_L1(n_est,n_est),N1_PhiL(n_est,n_est)
    real(dp) Win12(n_est), Win34(n_est), Win43(n_est)
    real(dp) f24(n_est), f13(n_est),f31(n_est), f42(n_est)
    integer file_id, PhiLix, Lix
    integer ij(2),pq(2), est1, est2
    real(dp) this13, this24
    real(dp), allocatable :: Matrix(:,:,:,:), MatrixL1(:,:,:)
    character(LEN=10) outtag

    outtag = 'N1_All'
    allocate(matrix((lmaxout-Lmin)/Lstep+1,nPhiSample, n_est,n_est))
    allocate(matrixL1(nPhiSample,n_est,n_est))
    matrix=0
    call loadNorm()

    call WriteRanges(outtag)
    open(file=trim(dir)//trim(outtag)//trim(vartag)//'.dat', newunit = file_id, form='formatted', status='replace')

    Lix=0
    do L=LMin, lmaxout, Lstep
        Lix=Lix+1
        Lvec(1) = L
        LVec(2)= 0
        N1=0
        do L1=max(lmin_filter,dL/2), lmax, dL
            N1_L1 = 0

            matrixL1=0

            nphi=(2*L1+1)
            if (L1>3*dL) nphi=2*nint(L1/real(2*dL))+1
            dphi=(2*Pi/nphi)

            !$OMP PARALLEL DO default(shared), private(PhiIx,phi,PhiL_nphi, PhiL_phi_dphi, PhiL_phi_ix, PhiL_phi,PhiLix, dPh), &
            !$OMP private(L1vec,L2,L2vec, L2int,  L3, L3vec, L3int, L4, L4vec, L4int),&
            !$OMP private(tmp, Win12, Win34, Win43, matrixfact,fact,phiL_dot_L1, phiL_dot_L2, phiL_dot_L3, phiL_dot_L4), &
            !$OMP private(f13, f31, f42, f24, ij, pq, est1, est2,this13,this24), &
            !$OMP private(PhiL, PhiLVec, N1_PhiL), schedule(STATIC), reduction(+:N1_L1)
            do phiIx=0,(nphi-1)/2
            phi= dphi*PhiIx
            L1vec(1)=L1*cos(phi)
            L1vec(2)=L1*sin(phi)
            L2vec = Lvec-L1vec
            L2=(sqrt(L2vec(1)**2+L2vec(2)**2))
            if (L2<lmin_filter .or. L2>lmax) cycle
            L2int=nint(L2)

            call getWins(L*L1vec(1),L*L2vec(1), L1vec,real(L1,dp),L1, L2vec,L2, L2int,  Win12)

            N1_PhiL=0
            do PhiLIx = 1, nPhiSample
                PhiL = Phi_Sample(PhiLIx)
                dPh = dPhi_Sample(PhiLIx)
                PhiL_nphi=(2*PhiL+1)
                if (phiL>20) PhiL_nphi=2*nint(real(PhiL_nphi)/dPh/2)+1
                PhiL_phi_dphi=(2*Pi/PhiL_nphi)
                tmp=0
                do PhiL_phi_ix=-(PhiL_nphi-1)/2, (PhiL_nphi-1)/2
                    PhiL_phi= PhiL_phi_dphi*PhiL_phi_ix
                    PhiLvec(1)=PhiL*cos(PhiL_phi)
                    PhiLvec(2)=PhiL*sin(PhiL_phi)
                    L3vec= PhiLvec - L1vec
                    L3 = sqrt(L3vec(1)**2+L3vec(2)**2)
                    if (L3>=lmin_filter .and. L3<=lmax) then
                        L3int = nint(L3)

                        L4vec = -Lvec-L3vec
                        L4 = sqrt(L4vec(1)**2+L4vec(2)**2)
                        L4int=nint(L4)
                        if (L4>=lmin_filter .and. L4<=lmax) then
                            call getWins(-L*L3vec(1),-L*L4vec(1), L3vec,L3,L3int, L4vec,L4, L4int,  Win34, Win43)

                            phiL_dot_L1=dot_product(PhiLVec,L1vec)
                            phiL_dot_L2=-dot_product(PhiLVec,L2vec)
                            phiL_dot_L3=dot_product(PhiLVec,L3vec)
                            phiL_dot_L4=-dot_product(PhiLVec,L4vec)

                            call getResponse(phiL_dot_L1,phiL_dot_L3, L1vec,real(L1,dp),L1, L3vec,L3, L3int,  f13, f31)
                            call getResponse(phiL_dot_L2,phiL_dot_L4, L2vec,L2,L2int, L4vec,L4, L4int,  f24, f42)

                            do est1=1,n_est
                                ij=lumped_indices(:,est1)
                                do est2=est1,n_est
                                    pq=lumped_indices(:,est2)
                                    this13 = responseFor(ij(1),pq(1),f13,f31)
                                    this24 = responseFor(ij(2),pq(2),f24,f42)
                                    tmp(est1,est2)=tmp(est1,est2)+this13*this24*Win34(est2)

                                    this13 = responseFor(ij(1),pq(2),f13,f31)
                                    this24 = responseFor(ij(2),pq(1),f24,f42)
                                    tmp(est1,est2)=tmp(est1,est2)+this13*this24*Win43(est2)
                                end do
                            end do
                        end if
                    end if
                end do
                if (phiIx/=0) tmp=tmp*2 !integrate 0-Pi for phi_L1
                fact = tmp* PhiL_phi_dphi* PhiL
                do est1=1,n_est
                    matrixfact(est1,:) = fact(est1,:)*Win12(est1)*dPh
                end do

                !$OMP CRITICAL
                matrixL1(phiLix,:,:)=matrixL1(phiLix,:,:) + matrixfact
                !$OMP END CRITICAL
                N1_PhiL= N1_PhiL + fact * Cphi(PhiL)*dPh
            end do
            do est1=1,n_est
                N1_PhiL(est1,:)=N1_PhiL(est1,:)*Win12(est1)
            end do
            N1_L1 = N1_L1+N1_PhiL
        end do
        !$OMP END PARALLEL DO

        matrix(Lix,:,:,:)=matrix(Lix,:,:,:) + matrixL1*dphi*L1*dL
        N1= N1 + N1_L1 * dphi* L1*dL

    end do !L1

    do est1=1,n_est
        do est2=est1,n_est
            matrix(Lix,:,est1,est2) = matrix(Lix,:,est1,est2)*norms(L,est1)*norms(L,est2) / (twopi**4)
            N1(est1,est2) = norms(L,est1)*norms(L,est2)*N1(est1,est2) / (twopi**4)
            N1(est2,est1) = N1(est1,est2)
        end do
    end do

    write(file_id,'(1I5)',advance='NO') L
    call WriteMatrixLine(file_id, N1)

    print *, 'N1 L, TTTT, EBEB: ',L, N1(i_TT,i_TT), N1(i_eb,i_eb)

    end do
    close(file_id)

    do est1=1,n_est
      do est2=est1,n_est
          outtag = 'N1_'//estnames(est1)//estnames(est2)
         call WriteMatrix(outtag, matrix(:,:,est1,est2))
      end do
    end do

    end subroutine GetN1MatrixGeneral

    subroutine saveN0matrix(matrix,filename)
    character(LEN=*) filename
    real(dp) matrix(:,:)
    integer L, L1, file_id

    open(file=filename, newunit = file_id, form='formatted', status='replace')
    do L=LMin, lmaxout, Lstep
        write (file_id, '(1I6)', advance='NO') L
        do L1=2, lmax
            write (file_id, '(1E16.6)', advance='NO') matrix(L,L1)
        end do
        write(file_id,'(a)') ''
    end do
    close(file_id)

    end  subroutine saveN0matrix

    subroutine getEBNormMatrix(isCurl)
    integer L, l1, nphi, phiIx, L2int
    logical isCurl
    real(dp) dphi
    real(dp) L1Vec(2), L2vec(2), LVec(2)
    real(dp) phi, cos2L1L2, sin2
    real(dP) L2
    real(dp) cosfac,f12(n_est),f21(n_est),Win21(n_est),Win12(n_est), fac
    integer i,j
    real(dp) N0(n_est,n_est)
    real(dp), allocatable :: matrix2(:,:,:)
    real(dp), allocatable :: matrix(:,:,:)
    character(LEN=:), allocatable :: tag

    !Order TT, EE, EB, TE, TB, BB

    allocate(matrix2(lmaxout, lmax,n_est))
    matrix2=0
    allocate(matrix(lmaxout, lmax,n_est))
    matrix=0

    do L=LMin, lmaxout, LStep
        print *,L
        Lvec(1) = L
        LVec(2)= 0
        !$OMP PARALLEL DO default(shared), private(L1,nphi,dphi,N0,PhiIx,phi,L1vec,L2vec,L2,L2int,cos2L1L2,sin2,cosfac,f12,f21,Win12,Win21,i,j)
        do L1=2, lmax
            nphi=(2*L1+1)
            dphi=(2*Pi/nphi)
            do phiIx=0,(nphi-1)/2
                !                  do phiIx=-(nphi-1)/2, (nphi-1)/2
                phi= dphi*PhiIx
                L1vec(1)=L1*cos(phi)
                L1vec(2)=L1*sin(phi)
                L2vec = Lvec-L1vec
                L2=(sqrt(L2vec(1)**2+L2vec(2)**2))
                if (L2<2 .or. L2>lmax) cycle
                L2int=nint(L2)

                if (isCurl) then
                    call getResponse(L1vec(2)*L,L2vec(2)*L, L1vec,real(L1,dp),L1, L2vec,L2, L2int,f12,f21)
                    call getWins(L1vec(2)*L,L2vec(2)*L, L1vec,real(L1,dp),L1, L2vec,L2, L2int,  Win12, Win21)
                else
                    call getResponse(L1vec(1)*L,L2vec(1)*L, L1vec,real(L1,dp),L1, L2vec,L2, L2int,f12,f21)
                    call getWins(L1vec(1)*L,L2vec(1)*L, L1vec,real(L1,dp),L1, L2vec,L2, L2int,  Win12, Win21)
                end if


                !call getResponse(L1vec(1)*L,L2vec(1)*L, L1vec,real(L1,dp),L1, L2vec,L2, L2int,  f12, f21)
                !call getWins(L1vec(1)*L,L2vec(1)*L, L1vec,real(L1,dp),L1, L2vec,L2, L2int,  Win12, Win21)

                fac = L1*dphi/(twopi**2)
                if (PhiIx/=0) fac=fac*2
                matrix(L, L1,i_EB) = matrix(L, L1,i_EB) + f12(I_EB)*Win12(i_EB)*fac
                matrix(L, L1,i_EE) = matrix(L, L1,i_EE) + f12(I_EE)*Win12(i_EE)*fac
                matrix(L, L1,i_TT) = matrix(L, L1,i_TT) + f12(I_TT)*Win12(i_TT)*fac
                matrix(L, L1,i_TE) = matrix(L, L1,i_TE) + f12(I_TE)*Win12(i_TE)*fac

                matrix2(L, L1,i_EB) = matrix2(L, L1,i_EB) + f21(I_EB)*Win21(i_EB)*fac
                matrix2(L, L1,i_TE) = matrix2(L, L1,i_TE) + f21(I_TE)*Win21(i_TE)*fac

            end do
        end do
        !$OMP END PARALLEL DO

    end do

    tag = trim(vartag)
    if (IsCurl) tag = tag //'_Curl'

    call saveN0Matrix(matrix(:,:,i_EB),trim(dir)//'EB_NO'//tag//'_matrix1.dat')
    call saveN0Matrix(matrix2(:,:,i_EB),trim(dir)//'EB_NO'//tag//'_matrix2.dat')
    call saveN0Matrix(matrix(:,:,i_EE),trim(dir)//'EE_NO'//tag//'_matrix.dat')
    call saveN0Matrix(matrix(:,:,i_TE),trim(dir)//'TE_NO'//tag//'_matrix1.dat')
    call saveN0Matrix(matrix2(:,:,i_TE),trim(dir)//'TE_NO'//tag//'_matrix2.dat')
    call saveN0Matrix(matrix(:,:,i_TT),trim(dir)//'TT_NO'//tag//'_matrix.dat')

    end subroutine getEBNormMatrix

    subroutine getRenormMatrix()
    integer L, l1, nphi, phiIx, L2int
    real(dp) dphi
    real(dp) L1Vec(2), L2vec(2), LVec(2)
    real(dp) phi
    real(dP) L2
    real(dp) Win21(n_est),Win12(n_est), fac
    integer file_id, Lix
    real(dp), allocatable :: matrix2(:,:,:)
    real(dp), allocatable :: matrix(:,:,:)
    real(dp) f12_a,f12_b
    integer L1max

    call loadNorm()

    allocate(matrix2(lmaxout, lmax,1))
    matrix2=0
    allocate(matrix(lmaxout, lmax,1))
    matrix=0

    Lix=0
    do L=LMin, lmaxout, Lstep
        Lix=Lix+1
        print *,L
        Lvec(1) = L
        LVec(2)= 0
        !!OMP PARALLEL DO default(shared), private(f12_a,f12_b,L1,nphi,dphi,N0,PhiIx,phi,L1vec,L2vec,L2,L2int,cos2L1L2,cosfac,fac,f12,f21,Win12,Win21,i,j)
        do L1=lmin_filter, lmax
            nphi=(2*L1+1)
            dphi=(2*Pi/nphi)
            do phiIx=0,(nphi-1)/2
                phi= dphi*PhiIx
                L1vec(1)=L1*cos(phi)
                L1vec(2)=L1*sin(phi)
                L2vec = Lvec-L1vec
                L2=(sqrt(L2vec(1)**2+L2vec(2)**2))
                if (L2<lmin_filter .or. L2>lmax) cycle
                L2int=nint(L2)

                f12_a = L1vec(1)*L
                f12_b = L2vec(1)*L

                call getWins(L1vec(1)*L,L2vec(1)*L, L1vec,real(L1,dp),L1, L2vec,L2, L2int,  Win12, Win21)

                fac = L1*dphi/(twopi**2)
                if (PhiIx/=0) fac=fac*2
                !$OMP ATOMIC
                matrix(Lix, L1,i_TT) = matrix(Lix, L1,i_TT) + f12_a*Win12(i_TT)*fac
                !$OMP ATOMIC
                matrix(Lix, L2int,i_TT) = matrix(Lix, L2int,i_TT) + f12_b*Win12(i_TT)*fac
            end do
        end do
        !!OMP END PARALLEL DO
        matrix(Lix, :,i_TT) = matrix(Lix, :,i_TT)*Norms(L,i_TT)
    end do

    open(file=trim(dir)//'like/'//trim(root)//'_renorm_TT_matrix.dat', newunit = file_id, form='formatted', status='replace')
    !Lix=0
    !do L=LMin, lmaxout, Lstep
    !    Lix=Lix+1
    !    write (file_id, '(1I6)', advance='NO') L
    !    do L1=2, lmax
    !        write (file_id, '(1E16.6)', advance='NO') matrix(Lix,L1,i_TT)
    !    end do
    !    write(file_id,'(a)') ''
    !end do
    !close(file_id)
    Lix=0
    write (file_id, '(1I6)', advance='NO') 0
    do L1=lmin_filter, lmax, renorm_deltaL
        write (file_id, '(1I16)', advance='NO') L1
    end do
    write(file_id,'(a)') ''
    do L=LMin, lmaxout, Lstep
        Lix=Lix+1
        write (file_id, '(1I6)', advance='NO') L
        do L1=lmin_filter, lmax, renorm_deltaL
            L1max = min(lmax, L1+(renorm_deltaL-1)/2)
            write (file_id, '(1E16.6)', advance='NO')  sum(matrix(Lix,L1-renorm_deltaL/2:L1max,i_TT))/(L1*(L1+1)/twopi) ! * (L*(L+1.d0))**2/twopi
        end do
        write(file_id,'(a)') ''
    end do
    close(file_id)

    end subroutine getRenormMatrix

    subroutine GetN1RenormMatrix()
    integer L, l1, nphi, phiIx, PhiL_nphi, PhiL_phi_ix, L2int, L3int, L4int
    integer PhiL
    real(dp) dphi, PhiL_phi_dphi
    real(dp) L1Vec(2), L2vec(2), LVec(2), L3Vec(2),L4Vec(2), phiLVec(2)
    real(dp) N1, phi, PhiL_phi, N1_L1, N1_PhiL
    real(dP) L2, L4, L3
    real(dp) g12, g34, tmp, fact, dPh
    real(dp), allocatable :: Matrix(:,:), MatrixL1(:)
    integer file_id,  Lix, PhiLix, L1max

    call loadNorm()

    allocate(matrix((lmaxout-Lmin)/Lstep+1,lmax))

    matrix=0
    Lix=0
    do L=LMin, lmaxout, Lstep
        Lix=Lix+1
        Lvec(1) = L
        LVec(2)= 0
        N1=0
        do L1=max(lmin_filter,dL/2), lmax, dL
            N1_L1 = 0
            matrixL1=0
            nphi=(2*L1+1)
            if (L1>3*dL) nphi=2*nint(L1/real(2*dL))+1
            dphi=(2*Pi/nphi)
            !!$OMP PARALLEL DO default(shared), private(PhiIx,phi, L1vec,L2vec, L2, L2int, g12, tmp, L3, PhiL_nphi, PhiL_phi_dphi, PhiL_phi_ix, PhiL_phi, L3vec, L4vec, L4), &
            !!$OMP private(L4int,  g34, PhiL, PhiLix, PhiLVec, fact,N1_PhiL, dPh), schedule(STATIC), reduction(+:N1_L1)
            do phiIx=0,(nphi-1)/2 ! -(nphi-1)/2, (nphi-1)/2
                phi= dphi*PhiIx
                L1vec(1)=L1*cos(phi)
                L1vec(2)=L1*sin(phi)
                L2vec = Lvec-L1vec
                L2=(sqrt(L2vec(1)**2+L2vec(2)**2))
                if (L2<lmin_filter .or. L2>lmax) cycle
                L2int=nint(L2)
                g12= (L1vec(1)*L*CTf(L1) + L2vec(1)*L*CTf(L2int))/(2*CTobs(L1)*CTobs(L2int))
                N1_PhiL=0
                do PhiLIx = 1, nPhiSample
                    PhiL = Phi_Sample(PhiLIx)
                    dPh = dPhi_Sample(PhiLIx)
                    PhiL_nphi=(2*PhiL+1)
                    if (phiL>20) PhiL_nphi=2*nint(real(PhiL_nphi)/dPh/2)+1
                    PhiL_phi_dphi=(2*Pi/PhiL_nphi)
                    tmp=0
                    do PhiL_phi_ix=-(PhiL_nphi-1)/2, (PhiL_nphi-1)/2
                        PhiL_phi= PhiL_phi_dphi*PhiL_phi_ix
                        PhiLvec(1)=PhiL*cos(PhiL_phi)
                        PhiLvec(2)=PhiL*sin(PhiL_phi)

                        L3vec= PhiLvec - L1vec
                        L3 = sqrt(L3vec(1)**2+L3vec(2)**2)
                        if (L3>=lmin_filter .and. L3<=lmax) then
                            L3int = nint(L3)
                            L4vec = -Lvec-L3vec
                            L4 = sqrt(L4vec(1)**2+L4vec(2)**2)
                            if (L4>=lmin_filter .and. L4<=lmax) then
                                L4int=nint(L4)
                                g34= (L3vec(1)*CTf(L3int) + L4vec(1)*CTf(L4int))/(2*CTobs(L3int)*CTobs(L4int))
                                if (phiIx/=0) g34=g34*2 !integrate 0-Pi for phi_L1
                                g34=g34* PhiL_phi_dphi* L*PhiL*g12* Cphi(PhiL)*dPh*L1*dphi*dL/2 !/2 because have twice number of terms below
                                !$OMP ATOMIC
                                matrix(Lix,L1) = matrix(Lix,L1) -g34*dot_product(L1vec,PhiLvec)*(dot_product(L2Vec,PhiLvec)*CT(L2int)+dot_product(L4Vec,PhiLvec)*CT(L4int))
                                !$OMP ATOMIC
                                matrix(Lix,L3int) = matrix(Lix,L3int) -g34*dot_product(L3vec,PhiLvec)*(dot_product(L2Vec,PhiLvec)*CT(L2int)+dot_product(L4Vec,PhiLvec)*CT(L4int))
                                !$OMP ATOMIC
                                matrix(Lix,L2int) = matrix(Lix,L2int) -g34*(dot_product(L1vec,PhiLvec)*CT(L1)+dot_product(L3vec,PhiLvec)*CT(L3int))*dot_product(L2Vec,PhiLvec)
                                !$OMP ATOMIC
                                matrix(Lix,L4int) = matrix(Lix,L4int) -g34*(dot_product(L1vec,PhiLvec)*CT(L1)+dot_product(L3vec,PhiLvec)*CT(L3int))*dot_product(L4Vec,PhiLvec)
                                !tmp = tmp -g34*(dot_product(L1vec,PhiLvec)*CT(L1)+dot_product(L3vec,PhiLvec)*CT(L3int))&
                                !*(dot_product(L2Vec,PhiLvec)*CT(L2int)+dot_product(L4Vec,PhiLvec)*CT(L4int))
                            end if
                        end if
                    end do
                end do
            end do
            !!$OMP END PARALLEL DO
        end do
        fact = -2*Norms(L,i_TT)**2 / (twopi**4)
        matrix(Lix,:)=matrix(Lix,:)*fact
        print *, '***', L
    end do

    open(file=trim(dir)//'like/'//trim(root)//'_renorm_N1_TT_matrix.dat', newunit = file_id, form='formatted', status='replace')
    Lix=0
    write (file_id, '(1I6)', advance='NO') 0
    do L1=lmin_filter+renorm_deltaL/2, lmax, renorm_deltaL
        write (file_id, '(1I16)', advance='NO') L1
    end do
    write(file_id,'(a)') ''
    do L=LMin, lmaxout, Lstep
        Lix=Lix+1
        write (file_id, '(1I6)', advance='NO') L
        do L1=lmin_filter+renorm_deltaL/2, lmax, renorm_deltaL
            L1max = min(lmax, L1+(renorm_deltaL-1)/2)
            write (file_id, '(1E16.6)', advance='NO') sum(matrix(Lix,L1-renorm_deltaL/2:L1max)) /(L1*(L1+1)/twopi) * (L*(L+1.d0))**2/twopi
        end do
        write(file_id,'(a)') ''
    end do
    close(file_id)

    end subroutine GetN1RenormMatrix

    end module LensingBiases


    program CalcBiases
    use  LensingBiases
    implicit none
    integer i, expt
    real(dp) frac
    real(dp) deltaj


    do expt = 10,10
        select case(expt)
        case (1)
            root = 'CMB-S4'
        case (2)
            root = 'Planck'
        case (3)
            root = 'PB1s3'
        case (4)
            root = 'PB1s1'
        case (5)
            root = 'CMB-S4p'
        case (6)
            root = 'CMB-S4m'
        case (7)
            root = 'CMB-S4mm'
        case (8)
            root = 'CMB-S4mmm'
        case (9)
            root = 'CMB-S4mmmm'
        case (10)
            root = 'SOm'
        end select
        if (root=='CMB-S4') then
            noise_fwhm_deg = 3.d0/60
            muKArcmin = 1.5
            NoiseVar =  ( muKArcmin * pi/ 180 / 60.) ** 2
            NoiseVarP=NoiseVar*2
            lmin_filter = 20
            LMin = lmin_filter
            lmaxout = 3250
            lmax = 3250 !2050*2
            lmax_TT = 3250 !the L_Max for TT, even if higher for pol
            dir = 'N1/CMB-S4/'
            call system('mkdir -p ./N1/CMB-S4/')
        else if (root=='CMB-S4p') then
            noise_fwhm_deg = 3.0/60
            muKArcmin = 0.75
            NoiseVar =  ( muKArcmin * pi/ 180 / 60.) ** 2
            NoiseVarP=NoiseVar*2
            lmin_filter = 20
            LMin = lmin_filter
            lmaxout = 3250
            lmax = 6000 !2050*2
            lmax_TT = 3250 !the L_Max for TT, even if higher for pol
            dir = 'N1/CMB-S4p/'
            call system('mkdir -p ./N1/CMB-S4p/')
        else if (root=='CMB-S4m') then
            noise_fwhm_deg = 3.0/60
            muKArcmin = 3.0
            NoiseVar =  ( muKArcmin * pi/ 180 / 60.) ** 2
            NoiseVarP=NoiseVar*2
            lmin_filter = 20
            LMin = lmin_filter
            lmaxout = 3250
            lmax = 6000 !2050*2
            lmax_TT = 3250 !the L_Max for TT, even if higher for pol
            dir = 'N1/CMB-S4m/'
            call system('mkdir -p ./N1/CMB-S4m/')
        else if (root=='CMB-S4mm') then
            noise_fwhm_deg = 3.d0/60
            muKArcmin = 6.0
            NoiseVar =  ( muKArcmin * pi/ 180 / 60.) ** 2
            NoiseVarP=NoiseVar*2
            lmin_filter = 20
            LMin = lmin_filter
            lmaxout = 3250
            lmax = 6000 !2050*2
            lmax_TT = 3250 !the L_Max for TT, even if higher for pol
            dir = 'N1/CMB-S4mm/'
            call system('mkdir -p ./N1/CMB-S4mm/')
        else if (root=='CMB-S4mmm') then
            noise_fwhm_deg = 3.0/60
            muKArcmin = 12.0
            NoiseVar =  ( muKArcmin * pi/ 180 / 60.) ** 2
            NoiseVarP=NoiseVar*2
            lmin_filter = 20
            LMin = lmin_filter
            lmaxout = 3250
            lmax = 6000 !2050*2
            lmax_TT = 3250 !the L_Max for TT, even if higher for pol
            dir = 'N1/CMB-S4mmm/'
            call system('mkdir -p ./N1/CMB-S4mmm/')
        else if (root=='CMB-S4mmmm') then
            noise_fwhm_deg = 9.0/60
            muKArcmin = 1.5
            NoiseVar =  ( muKArcmin * pi/ 180 / 60.) ** 2
            NoiseVarP=NoiseVar*2
            lmin_filter = 20
            LMin = lmin_filter
            lmaxout = 3050
            lmax = 6000 !2050*2
            lmax_TT = 6000 !the L_Max for TT, even if higher for pol
            dir = 'N1/CMB-S4mmmm/'
            call system('mkdir -p ./N1/CMB-S4mmmm/')
        else if (root=='Planck') then
            noise_fwhm_deg = 7.d0/60
            muKArcmin = 27
            NoiseVar =  ( muKArcmin * pi/ 180 / 60.) ** 2
            NoiseVarP=NoiseVar*2
            lmin_filter = 2
            LMin = lmin_filter
            lmaxout = 2150
            lmax = 4000 !2050*2
            lmax_TT = 4000 !the L_Max for TT, even if higher for pol
            dir = 'N1/Planck/'
            call system('mkdir -p ./N1/Planck/')
        else if (root=='PB1s3') then
            noise_fwhm_deg = 3.5d0/60
            muKArcmin = 35.5/2.0
            NoiseVar =  ( muKArcmin * pi/ 180 / 60.) ** 2
            NoiseVarP=NoiseVar*2
            lmin_filter = 30
            LMin = lmin_filter
            lmaxout = 3150
            lmax = 6000 !2050*2
            lmax_TT = 6000 !the L_Max for TT, even if higher for pol
            dir = 'N1/PB1s3/'
            call system('mkdir -p ./N1/PB1s3/')
        else if (root=='PB1s1') then
            noise_fwhm_deg = 3.5d0/60
            muKArcmin = 7.0
            NoiseVar =  ( muKArcmin * pi/ 180 / 60.) ** 2
            NoiseVarP=NoiseVar*2
            lmin_filter = 150
            LMin = lmin_filter
            lmaxout = 2150
            lmax = 3000 !2050*2
            lmax_TT = 3000 !the L_Max for TT, even if higher for pol
            dir = 'N1/PB1s1/'
            call system('mkdir -p ./N1/PB1s1/')
        else if (root=='SOm') then
            noise_fwhm_deg = 1.0/60
            muKArcmin = 0.5
            NoiseVar =  ( muKArcmin * pi/ 180 / 60.) ** 2
            NoiseVarP=NoiseVar*2
            lmin_filter = 30
            LMin = lmin_filter
            lmaxout = 5000
            lmax = 5000 !2050*2
            lmax_TT = 3000 !the L_Max for TT, even if higher for pol
            dir = 'N1/SOm/'
            call system('mkdir -p ./N1/SOm/')
        end if

        if (.false.) then
            !planck
            root = 'g60_full_pttptt'

            !        call ReadPhiPhi(dir//'like/'//root//'_fid_phi.dat')
            call ReadPower(dir//'like/'//root//'_fid_cl.dat')
            call NoiseInit(NoiseVar,NoiseVarP) !for pol
            CEobs = CE + NP
            CBobs = CB + NP
            call ReadFilter(dir//'like/'//root//'_TT_Filter.dat')
        else if (.true.) then
            !planck
            !call ReadPhiPhi('HighLExtrapTemplate_lenspotentialCls.dat')
            !call ReadPower('HighL_lensedCls.dat')
            !C:\Tmp\Planck\FFP9\
            call ReadPhiPhi('additional_files/planck_lensing_wp_highL_bestFit_20130627_lenspotentialCls.dat')
            call ReadPower('additional_files/planck_lensing_wp_highL_bestFit_20130627_lensedCls.dat')

            call NoiseInit(NoiseVar, NoiseVarP)
            CTobs = CT + NT
            CEobs = CE + NP
            CBobs = CB + NP
        end if
        vartag = '_'//root

        if (.false.) then
            !Get N0 for various noise levels, e.g. for estimating weighted S/N
            deltaj = 0.03
            do i=-15,15
                frac = 1 + i*deltaj
                call NoiseInit(frac*NoiseVar, frac*NoiseVarP)
                CTobs = CT + NT
                CEobs = CE + NP
                CBobs = CB + NP
                write(vartag,*) i
                vartag = trim(adjustl(vartag))
                call getNorm(.true., .true.)

            end do
            stop
        end if


        call getNorm(.false., .false.)

        call SetPhiSampling(lmax, UsePhiSampling)

        !call GetN1Temp()
        call GetN1MatrixGeneral()

        !  call GetN1General()
        !call getRenormMatrix()
        !       call getN1RenormMatrix()

        !call getEBNormMatrix(.false.)
        !call getEBNormMatrix(.true.)
    end do

    end program CalcBiases
