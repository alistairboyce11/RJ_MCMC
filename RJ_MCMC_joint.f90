!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! This progranm does a Transdimensional inversion of Receiver functions
! with the reversible jump algorithm

! Thomas Bodin, ANU, December 2011
!*********************************************************************
!/////////////////////////////////////////////////////////////////////
!********************************************************************

program RJ_MCMC
    use mpi
    implicit none
    include 'params.h'
    include 'data_joint.h'

    !***************************************************************

    ! DECLARATIONS OF VARIABLES

    !****************************************************************

    real , EXTERNAL    ::    gasdev,ran3,interp
    real log, sqrt

    integer i,ii,sample,ind,ount,k ! various indices
    integer th,th_all
    integer nlims_cur,nlims_cur_diff
    logical ier, error_flag ! error flags for dispersion curves and minos_bran
    integer npt,npt_prop,npt_iso,npt_ani ! numbers of points, isotropic, anisotropic
    logical accept,tes,birth,birtha,death,deatha,move,value,noisd_R,noisd_L,ani,change_vp !which parameter did we change?
    integer ind2,j !more indices, mostly for posterior. j sometimes is the number of processors, careful at the end!
    real peri_R(ndatadmax),peri_L(ndatadmax) !periods of data
    real peri_R_tmp(ndatadmax),peri_L_tmp(ndatadmax) !periods of data
    integer n_R(ndatadmax),n_L(ndatadmax),nmax_R,nmax_L,nmin_R,nmin_L !harmonic number of modes
    integer n_R_tmp(ndatadmax),n_L_tmp(ndatadmax)
    integer nlims_R(2,nharmo_max),nlims_L(2,nharmo_max),numharm_count,numharm_R(nharmo_max),numharm_L(nharmo_max)

    real PrB,AcB,PrD,AcD,Prnd_R,Acnd_R,Prnd_L,Acnd_L,&
        Acba,Prba,Acda,Prda,out,Prxi,Acxi,Pr_vp,&
        Ac_vp,PrP,PrV,AcP,AcV ! to adjust acceptance rates for all different parameters
    real lsd_L,lsd_L_prop,lsd_R,lsd_R_prop !logprob without uncertainties
    real logprob_vsv,logprob_vp !for reversible jump
    real voro(malay,4)
    real voro_prop(malay,4)
    real t1,t2 !timers
    real like,like_prop,like_w,like_prop_w,u,& !log-likelyhoods
        liked_R_prop,liked_R,liked_L_prop,liked_L
    double precision logrsig,Ad_R,Ad_R_prop,Ad_L,Ad_L_prop !uncertainty parameters
    real sigmav,sigmav_old,sigmav_new,AR_birth_old   !proposal on velocity when Birth move
    real d_obsdcR(ndatadmax),d_obsdcL(ndatadmax),d_obsdcRe(ndatadmax),d_obsdcLe(ndatadmax) !observed data
    integer inorout(21000),inorouts(21000),members(21000),io ! mpi related variables
    integer ndatad_R,ndatad_L,ndatad_R_tmp,ndatad_L_tmp !number of observed data points
    real pxi,p_vp,pd,pv,pAd_R,pAd_L ! related to observed data for azimuthal anisotropy
    integer nptref,malayv ! number of points in reference model, number of layers in voro

    logical isoflag(malay),isoflag_prop(malay) ! is this layer isotropic?
    real d_cR(ndatadmax),d_cL(ndatadmax) ! phase velocity as simulated by forward modelling
    real d_cR_tmp(ndatadmax),d_cL_tmp(ndatadmax)
    real rq_R(ndatadmax),rq_L(ndatadmax) ! errors from modelling
    real d_cR_prop(ndatadmax),d_cL_prop(ndatadmax) ! phase velocity as simulated by forward modelling
    real d_cR_prop_tmp(ndatadmax),d_cL_prop_tmp(ndatadmax) ! phase velocity as simulated by forward mod
    real rq_R_prop(ndatadmax),rq_L_prop(ndatadmax)

    !for MPI
    integer ra,ran,rank, nbproc, ierror ,status(MPI_STATUS_SIZE),group_world,good_group,MPI_COMM_small,flag
    ! ierror: MPI-related error status, nbproc: MPI-related, number of processes, rank, group_world: MPI-related
    character filename*130, number*4, filenametmp*300
    character filenamemax*30  !filename for writing the best model
    real model_ref(mk,9) ! reference model, all models are deviations from it
    real,dimension(mk) :: r,rho,vpv,vph,vsv,vsh,qkappa,qshear,eta,xi,vp_data,vs_data ! input for forward modelling
    integer nptfinal,nic,noc,nic_ref,noc_ref,jcom ! more inputs
    integer nmodes,n_mode(nmodes_max),l_mode(nmodes_max) ! outputs of forward modelling
    real c_ph(nmodes_max),period(nmodes_max),raylquo(nmodes_max),tref ! more outputs, rayquo: error measurement, should be of the order of eps
    real wmin_R(nharmo_max),wmax_R(nharmo_max),wmin_L(nharmo_max),wmax_L(nharmo_max) ! more inputs f

    integer nharm_R,nharm_L,iharm
    integer numdis ! additional stuff for joint inversion

    integer i_w
    real widening
    logical burnin_in_progress

    ! covergence plots
    real ncell(burn_in_widening+nsample_widening)
    real ncells(burn_in_widening+nsample_widening)
    real convBs(burn_in_widening+nsample_widening)
    real convB(burn_in_widening+nsample_widening)
    real convPs(burn_in_widening+nsample_widening)
    real convP(burn_in_widening+nsample_widening) ! birth rates, vp change rates
    real convvs(burn_in_widening+nsample_widening)
    real convvss(burn_in_widening+nsample_widening)
    real convdp(burn_in_widening+nsample_widening)
    real convdps(burn_in_widening+nsample_widening)
    real convvp(burn_in_widening+nsample_widening)
    real convvps(burn_in_widening+nsample_widening)
    real convxis(burn_in_widening+nsample_widening)
    real convxi(burn_in_widening+nsample_widening)
    real convBas(burn_in_widening+nsample_widening)
    real convBa(burn_in_widening+nsample_widening) ! anisotropic birth rates
    real convDs(burn_in_widening+nsample_widening)
    real convD(burn_in_widening+nsample_widening)
    real convDas(burn_in_widening+nsample_widening)
    real convDa(burn_in_widening+nsample_widening)
    real convd_R(burn_in_widening+nsample_widening)
    real convd_L(burn_in_widening+nsample_widening)
    real convd_Rs(burn_in_widening+nsample_widening)
    real convd_Ls(burn_in_widening+nsample_widening) ! the proposed model, history of lsd
    real convAd_R(burn_in_widening+nsample_widening)
    real convAd_Rs(burn_in_widening+nsample_widening)
    real convAd_L(burn_in_widening+nsample_widening)
    real convAd_Ls(burn_in_widening+nsample_widening) !variations of uncertainty parameters
    real convsigmav(burn_in_widening+nsample_widening)
    real convsigmavs(burn_in_widening+nsample_widening)
    real convsigmavp(burn_in_widening+nsample_widening)
    real convsigmavps(burn_in_widening+nsample_widening)

    real ncell_bi(burn_in)
    real ncells_bi(burn_in)
    real convBs_bi(burn_in)
    real convB_bi(burn_in)
    real convPs_bi(burn_in)
    real convP_bi(burn_in) ! birth rates, vp change rates
    real convvs_bi(burn_in)
    real convvss_bi(burn_in)
    real convdp_bi(burn_in)
    real convdps_bi(burn_in)
    real convvp_bi(burn_in)
    real convvps_bi(burn_in)
    real convxis_bi(burn_in)
    real convxi_bi(burn_in)
    real convBas_bi(burn_in)
    real convBa_bi(burn_in) ! anisotropic birth rates
    real convDs_bi(burn_in)
    real convD_bi(burn_in)
    real convDas_bi(burn_in)
    real convDa_bi(burn_in)
    real convd_R_bi(burn_in)
    real convd_L_bi(burn_in)
    real convd_Rs_bi(burn_in)
    real convd_Ls_bi(burn_in) ! the proposed model, history of lsd
    real convAd_R_bi(burn_in)
    real convAd_Rs_bi(burn_in)
    real convAd_L_bi(burn_in)
    real convAd_Ls_bi(burn_in) !variations of uncertainty parameters
    real convsigmav_bi(burn_in)
    real convsigmavs_bi(burn_in)
    real convsigmavp_bi(burn_in)
    real convsigmavps_bi(burn_in)

    ! for mpi write purposes
    integer filehandle
    integer, dimension(MPI_STATUS_SIZE) :: status_write_mpi
    integer nb_bytes_integer,nb_bytes_real,nb_bytes_double
    integer nb_bytes_header,nb_bytes_model
    integer(kind=MPI_OFFSET_KIND) :: position_file
    integer :: num_file=-1

1000 format(I4)


    !***********************************************************************

    CALL cpu_time(t1)  !Tic. start counting time

    call MPI_INIT(ierror)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nbproc, ierror) ! MPI_COMM_WORLD: communicator (https://www.open-mpi.org/doc/v4.0/man3/MPI_Comm_size.3.php)
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror) ! https://www.open-mpi.org/doc/v4.0/man3/MPI_Comm_rank.3.php
    call MPI_Comm_group(MPI_COMM_WORLD,group_world,ierror) ! https://www.open-mpi.org/doc/v3.0/man3/MPI_Comm_group.3.php

    ! Get byte size of different types
    call MPI_TYPE_SIZE(MPI_INTEGER,nb_bytes_integer,ierror)
    call MPI_TYPE_SIZE(MPI_REAL,nb_bytes_real,ierror)
    call MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION,nb_bytes_double,ierror)
    nb_bytes_header=8*nb_bytes_integer+8*nb_bytes_real+4*nb_bytes_double ! number of bytes of the header
    write(*,*)'nb bytes header',nb_bytes_header
    nb_bytes_model=7*nb_bytes_integer+(1+4*mk+2*ndatadmax)*nb_bytes_real+2*nb_bytes_double ! number of bytes in each model
    write(*,*)'nb_bytes_model',nb_bytes_model

    !Start Parralelization of the code. From now on, the code is run on each
    !processor independently, with ra = the number of the proc.

    ra=0

    !************************************************************
    !                READ PREM
    !************************************************************

    nlims_cur=1
    open(7,file="Modified_PREM_GLOBAL.in",status='old',form='formatted')
    !open(7,file="Model_PREM_SIMPLE.in",status='old',form='formatted')
    !  110 format(20a4)
    read(7,*) nptref,nic_ref,noc_ref
    do i = 1,nptref
         read(7,*) model_ref(i,1),&
                    model_ref(i,2),model_ref(i,3),model_ref(i,4),model_ref(i,5),&
                    model_ref(i,6),model_ref(i,7),model_ref(i,8),model_ref(i,9)
    end do
    close(7)
    j=1
    do while (model_ref(j,1)<rearth-d_max*1000.)
        j=j+1
    end do
    j=j-1

    write(*,*)'read prem'
    open(65,file=dirname//'/dispersion_all.in',status='old')! 65: name of the opened file in memory (unit identifier)
    read(65,*,IOSTAT=io)numdis
    read(65,*,IOSTAT=io) ! I don't need to read the individual dispersion curves anymore
    read(65,*,IOSTAT=io)
    read(65,*,IOSTAT=io)
    read(65,*,IOSTAT=io)ndatad_R
    read(65,*,IOSTAT=io)nharm_R
    do k=1,nharm_R
        read(65,*,IOSTAT=io)numharm_R(k)
        read(65,*,IOSTAT=io)nlims_cur_diff
        nlims_R(1,k)=nlims_cur
        nlims_R(2,k)=nlims_R(1,k)+nlims_cur_diff
        do i=nlims_cur,nlims_cur+nlims_cur_diff
            read(65,*,IOSTAT=io)n_R(i),peri_R(i),d_obsdcR(i),d_obsdCRe(i)
            do j=1,numdis
                read(65,*,IOSTAT=io) ! I don't need to read the individual dispersion curves anymore
            enddo
        enddo
        nlims_cur=nlims_R(2,k)+1
    enddo

    nlims_cur=1
    read(65,*,IOSTAT=io)ndatad_L
    read(65,*,IOSTAT=io)nharm_L
    do k=1,nharm_L
        read(65,*,IOSTAT=io)numharm_L(k)
        read(65,*,IOSTAT=io)nlims_cur_diff
        nlims_L(1,k)=nlims_cur
        nlims_L(2,k)=nlims_L(1,k)+nlims_cur_diff
        do i=nlims_cur,nlims_cur+nlims_cur_diff
            read(65,*,IOSTAT=io)n_L(i),peri_L(i),d_obsdcL(i),d_obsdCLe(i)
            do j=1,numdis
                read(65,*,IOSTAT=io) ! I don't need to read the individual dispersion curves anymore
            enddo
        enddo
        nlims_cur=nlims_L(2,k)+1
    enddo
    close(65)

    do k=1,nharm_R
        wmin_R(k)=1000./(maxval(peri_R(nlims_R(1,k):nlims_R(2,k)))+20)
        wmax_R(k)=1000./(minval(peri_R(nlims_R(1,k):nlims_R(2,k)))-2)
    enddo

    do k=1,nharm_L
        wmin_L(k)=1000./(maxval(peri_L(nlims_L(1,k):nlims_L(2,k)))+20)
        wmax_L(k)=1000./(minval(peri_L(nlims_L(1,k):nlims_L(2,k)))-2)
    enddo

    nmin_R=minval(n_R(:ndatad_R))
    nmax_R=maxval(n_R(:ndatad_R))

    nmin_L=minval(n_L(:ndatad_L))
    nmax_L=maxval(n_L(:ndatad_L))

    if (ndatad_R>=ndatad_L) tref=sum(peri_R(:ndatad_R))/ndatad_R
    if (ndatad_L>ndatad_R) tref=sum(peri_L(:ndatad_L))/ndatad_L

    write(*,*)'DONE INITIALIZING'

    !**************************************************************

    !    Draw the first model randomly from the prior distribution

    !**************************************************************

    write(*,*)dirname

    pxi = 0.4             ! proposal for change in xi
    p_vp = 0.1           ! proposal for change in vp/vsv
    pd = 10!0.2         ! proposal on change in position
    pv = 0.1!0.04     ! proposal on velocity
    pAd_R = 0.5        ! proposal for change in R noise
    pAd_L = 0.5        ! proposal for change in L noise
    sigmav=0.15         ! proposal for vsv when creating a new layer



    ra=rank !seed for RNG
    ran=rank
    inorout=0
    inorouts=0

    ier=.false.

    ! Initilalise sigma (noise parameter)
    Ad_R = Ad_R_max-pAd_R! Ad_min+ran3(ra)*(Ad_R_max-Ad_min)
    Ad_L = Ad_L_max-pAd_L

    !! start from true Ad
    !Ad_R = 10
    !Ad_L = 5

    j=0
    tes=.false.
    do while(.not.tes)
150      tes=.true.
        ! we want to use the real model as a starting model
        ! create a starting model randomly
        !npt = milay+ran3(ra)*(malay-milay)!(12-milay)!(maxlay-minlay)
        npt=15 ! not too many or too few initial layers to avoid getting stuck right away
        if ((npt>malay).or.(npt<milay)) goto 150
        j=j+1 ! première couche à 0 km/1km qui ne change jamais ???

        ! first layer always at 0
        voro(1,1)= 0.
        voro(1,2)= (-0.1+2*0.1*ran3(ra)) ! have an initial model reasonably close to the reference
          if (ran3(ra)<0.5) then
              voro(1,3)= 0.8+(1.2-0.8)*ran3(ra)
          else
              voro(1,3)= -1
          endif
        !voro(i,3)=-1
        voro(1,4)= vp_min+(vp_max-vp_min)*ran3(ra)


        do i=2,npt

            voro(i,1)= d_min+ran3(ra)*(d_max-d_min)
            voro(i,2)= (-0.1+2*0.1*ran3(ra)) ! have an initial model reasonably close to the reference
              if (ran3(ra)<0.5) then
                  voro(i,3)= 0.8+(1.2-0.8)*ran3(ra)
              else
                  voro(i,3)= -1
              endif
            !voro(i,3)=-1
            voro(i,4)= vp_min+(vp_max-vp_min)*ran3(ra)
        enddo

        !write(*,*)'before combine_linear'
        call combine_linear_vp(model_ref,nptref,nic_ref,noc_ref,voro,npt,d_max,&
            r,rho,vpv,vph,vsv,vsh,qkappa,qshear,eta,nptfinal,nic,noc,vs_data,xi,vp_data)

        !write(*,*)'after combine_linear'

        if (ndatad_R>0) then

            jcom=3 !rayleigh waves

            do iharm=1,nharm_R
                nmodes=0
                call minos_bran(1,tref,nptfinal,nic,noc,r,rho,vpv,vph,vsv,vsh,&
                    qkappa,qshear,eta,eps,wgrav,jcom,lmin,lmax,&
                    wmin_R(iharm),wmax_R(iharm),numharm_R(iharm),numharm_R(iharm),&
                    nmodes_max,nmodes,n_mode,l_mode,c_ph,period,raylquo,error_flag) ! calculate normal modes
                if (error_flag) then
                    tes=.false.
                    write(*,*)"Minos_bran FAILED for RAYLEIGH 004"
                end if
                peri_R_tmp(:nlims_R(2,iharm)-nlims_R(1,iharm)+1)=peri_R(nlims_R(1,iharm):nlims_R(2,iharm))
                n_R_tmp(:nlims_R(2,iharm)-nlims_R(1,iharm)+1)=n_R(nlims_R(1,iharm):nlims_R(2,iharm))
                ndatad_R_tmp=nlims_R(2,iharm)-nlims_R(1,iharm)+1 ! fortran slices take the first and the last element
                call dispersion_minos(nmodes_max,nmodes,n_mode,c_ph,period,raylquo,&
                    peri_R_tmp,n_R_tmp,d_cR_tmp,rq_R,ndatad_R_tmp,ier) ! extract phase velocities from minos output (pretty ugly)
                d_cR(nlims_R(1,iharm):nlims_R(2,iharm))=d_cR_tmp(:nlims_R(2,iharm)-nlims_R(1,iharm)+1)
                if (ier) then
                    tes=.false.
                    write(*,*)"Minos_bran FAILED for RAYLEIGH 005"
                endif
                if (maxval(abs(rq_R(:ndatad_R_tmp)))>maxrq*eps) tes=.false.
            enddo
        endif

        if (ndatad_L>0) then

            jcom=2 !love waves
            do iharm=1,nharm_L
                nmodes=0
                call minos_bran(1,tref,nptfinal,nic,noc,r,rho,vpv,vph,vsv,vsh,&
                    qkappa,qshear,eta,eps,wgrav,jcom,lmin,lmax,&
                    wmin_L(iharm),wmax_L(iharm),numharm_L(iharm),numharm_L(iharm),&
                    nmodes_max,nmodes,n_mode,l_mode,c_ph,period,raylquo,error_flag) ! calculate normal modes
                if (error_flag) then
                    tes=.false.
                    write(*,*)"Minos_bran FAILED for LOVE 004"
                end if
                peri_L_tmp(:nlims_L(2,iharm)-nlims_L(1,iharm)+1)=peri_L(nlims_L(1,iharm):nlims_L(2,iharm))
                n_L_tmp(:nlims_L(2,iharm)-nlims_L(1,iharm)+1)=n_L(nlims_L(1,iharm):nlims_L(2,iharm))
                ndatad_L_tmp=nlims_L(2,iharm)-nlims_L(1,iharm)+1 ! fortran slices take the first and the last element
                call dispersion_minos(nmodes_max,nmodes,n_mode,c_ph,period,raylquo,&
                peri_L_tmp,n_L_tmp,d_cL_tmp,rq_L,ndatad_L_tmp,ier) ! extract phase velocities from minos output (pretty ugly)
                d_cL(nlims_L(1,iharm):nlims_L(2,iharm))=d_cL_tmp(:nlims_R(2,iharm)-nlims_R(1,iharm)+1)
                if (ier) then
                    tes=.false.
                    write(*,*)"Minos_bran FAILED for LOVE 005"
                endif
                if (maxval(abs(rq_L(:ndatad_L_tmp)))>maxrq*eps) tes=.false.
            enddo

        endif

        if (j>250) stop "CAN'T INITIALIZE MODEL" ! if it doesn't work after 250 tries, give up

    enddo

    write(*,*)'made initial model'


    ! isoflag says if a layer is isotropic
    do i=1,npt
        if (voro(i,3)==-1) then
            isoflag(i)=.true.
        else
            isoflag(i)=.false.
        endif
    enddo

    !***********************************************************

    !                 Get initial likelihood

    !***********************************************************

    lsd_R=0
    lsd_L=0
    liked_R=0
    liked_L=0
    do i=1,ndatad_R
        lsd_R=lsd_R+(d_obsdCR(i)-d_cR(i))**2
        liked_R=liked_R+(d_obsdCR(i)-d_cR(i))**2/(2*(Ad_R*d_obsdCRe(i)*(d_obsdCR(i)/100))**2) ! gaussian errors
    enddo
    do i=1,ndatad_L
        lsd_L=lsd_L+(d_obsdCL(i)-d_cL(i))**2
        liked_L=liked_L+(d_obsdCL(i)-d_cL(i))**2/(2*(Ad_L*d_obsdCLe(i)*(d_obsdCL(i)/100))**2)
    enddo

    like= (liked_R + liked_L)
    like_w=like/widening

    if (ran==0) write(*,*)widening
    if (ran==0) write(*,*)like,like_w

    lsd_R=0
    lsd_L=0

    ncell=0
    ncells=0
    convd_R=0
    convd_L=0
    convd_Rs=0
    convd_Ls=0
    convB=0
    convBs=0
    convBa=0
    convBas=0
    convD=0
    convDs=0
    convvs=0
    convvss=0
    convdp=0
    convdps=0
    convDa=0
    convDas=0
    convP=0
    convPs=0
    ier=.false.

    write(*,*)'DONE INITIALIZING'

    write(*,*)'all initialized'

    burnin_in_progress=.true.

    widening=widening_start
    do i_w=1,n_w

        num_file=-1

        inorout=0
        inorouts=0

        ncell=0
        ncells=0
        convd_R=0
        convd_L=0
        convd_Rs=0
        convd_Ls=0
        convB=0
        convBs=0
        convBa=0
        convBas=0
        convD=0
        convDs=0
        convvs=0
        convvss=0
        convdp=0
        convdps=0
        convDa=0
        convDas=0
        convP=0
        convPs=0
        convAd_R=0
        convAd_Rs=0
        convAd_L=0
        convAd_Ls=0
        convsigmav=0
        convsigmavs=0

        ier=.false.

        sample=0
        th=0
        ount=0
        PrP=0
        PrV=0
        PrB=0
        PrD=0
        PrBa=0
        PrDa=0
        AcP=0
        AcV=0
        AcB=0
        AcD=0
        AcBa=0
        AcDa=0
        Acnd_R=0
        Prnd_R=0
        Prxi=0
        Acxi=0
        Pr_vp=0
        Ac_vp=0
        Prnd_L=0
        Acnd_L=0

        sigmav_old=0
        Ar_birth_old=0

        write(*,*)voro

        do while ((sample<nsample_widening).or.burnin_in_progress) ! main loop, sample: number of sample post burn-in
            ount=ount+1
            malayv=malay
            if (mod(ount,every)==0) then ! check regularly if acceptance rates are in an acceptable range, else change proposal width

                if ((Ac_vp/(Pr_vp+1))>0.54) p_vp=p_vp*(1+perturb)! if not increase width of proposition density
                if ((Ac_vp/(Pr_vp+1))<0.34) p_vp=p_vp*(1-perturb)! or decrease it
                if ((Acxi/(Prxi+1))>0.54) pxi=pxi*(1+perturb)
                if ((Acxi/(Prxi+1))<0.34) pxi=pxi*(1-perturb)
                if ((Acnd_R/(Prnd_R+1))>0.54) pAd_R=pAd_R*(1+perturb) ! for rayleigh waves
                if ((Acnd_R/(Prnd_R+1))<0.34) pAd_R=pAd_R*(1-perturb)
                if ((Acnd_L/(Prnd_L+1))>0.54) pAd_L=pAd_L*(1+perturb) ! for love waves
                if ((Acnd_L/(Prnd_L+1))<0.34) pAd_L=pAd_L*(1-perturb)

                if ((Acv/(Prv+1))>0.54) pv=pv*(1+perturb) ! 2 layers for vsv
                if ((Acv/(Prv+1))<0.34) pv=pv*(1-perturb)

                if ((Acp/(Prp+1))>0.54) pd=pd*(1+perturb) ! 2 layers for changing depth
                if ((Acp/(Prp+1))<0.34) pd=pd*(1-perturb)

                !------------------------------------------------
                ! UPDATE SIGMAV

                !write(*,*) '****',AcB/PrB,(Ar_birth_old/100.),0.5*((AcB/PrB)-(AcD/PrD))
                if ((abs((AcB/PrB)-Ar_birth_old)>2*abs((AcB/PrB)-(AcD/PrD))).and.&
                    (AcB.ne.0)) then !special treatement for adding/removing layers
                    !write(*,*)
                    !write(*,*)'update ---------------------'
                    !write(*,*)Ar_birth_old,100*AcB/PrB
                    !write(*,*)sigmav_old,sigmav

                    if ((AcB/PrB)>Ar_birth_old) then! Going in the right direction
                        if (sigmav>sigmav_old) then !Going up
                            sigmav_new=sigmav*(1+perturb)
                        else !Going down
                            sigmav_new=sigmav*(1-perturb)
                        endif
                    else ! Going in the wrong direction
                        if (sigmav>sigmav_old) then !Going up
                            sigmav_new=sigmav*(1-perturb)
                        else
                            sigmav_new=sigmav*(1+perturb)
                        endif
                    endif
                    !write(*,*)sigmav_new
                    !write(*,*)
                    sigmav_old=sigmav
                    sigmav=sigmav_new
                    Ar_birth_old=AcB/PrB
                    PrB=0
                    PrD=0
                    AcB=0
                    AcD=0
                endif

                PrP=0
                PrV=0
                PrBa=0
                PrDa=0
                AcP=0
                AcV=0
                AcBa=0
                AcDa=0
                Acnd_R=0
                Prnd_R=0
                Acxi=0
                Prxi=0
                Ac_vp=0
                Pr_vp=0
                Prnd_L=0
                Acnd_L=0
            endif

            isoflag_prop = isoflag
            voro_prop=voro
            like_prop =like
            liked_R_prop=liked_R
            liked_L_prop=liked_L


            npt_prop = npt

            lsd_R_prop = lsd_R
            lsd_L_prop =lsd_L
            Ad_R_prop = Ad_R
            Ad_L_prop = Ad_L


            u=ran3(ra)
            ind=1
            out=1 ! indicates if change is acceptable, i. e. not out of bounds etc.
            move=.false.
            value=.false.
            birth=.false.
            birtha=.false.
            deatha=.false.
            death=.false.
            noisd_R=.false.
            noisd_L=.false.
            logprob_vsv=0
            !logprob_vp=0
            ani=.false.
            change_vp=.false.


            npt_ani=0
            npt_iso=0

            do i=1,npt
                if (isoflag(i).eqv..true.) npt_iso=npt_iso+1
                if (isoflag(i).eqv..false.) npt_ani=npt_ani+1
            enddo
            if ((npt_iso+npt_ani).ne.npt) stop "Error here"

            !*************************************************************

            !                   Propose a new model

            !*************************************************************
            if (u<0) then
                continue

            elseif (u<0.1) then !change xi--------------------------------------------
                ani=.true.
                Prxi=Prxi+1 !increase counter to calculate acceptance rates

                if (npt_ani.ne.0) then

                    ! Choose randomly an anisotropic cell among the npt_ani possibilities
                    ind2 = ceiling(ran3(ra)*(npt_ani))
                    if (ind2==0) ind2=1 !number of the anisotropic cell

                    j=0
                    do ii=1,npt
                        if (isoflag(ii).eqv..false.) j=j+1
                        if (j==ind2) then
                            ind=ii ! number of the cell in voro
                            exit
                        endif
                    enddo

                    if (ind>npt) stop "684"
                    if (voro(ind,3)==-1) stop "874"
                    voro_prop(ind,3)=voro(ind,3)+gasdev(ra)*pxi
                    if (isoflag(ind).eqv..true.) stop "isoflag1"

                    !Check if oustide bounds of prior
                    if ((voro_prop(ind,3)<=xi_min).or.(voro_prop(ind,3)>=xi_max)) then
                        out=0
                    endif
                else
                    out=0
                endif
            elseif (u<0.2) then !change vp --------------------------------------------
                change_vp=.true.

                ! Choose a random cell
                ind=ceiling(ran3(ra)*npt)
                if (ind==0) ind=1

                Pr_vp=Pr_vp+1
                if (ind>npt) then
                    out=0
                else
                    voro_prop(ind,4)=voro(ind,4)+gasdev(ra)*p_vp

                    !Check if oustide bounds of prior
                    if ((voro_prop(ind,4)<=vp_min).or.(voro_prop(ind,4)>=vp_max)) out=0

                endif
            elseif (u<0.3) then !change position--------------------------------------------
                move=.true.
                ! layer 1 always stays at depth 0
                ind=1+ceiling(ran3(ra)*(npt-1))
                if (ind==1) ind=2

                PrP=PrP+1
                voro_prop(ind,1)=voro(ind,1)+gasdev(ra)*pd

                if ((voro_prop(ind,1)<=d_min).or.(voro_prop(ind,1)>=d_max)) then
                    out=0
                endif
                if ((voro_prop(ind,2)<=-width).or.(voro_prop(ind,2)>=width)) then
                    write(*,*)'bad vx when changing depth'
                    out=0
                endif

            elseif (u<0.4) then ! Change noise parameter for rayleigh waves
                noisd_R=.true.
                Prnd_R = Prnd_R + 1
                Ad_R_prop = Ad_R+gasdev(ra)*pAd_R
                !Check if oustide bounds of prior
                if ((Ad_R_prop<=Ad_R_min).or.(Ad_R_prop>=Ad_R_max)) then
                    out=0
                endif
            elseif (u<0.5) then ! Change noise parameter for love waves
                noisd_L=.true.
                Prnd_L = Prnd_L + 1
                Ad_L_prop = Ad_L+gasdev(ra)*pAd_L
                !Check if oustide bounds of prior
                if ((Ad_L_prop<=Ad_L_min).or.(Ad_L_prop>=Ad_L_max)) then
                    out=0
                endif

            elseif (u<.6) then ! change vsv-----------------------------------

                value=.true.
                ind=ceiling(ran3(ra)*npt)
                if (ind==0) ind=1

                !if (ount.GT.burn_in) then
                !write(*,*)voro(ind,1),d_max/2
                PrV=PrV+1
                voro_prop(ind,2)=voro(ind,2)+gasdev(ra)*pv

                if (ind>npt) then
                    write(*,*)npt,ind
                    stop "763"
                endif

                !Check if oustide bounds of prior
                if ((voro_prop(ind,2)<=-width).or.(voro_prop(ind,2)>=width)) then
                    out=0
                endif
                !if (out==1) write(*,*)ind


            elseif (u<0.7) then !Birth of an isotropic cell -------------------------------------
                birth = .true.
                PrB = PrB + 1
                npt_prop = npt + 1
                if (npt_prop>malayv) then  !MODIF
                    out=0
                else
                    voro_prop(npt_prop,1) = d_min+ran3(ra)*(d_max-d_min)
                    call whichcell_d(voro_prop(npt_prop,1),voro,npt,ind)!

                    voro_prop(npt_prop,2) = voro(ind,2)+gasdev(ra)*sigmav ! sigmav: special width for new layers
                    !voro_prop(npt_prop,2) = -width+2*width*ran3(ra) ! use completely random new value
                    voro_prop(npt_prop,3) = -1
                    !voro_prop(npt_prop,4) = voro(ind,4)+gasdev(ra)*p_vp
                    voro_prop(npt_prop,4) = vp_min+(vp_max-vp_min)*ran3(ra)
                    isoflag_prop(npt_prop) = .true.

                    logprob_vsv=log(1/(sigmav*sqrt(2*pi)))-((voro(ind,2)-voro_prop(npt_prop,2))**2)/(2*sigmav**2) ! correct acceptance rates because transdimensional
                    !logprob_vp=log(1/(p_vp*sqrt(2*pi)))-((voro(ind,4)-voro_prop(npt_prop,4))**2)/(2*p_vp**2)

                    !Check bounds
                    if ((voro_prop(npt_prop,2)<=-width).or.(voro_prop(npt_prop,2)>=width)) then
                        out=0
                    end if
                    if ((voro_prop(npt_prop,4)<=vp_min).or.(voro_prop(npt_prop,4)>=vp_max)) then
                        out=0
                    end if
                endif
            elseif (u<0.8) then !death an isotropic cell !---------------------------------------    !

                death = .true.
                PrD = PrD + 1
                if(npt_iso==0) then
                    out=0
                else

                    ind2 = ceiling(ran3(ra)*(npt_iso))
                    if (ind2==0) ind2=1

                    j=0
                    do ii=1,npt
                        if (isoflag(ii).eqv..true.) j=j+1
                        if (j==ind2) then
                            ind=ii
                            exit
                        endif
                    enddo
                    if (voro(ind,3).NE.-1) stop "1092"
                endif

                npt_prop=npt-1

                if ((npt_prop<milay).or.(ind==1)) then
                    out=0
                else
                    voro_prop(ind,:)=voro(npt,:)
                    isoflag_prop(ind)=isoflag(npt)

                    call whichcell_d(voro(ind,1),voro_prop,npt_prop,ind2)
                    logprob_vsv=log(1/(sigmav*sqrt(2*pi)))-((voro(ind,2)-voro_prop(ind2,2))**2)/(2*sigmav**2) ! same as for birth
                    !logprob_vp=log(1/(p_vp*sqrt(2*pi)))-((voro(ind,4)-voro_prop(ind2,4))**2)/(2*p_vp**2)

                endif
            elseif (u<0.9) then !Birth an anisotropic layer----------------------------------------
                birtha = .true.
                PrBa = PrBa + 1
                if (npt_iso==0) then
                    out=0
                else
                    ! Choose randomly an isotropic cell among the npt_iso possibilities
                    ind2 = ceiling(ran3(ra)*(npt_iso))
                    if (ind2==0) ind2=1
                    j=0
                    do ii=1,npt
                        if (isoflag(ii).eqv..true.) j=j+1
                        if (j==ind2) then
                            ind=ii
                            exit
                        endif
                    enddo
                    voro_prop(ind,3) =  xi_min+(xi_max-xi_min)*ran3(ra)
                    if (voro(ind,3).NE.-1) stop "1130"
                    if ((voro_prop(ind,3)<=xi_min).or.(voro_prop(ind,3)>=xi_max)) then
                        write(*,*)'anisotropy out of bounds'
                        out=0
                    endif
                    isoflag_prop(ind)=.false.
                endif

            !else
            elseif (u<1.) then !death of an anisotropic layer!---------------------------------------
                deatha = .true.
                PrDa = PrDa + 1
                if (npt_ani==0) then
                    out=0
                else
                    ! Choose randomly an anisotropic cell among the npt_ani possibilities
                    ind2 = ceiling(ran3(ra)*(npt_ani))
                    if (ind2==0) ind2=1
                    j=0
                    do ii=1,npt
                        if (isoflag(ii).eqv..false.) j=j+1
                        if (j==ind2) then
                            ind=ii
                            exit
                        endif
                    enddo
                    voro_prop(ind,3)=-1
                    isoflag_prop(ind)=.true.
                endif
            else
                out=0

            endif

            !**************************************************************************

            !         After  moving a cell, Get the misfit of the proposed model

            !**************************************************************************
            if (out==1) then ! maybe we don't need a forward calculation for changes in noise parameters

                call combine_linear_vp(model_ref,nptref,nic_ref,noc_ref,voro_prop,npt_prop,d_max,&
                    r,rho,vpv,vph,vsv,vsh,qkappa,qshear,eta,nptfinal,nic,noc,vs_data,xi,vp_data)

                if (ndatad_R>0) then

                    jcom=3 !rayleigh waves

                    do iharm=1,nharm_R
                        nmodes=0
                        call minos_bran(1,tref,nptfinal,nic,noc,r,rho,vpv,vph,vsv,vsh,&
                            qkappa,qshear,eta,eps,wgrav,jcom,lmin,lmax,&
                            wmin_R(iharm),wmax_R(iharm),numharm_R(iharm),numharm_R(iharm),&
                            nmodes_max,nmodes,n_mode,l_mode,c_ph,period,raylquo,error_flag) ! calculate normal modes
                        if (error_flag) then
                            out=0
                            !write(*,*)"INVALID PROPOSED MODEL - RAYLEIGH - minos_bran.f FAIL 005"
                            goto 1142
                        endif

                        peri_R_tmp(:nlims_R(2,iharm)-nlims_R(1,iharm)+1)=peri_R(nlims_R(1,iharm):nlims_R(2,iharm))
                        n_R_tmp(:nlims_R(2,iharm)-nlims_R(1,iharm)+1)=n_R(nlims_R(1,iharm):nlims_R(2,iharm))
                        ndatad_R_tmp=nlims_R(2,iharm)-nlims_R(1,iharm)+1 ! fortran slices take the first and the last element
                        call dispersion_minos(nmodes_max,nmodes,n_mode,c_ph,period,raylquo,&
                            peri_R_tmp,n_R_tmp,d_cR_prop_tmp,rq_R,ndatad_R_tmp,ier) ! extract phase velocities from minos output (pretty ugly)
                        !write(*,*)"dispersion_minos for RAYLEIGH"
                        d_cR_prop(nlims_R(1,iharm):nlims_R(2,iharm))=d_cR_prop_tmp(:nlims_R(2,iharm)-nlims_R(1,iharm)+1)

                        if (ier) then
                            out=0
                            !write(*,*)"INVALID PROPOSED MODEL - RAYLEIGH..."
                            goto 1142
                        endif

                        if (maxval(abs(rq_R(:ndatad_R_tmp)))>maxrq*eps) then
                            out=0
                            !write(*,*)"BAD UNDERTAINTIES - RAYLEIGH..."
                            goto 1142
                        endif
                    enddo
                endif

                if (ndatad_L>0) then

                    jcom=2 !love waves
                    do iharm=1,nharm_L
                        nmodes=0
                        call minos_bran(1,tref,nptfinal,nic,noc,r,rho,vpv,vph,vsv,vsh,&
                            qkappa,qshear,eta,eps,wgrav,jcom,lmin,lmax,&
                            wmin_L(iharm),wmax_L(iharm),numharm_L(iharm),numharm_L(iharm),&
                            nmodes_max,nmodes,n_mode,l_mode,c_ph,period,raylquo,error_flag) ! calculate normal modes
                        if (error_flag) then
                            out=0
                            !write(*,*)"INVALID PROPOSED MODEL - LOVE - minos_bran.f FAIL 006"
                            goto 1142
                        endif

                        peri_L_tmp(:nlims_L(2,iharm)-nlims_L(1,iharm)+1)=peri_L(nlims_L(1,iharm):nlims_L(2,iharm))
                        n_L_tmp(:nlims_L(2,iharm)-nlims_L(1,iharm)+1)=n_L(nlims_L(1,iharm):nlims_L(2,iharm))
                        ndatad_L_tmp=nlims_L(2,iharm)-nlims_L(1,iharm)+1 ! fortran slices take the first and the last element
                        call dispersion_minos(nmodes_max,nmodes,n_mode,c_ph,period,raylquo,&
                            peri_L_tmp,n_L_tmp,d_cL_prop_tmp,rq_L,ndatad_L_tmp,ier) ! extract phase velocities from minos output (pretty ugly)
                        !write(*,*)"dispersion_minos for LOVE"
                        d_cL_prop(nlims_L(1,iharm):nlims_L(2,iharm))=d_cL_prop_tmp(:nlims_L(2,iharm)-nlims_L(1,iharm)+1)

                        if (ier) then
                            out=0
                            !write(*,*)"INVALID PROPOSED MODEL - LOVE..."
                            goto 1142
                        endif

                        if (maxval(abs(rq_L(:ndatad_L_tmp)))>maxrq*eps) then
                            out=0
                            !write(*,*)"BAD UNDERTAINTIES - LOVE..."
                            goto 1142
                        endif
                    enddo
                endif

                lsd_R_prop=0
                lsd_L_prop=0
                liked_R_prop=0
                liked_L_prop=0
                do i=1,ndatad_R
                    lsd_R_prop=lsd_R_prop+(d_obsdCR(i)-d_cR_prop(i))**2
                    liked_R_prop=liked_R_prop+(d_obsdCR(i)-d_cR_prop(i))**2/(2*(Ad_R_prop*d_obsdCRe(i)*(d_obsdCR(i)/100))**2) ! prendre en compte erreurs mesurées
                end do
                do i=1,ndatad_L
                    lsd_L_prop=lsd_L_prop+(d_obsdCL(i)-d_cL_prop(i))**2
                    liked_L_prop=liked_L_prop+(d_obsdCL(i)-d_cL_prop(i))**2/(2*(Ad_L_prop*d_obsdCLe(i)*(d_obsdCL(i)/100))**2) ! prendre en compte erreurs mesurées
                end do

            endif

            like_prop=(liked_R_prop+liked_L_prop) !log-likelihood of the proposed model
            like_prop_w=like_prop/widening
            like_w=like/widening

            !like_prop=1
            !like_prop_w=like_prop/widening

       1142 Accept = .false.

            ! now check if we accept the new model - different treatement depending on the change
            if (birth) then!------------------------------------------------------------------

                if (log(ran3(ra))<log(out) + log(real(npt)+1)-log(real(npt_prop)+1) -&
                    log(2*width)-logprob_vsv&!-log(vpvsv_max-vpvsv_min)-logprob_vp&
                    -like_prop_w+like_w) then ! transdimensional case
                    accept=.true.
                    AcB=AcB+1
                endif

            elseif (death) then!-------------------------------------------------------------

                if (log(ran3(ra))<log(out) + log(real(npt)+1)&
                    -log(real(npt_prop)+1) + &
                    log(2*width)+logprob_vsv&!+log(vpvsv_max-vpvsv_min)+logprob_vp&
                    -like_prop_w+like_w) then! transdimensional case
                    accept=.true.
                    AcD=AcD+1
                endif

            elseif (noisd_R) then !@@@@@@@@@@@@@@@@@@@@@@@@@@@@ logrsig  @@@@@@@@@@@@
                logrsig=ndatad_R*log(Ad_R/Ad_R_prop)/widening ! ATTENTION avc ld 2
                if (log(ran3(ra))<logrsig+log(out)-like_prop_w+like_w) then ! hierarchical case
                    accept=.true.
                    Acnd_R=Acnd_R+1
                endif
            elseif (noisd_L) then !@@@@@@@@@@@@@@@@@@@@@@@@@@@@ logrsig  @@@@@@@@@@@@
                logrsig=ndatad_L*log(Ad_L/Ad_L_prop)/widening ! ATTENTION avc ld 2
                if (log(ran3(ra))<logrsig+log(out)-like_prop_w+like_w) then ! hierarchical case
                    accept=.true.
                    Acnd_L=Acnd_L+1
                endif
            else !NO JUMP-------------------------------------------------------------------

                if (log(ran3(ra))<log(out)-like_prop_w+like_w)then
                    if (ind>malayv) stop  '1082'
                    accept=.true.
                    if (value) then
                        AcV=AcV+1
                    elseif (move) then
                        AcP=AcP+1
                    elseif(ani)then
                        Acxi=Acxi+1
                    elseif(change_vp)then
                        Ac_vp=Ac_vp+1
                    endif
                    if (birtha) then
                        Acba=Acba+1
                    elseif (deatha) then
                        Acda=Acda+1
                    endif
                endif!accept
            endif

            !***********************************************************************************
            !   If we accept the proposed model, update the status of the Markov Chain
            !***********************************************************************************

            if (accept) then

                isoflag=isoflag_prop
                voro=voro_prop
                like=like_prop
                like_w=like_prop_w

                liked_L=liked_L_prop
                liked_R=liked_R_prop
                lsd_L=lsd_L_prop
                lsd_R=lsd_R_prop
                npt=npt_prop
                Ad_R=Ad_R_prop
                Ad_L=Ad_L_prop

                d_cR=d_cR_prop
                d_cL=d_cL_prop

            endif

            npt_ani=0
            npt_iso=0

            do i=1,npt
                if (isoflag(i).eqv..true.) npt_iso=npt_iso+1
                if (isoflag(i).eqv..false.) npt_ani=npt_ani+1
            enddo
            if (npt_iso+npt_ani.ne.npt) stop "Error here"

            ! check for burnin
            if (burnin_in_progress) then
                !ount=ount+1


                IF ((mod(ount,display).EQ.0)) THEN !.and.(mod(ran,50).EQ.0)

                    write(*,*)'processor number',ran+1,'/',nbproc
                    write(*,*)'widening step:','burn-in'
                    write(*,*)'sample:',ount,'/',burn_in
                    write(*,*)'number of cells:',npt
                    write(*,*)'Ad_R',Ad_R,'Ad_L',Ad_L
                    write(*,*)'Acceptance rates'
                    write(*,*)'AR_move',100*AcP/PrP
                    write(*,*)'AR_value',100*AcV/PrV

                    write(*,*)'AR_Birth',100*AcB/PrB,'AR_Death',100*AcD/PrD
                    write(*,*)'sigmav',sigmav
                    write(*,*)'AR_Birtha',100*AcBa/PrBa,'AR_Deatha',100*AcDa/PrDa
                    write(*,*)'AR_xi',100*Acxi/Prxi,'pxi',pxi
                    write(*,*)'AR__vp',100*Ac_vp/Pr_vp,'p_vp',p_vp
                    write(*,*)'AR_Ad_R',100*Acnd_R/Prnd_R,'pAd_R',pAd_R
                    write(*,*)'AR_Ad_L',100*Acnd_L/Prnd_L,'pAd_L',pAd_L
                    write(*,*)'npt_iso',npt_iso,'npt_ani',npt_ani
                    write(*,*)'widening',widening
                    !write(*,*)Ar_birth_old,sigmav_old,sigmav
                    write(*,*)'-----------------------------------------'
                    write(*,*)
                    write(*,*)


                END IF

                !IF (mod(ount,thin)==0) THEN
                convBa_bi(ount)=100*AcBa/PrBa ! acceptance rates in percent
                convB_bi(ount)=100*AcB/PrB
                convDa_bi(ount)=100*AcDa/PrDa
                convD_bi(ount)=100*AcD/PrD
                convvp_bi(ount)=100*Ac_vp/Pr_vp
                convvs_bi(ount)=100*Acv/Prv
                convdp_bi(ount)=100*Acp/Prp
                convxi_bi(ount)=100*Acxi/Prxi
                convd_R_bi(ount)=lsd_R
                convd_L_bi(ount)=lsd_L
                ncell_bi(ount)=npt
                convAd_R_bi(ount)=Ad_R
                convAd_L_bi(ount)=Ad_L
                !endif

                if (ount>burn_in) then
                    write(*,*)'burn-in done',rank
                    burnin_in_progress=.false.
                    ount=0
                    call MPI_BARRIER(MPI_COMM_WORLD, ierror)
                    write(*,*)'barrier reached',rank
                    call MPI_REDUCE(convBa_bi,convBas_bi,burn_in,&
                    MPI_Real,MPI_Sum,0,MPI_COMM_WORLD,ierror)
                    call MPI_REDUCE(convB_bi,convBs_bi,burn_in,&
                    MPI_Real,MPI_Sum,0,MPI_COMM_WORLD,ierror)
                    call MPI_REDUCE(convDa_bi,convDas_bi,burn_in,&
                    MPI_Real,MPI_Sum,0,MPI_COMM_WORLD,ierror)
                    call MPI_REDUCE(convD_bi,convDs_bi,burn_in,&
                    MPI_Real,MPI_Sum,0,MPI_COMM_WORLD,ierror)
                    call MPI_REDUCE(convvp_bi,convvps_bi,burn_in,&
                    MPI_Real,MPI_Sum,0,MPI_COMM_WORLD,ierror)
                    call MPI_REDUCE(convvs_bi,convvss_bi,burn_in,&
                    MPI_Real,MPI_Sum,0,MPI_COMM_WORLD,ierror)
                    call MPI_REDUCE(convdp_bi,convdps_bi,burn_in,&
                    MPI_Real,MPI_Sum,0,MPI_COMM_WORLD,ierror)
                    call MPI_REDUCE(convxi_bi,convxis_bi,burn_in,&
                    MPI_Real,MPI_Sum,0,MPI_COMM_WORLD,ierror)
                    call MPI_REDUCE(convd_R_bi,convd_Rs_bi,burn_in,&
                    MPI_Real,MPI_Sum,0,MPI_COMM_WORLD,ierror)
                    call MPI_REDUCE(convd_L_bi,convd_Ls_bi,burn_in,&
                    MPI_Real,MPI_Sum,0,MPI_COMM_WORLD,ierror)
                    call MPI_REDUCE(convAd_R_bi,convAd_Rs_bi,burn_in,&
                    MPI_Real,MPI_Sum,0,MPI_COMM_WORLD,ierror)
                    call MPI_REDUCE(convAd_L_bi,convAd_Ls_bi,burn_in,&
                    MPI_Real,MPI_Sum,0,MPI_COMM_WORLD,ierror)
                    call MPI_REDUCE(ncell_bi,ncells_bi,burn_in,&
                    MPI_Real,MPI_Sum,0,MPI_COMM_WORLD,ierror)
                    call MPI_REDUCE(convsigmav_bi,convsigmavs_bi,burn_in,&
                    MPI_Real,MPI_Sum,0,MPI_COMM_WORLD,ierror)

                    if (rank==0) then
                        convPs_bi=convPs_bi/nbproc
                        convBs_bi=convBs_bi/nbproc
                        convBas_bi=convBas_bi/nbproc
                        convDs_bi=convDs_bi/nbproc
                        convDas_bi=convDas_bi/nbproc
                        convvps_bi=convvps_bi/nbproc
                        convvss_bi=convvss_bi/nbproc
                        convdps_bi=convdps_bi/nbproc
                        convxis_bi=convxis_bi/nbproc
                        convd_Rs_bi=convd_Rs_bi/nbproc
                        convd_Ls_bi=convd_Ls_bi/nbproc
                        convAd_Rs_bi=convAd_Rs_bi/nbproc
                        convAd_Ls_bi=convAd_Ls_bi/nbproc
                        convsigmavs_bi=convsigmavs_bi/nbproc
                        ncells_bi=ncells_bi/nbproc
                    endif

                    call MPI_BARRIER(MPI_COMM_WORLD, ierror)
                    write(*,*)'reducing done',rank
                    IF (rank==0) THEN

                        write(filenametmp,"('/Convergence_Birth_burn-in_',I3.3,'.out')")rank
                        open(54,file=dirname//filenametmp,status='replace')
                        write(54,*)burn_in,n_w,burn_in_widening,nsample_widening
                        do i=1,burn_in
                            write(54,*)convB_bi(i),convBa_bi(i)
                        enddo
                        close(54)

                        write(filenametmp,"('/Convergence_Death_burn-in_',I3.3,'.out')")rank
                        open(54,file=dirname//filenametmp,status='replace')
                        write(54,*)burn_in,n_w,burn_in_widening,nsample_widening
                        do i=1,burn_in
                            write(54,*)convD_bi(i),convDa_bi(i)
                        enddo
                        close(54)

                        write(filenametmp,"('/Convergence_vp_burn-in_',I3.3,'.out')")rank
                        open(54,file=dirname//filenametmp,status='replace')
                        write(54,*)burn_in,n_w,burn_in_widening,nsample_widening
                        do i=1,burn_in
                            write(54,*)convvp_bi(i)
                        enddo
                        close(54)

                        write(filenametmp,"('/Convergence_vs_burn-in_',I3.3,'.out')")rank
                        open(54,file=dirname//filenametmp,status='replace')
                        write(54,*)burn_in,n_w,burn_in_widening,nsample_widening
                        do i=1,burn_in
                            write(54,*)convvs_bi(i)
                        enddo
                        close(54)

                        write(filenametmp,"('/Convergence_dp_burn-in_',I3.3,'.out')")rank
                        open(54,file=dirname//filenametmp,status='replace')
                        write(54,*)burn_in,n_w,burn_in_widening,nsample_widening
                        do i=1,burn_in
                            write(54,*)convdp(i)
                        enddo
                        close(54)

                        write(filenametmp,"('/Convergence_xi_burn-in_',I3.3,'.out')")rank
                        open(54,file=dirname//filenametmp,status='replace')
                        write(54,*)burn_in,n_w,burn_in_widening,nsample_widening
                        do i=1,burn_in
                            write(54,*)convxi(i)
                        enddo
                        close(54)

                        write(filenametmp,"('/Convergence_misfit_burn-in_',I3.3,'.out')")rank
                        open(54,file=dirname//filenametmp,status='replace')
                        write(54,*)burn_in,n_w,burn_in_widening,nsample_widening
                        do i=1,burn_in
                            write(54,*)convd_R_bi(i),convd_L_bi(i)
                        enddo
                        close(54)

                        write(filenametmp,"('/Convergence_nb_burn-in_',I3.3,'.out')")rank
                        open(53,file=dirname//filenametmp,status='replace')
                        write(53,*)burn_in,n_w,burn_in_widening,nsample_widening
                        do i=1,burn_in
                            write(53,*)ncell_bi(i)
                        enddo
                        close(53)

                        write(filenametmp,"('/Convergence_sigma_R_burn-in_',I3.3,'.out')")rank
                        open(52,file=dirname//filenametmp,status='replace')
                        write(52,*)burn_in,n_w,burn_in_widening,nsample_widening
                        do i=1,burn_in
                            write(52,*)convAd_R_bi(i)
                        enddo
                        close(52)

                        write(filenametmp,"('/Convergence_sigma_L_burn-in_',I3.3,'.out')")rank
                        open(53,file=dirname//filenametmp,status='replace')
                        write(53,*)burn_in,n_w,burn_in_widening,nsample_widening
                        do i=1,burn_in
                            write(53,*)convAd_L_bi(i)
                        enddo
                        close(53)

                        write(filenametmp,"('/Convergence_sigmav_burn-in_',I3.3,'.out')")rank
                        open(53,file=dirname//filenametmp,status='replace')
                        write(53,*)burn_in,n_w,burn_in_widening,nsample_widening
                        do i=1,burn_in
                            write(53,*)convsigmav_bi(i)
                        enddo
                        close(53)

                        write(filenametmp,"('/Convergence_Birth_burn-in_ref.out')")
                        open(54,file=dirname//filenametmp,status='replace')
                        write(54,*)burn_in,n_w,burn_in_widening,nsample_widening
                        do i=1,burn_in
                            write(54,*)convBs_bi(i),convBas_bi(i)
                        enddo
                        close(54)

                        write(filenametmp,"('/Convergence_Death_burn-in_ref.out')")
                        open(54,file=dirname//filenametmp,status='replace')
                        write(54,*)burn_in,n_w,burn_in_widening,nsample_widening
                        do i=1,burn_in
                            write(54,*)convDs_bi(i),convDas_bi(i)
                        enddo
                        close(54)

                        write(filenametmp,"('/Convergence_vp_burn-in_ref.out')")
                        open(54,file=dirname//filenametmp,status='replace')
                        write(54,*)burn_in,n_w,burn_in_widening,nsample_widening
                        do i=1,burn_in
                            write(54,*)convvps_bi(i)
                        enddo
                        close(54)

                        write(filenametmp,"('/Convergence_vs_burn-in_ref.out')")
                        open(54,file=dirname//filenametmp,status='replace')
                        write(54,*)burn_in,n_w,burn_in_widening,nsample_widening
                        do i=1,burn_in
                            write(54,*)convvss_bi(i)
                        enddo
                        close(54)

                        write(filenametmp,"('/Convergence_dp_burn-in_ref.out')")
                        open(54,file=dirname//filenametmp,status='replace')
                        write(54,*)burn_in,n_w,burn_in_widening,nsample_widening
                        do i=1,burn_in
                            write(54,*)convdps(i)
                        enddo
                        close(54)

                        write(filenametmp,"('/Convergence_xi_burn-in_ref.out')")
                        open(54,file=dirname//filenametmp,status='replace')
                        write(54,*)burn_in,n_w,burn_in_widening,nsample_widening
                        do i=1,burn_in
                            write(54,*)convxis(i)
                        enddo
                        close(54)

                        write(filenametmp,"('/Convergence_misfit_burn-in_ref.out')")
                        open(54,file=dirname//filenametmp,status='replace')
                        write(54,*)burn_in,n_w,burn_in_widening,nsample_widening
                        do i=1,burn_in
                            write(54,*)convd_Rs_bi(i),convd_Ls_bi(i)
                        enddo
                        close(54)

                        write(filenametmp,"('/Convergence_nb_burn-in_ref.out')")
                        open(53,file=dirname//filenametmp,status='replace')
                        write(53,*)burn_in,n_w,burn_in_widening,nsample_widening
                        do i=1,burn_in
                            write(53,*)ncells_bi(i)
                        enddo
                        close(53)

                        write(filenametmp,"('/Convergence_sigma_R_burn-in_ref.out')")
                        open(52,file=dirname//filenametmp,status='replace')
                        write(52,*)burn_in,n_w,burn_in_widening,nsample_widening
                        do i=1,burn_in
                            write(52,*)convAd_Rs_bi(i)
                        enddo
                        close(52)

                        write(filenametmp,"('/Convergence_sigma_L_burn-in_ref.out')")
                        open(53,file=dirname//filenametmp,status='replace')
                        write(53,*)burn_in,n_w,burn_in_widening,nsample_widening
                        do i=1,burn_in
                            write(53,*)convAd_Ls_bi(i)
                        enddo
                        close(53)

                        write(filenametmp,"('/Convergence_sigmav_burn-in_ref.out')")
                        open(53,file=dirname//filenametmp,status='replace')
                        write(53,*)burn_in,n_w,burn_in_widening,nsample_widening
                        do i=1,burn_in
                            write(53,*)convsigmavs_bi(i)
                        enddo
                        close(53)

                    else ! only print convergence diagnostics for the core
                        write(filenametmp,"('/Convergence_Birth_burn-in_',I3.3,'.out')")rank
                        open(54,file=dirname//filenametmp,status='replace')
                        write(54,*)burn_in,n_w,burn_in_widening,nsample_widening
                        do i=1,burn_in
                            write(54,*)convB_bi(i),convBa_bi(i)
                        enddo
                        close(54)

                        write(filenametmp,"('/Convergence_Death_burn-in_',I3.3,'.out')")rank
                        open(54,file=dirname//filenametmp,status='replace')
                        write(54,*)burn_in,n_w,burn_in_widening,nsample_widening
                        do i=1,burn_in
                            write(54,*)convD_bi(i),convDa_bi(i)
                        enddo
                        close(54)

                        write(filenametmp,"('/Convergence_vp_burn-in_',I3.3,'.out')")rank
                        open(54,file=dirname//filenametmp,status='replace')
                        write(54,*)burn_in,n_w,burn_in_widening,nsample_widening
                        do i=1,burn_in
                            write(54,*)convvp_bi(i)
                        enddo
                        close(54)

                        write(filenametmp,"('/Convergence_vs_burn-in_',I3.3,'.out')")rank
                        open(54,file=dirname//filenametmp,status='replace')
                        write(54,*)burn_in,n_w,burn_in_widening,nsample_widening
                        do i=1,burn_in
                            write(54,*)convvs_bi(i)
                        enddo
                        close(54)

                        write(filenametmp,"('/Convergence_dp_burn-in_',I3.3,'.out')")rank
                        open(54,file=dirname//filenametmp,status='replace')
                        write(54,*)burn_in,n_w,burn_in_widening,nsample_widening
                        do i=1,burn_in
                            write(54,*)convdp(i)
                        enddo
                        close(54)

                        write(filenametmp,"('/Convergence_xi_burn-in_',I3.3,'.out')")rank
                        open(54,file=dirname//filenametmp,status='replace')
                        write(54,*)burn_in,n_w,burn_in_widening,nsample_widening
                        do i=1,burn_in
                            write(54,*)convxi(i)
                        enddo
                        close(54)

                        write(filenametmp,"('/Convergence_misfit_burn-in_',I3.3,'.out')")rank
                        open(54,file=dirname//filenametmp,status='replace')
                        write(54,*)burn_in,n_w,burn_in_widening,nsample_widening
                        do i=1,burn_in
                            write(54,*)convd_R_bi(i),convd_L_bi(i)
                        enddo
                        close(54)

                        write(filenametmp,"('/Convergence_nb_burn-in_',I3.3,'.out')")rank
                        open(53,file=dirname//filenametmp,status='replace')
                        write(53,*)burn_in,n_w,burn_in_widening,nsample_widening
                        do i=1,burn_in
                            write(53,*)ncell_bi(i)
                        enddo
                        close(53)

                        write(filenametmp,"('/Convergence_sigma_R_burn-in_',I3.3,'.out')")rank
                        open(52,file=dirname//filenametmp,status='replace')
                        write(52,*)burn_in,n_w,burn_in_widening,nsample_widening
                        do i=1,burn_in
                            write(52,*)convAd_R_bi(i)
                        enddo
                        close(52)

                        write(filenametmp,"('/Convergence_sigma_L_burn-in_',I3.3,'.out')")rank
                        open(53,file=dirname//filenametmp,status='replace')
                        write(53,*)burn_in,n_w,burn_in_widening,nsample_widening
                        do i=1,burn_in
                            write(53,*)convAd_L_bi(i)
                        enddo
                        close(53)

                        write(filenametmp,"('/Convergence_sigmav_burn-in_',I3.3,'.out')")rank
                        open(53,file=dirname//filenametmp,status='replace')
                        write(53,*)burn_in,n_w,burn_in_widening,nsample_widening
                        do i=1,burn_in
                            write(53,*)convsigmav_bi(i)
                        enddo
                        close(53)

                    endif
                    write(*,*)'writing done',rank
                endif

                cycle ! continue the loop skipping anything below
            endif
            !****************************************************************

            !                  Store models for ensemble solution

            !****************************************************************
            ! Store models for restart
            IF ((mod(ount,store)==0)) THEN
                write(number,1000)rank
                filename=number//'model.dat'
                open(rank*100,file=storename//'/'//filename,status='replace')
                write(rank*100,*)nbproc
                write(rank*100,*)npt
                do i=1,npt
                    write(rank*100,*)voro(i,:)
                enddo
                close(rank*100)
            ENDIF

            IF ((ount.GT.burn_in_widening).and.(.not.(burnin_in_progress))) THEN
                sample=sample+1

                IF (mod(ount,thin)==0) THEN

                    th = th + 1
                    th_all=th_all+1

                    call combine_linear_vp(model_ref,nptref,nic_ref,noc_ref,voro,npt,d_max,&
                        r,rho,vpv,vph,vsv,vsh,qkappa,qshear,eta,nptfinal,nic,noc,vs_data,xi,vp_data)

                    !write(*,*)th-1,nbproc,rank,(th-1)*nbproc,everyall,filenamemax

                    if (((th-1)*nbproc>=(num_file+1)*everyall)) then !.and.(th>0) ! if enough models in current files, create a new file
                    call MPI_BARRIER(MPI_COMM_WORLD, ierror)
                    call MPI_FILE_CLOSE(filehandle,ierror)

                    num_file=num_file+1

                    write(*,*)num_file,widening
                    write(filenamemax,"('/All_models_prepare_',I3.3,'_',f5.2,'.out')")num_file,widening
                    write(*,*)filenamemax

                    call MPI_FILE_OPEN(MPI_COMM_WORLD,dirname//filenamemax,MPI_MODE_WRONLY + MPI_MODE_CREATE, &
                    MPI_INFO_NULL,filehandle,ierror) ! check MPI_INFO_NULL: might lead to conflicting filehandles
                    call MPI_BARRIER(MPI_COMM_WORLD, ierror)

                    if (rank==0) then ! create header

                        position_file=0

                        write(*,*)'writing header of file',filenamemax

                        call MPI_FILE_WRITE_AT(filehandle,position_file,nbproc,1,MPI_INTEGER, &
                        status_write_mpi,ierror)!sample/everyall/thin*nbproc
                        position_file=position_file+nb_bytes_integer

                        call MPI_FILE_WRITE_AT(filehandle,position_file,everyall,1,MPI_INTEGER, &
                        status_write_mpi,ierror)
                        position_file=position_file+nb_bytes_integer

                        call MPI_FILE_WRITE_AT(filehandle,position_file,burn_in,1,MPI_INTEGER, &
                        status_write_mpi,ierror)
                        position_file=position_file+nb_bytes_integer
                        call MPI_FILE_WRITE_AT(filehandle,position_file,widening,1,MPI_REAL, &
                        status_write_mpi,ierror)
                        position_file=position_file+nb_bytes_real
                        call MPI_FILE_WRITE_AT(filehandle,position_file,thin,1,MPI_INTEGER, &
                        status_write_mpi,ierror)
                        position_file=position_file+nb_bytes_integer

                        call MPI_FILE_WRITE_AT(filehandle,position_file,d_min,1,MPI_REAL, &
                        status_write_mpi,ierror)
                        position_file=position_file+nb_bytes_real
                        call MPI_FILE_WRITE_AT(filehandle,position_file,d_max,1,MPI_REAL, &
                        status_write_mpi,ierror)
                        position_file=position_file+nb_bytes_real

                        call MPI_FILE_WRITE_AT(filehandle,position_file,width,1,MPI_REAL, &
                        status_write_mpi,ierror)
                        position_file=position_file+nb_bytes_real

                        call MPI_FILE_WRITE_AT(filehandle,position_file,xi_min,1,MPI_REAL, &
                        status_write_mpi,ierror)
                        position_file=position_file+nb_bytes_real
                        call MPI_FILE_WRITE_AT(filehandle,position_file,xi_max,1,MPI_REAL, &
                        status_write_mpi,ierror)
                        position_file=position_file+nb_bytes_real

                        call MPI_FILE_WRITE_AT(filehandle,position_file,vp_min,1,MPI_REAL, &
                        status_write_mpi,ierror)
                        position_file=position_file+nb_bytes_real
                        call MPI_FILE_WRITE_AT(filehandle,position_file,vp_max,1,MPI_REAL, &
                        status_write_mpi,ierror)
                        position_file=position_file+nb_bytes_real

                        call MPI_FILE_WRITE_AT(filehandle,position_file,Ad_R_min,1,MPI_DOUBLE_PRECISION, &
                        status_write_mpi,ierror)
                        position_file=position_file+nb_bytes_double
                        call MPI_FILE_WRITE_AT(filehandle,position_file,Ad_R_max,1,MPI_DOUBLE_PRECISION, &
                        status_write_mpi,ierror)
                        position_file=position_file+nb_bytes_double

                        call MPI_FILE_WRITE_AT(filehandle,position_file,Ad_L_min,1,MPI_DOUBLE_PRECISION, &
                        status_write_mpi,ierror)
                        position_file=position_file+nb_bytes_double
                        call MPI_FILE_WRITE_AT(filehandle,position_file,Ad_L_max,1,MPI_DOUBLE_PRECISION, &
                        status_write_mpi,ierror)
                        position_file=position_file+nb_bytes_double

                        call MPI_FILE_WRITE_AT(filehandle,position_file,milay,1,MPI_INTEGER, &
                        status_write_mpi,ierror)
                        position_file=position_file+nb_bytes_integer
                        call MPI_FILE_WRITE_AT(filehandle,position_file,malay,1,MPI_INTEGER, &
                        status_write_mpi,ierror)
                        position_file=position_file+nb_bytes_integer

                        call MPI_FILE_WRITE_AT(filehandle,position_file,mk,1,MPI_INTEGER, &
                        status_write_mpi,ierror)
                        position_file=position_file+nb_bytes_integer
                        call MPI_FILE_WRITE_AT(filehandle,position_file,ndatadmax,1,MPI_INTEGER, &
                        status_write_mpi,ierror)
                        position_file=position_file+nb_bytes_integer
                    endif

                    call MPI_BARRIER(MPI_COMM_WORLD, ierror)

                    endif

                    position_file=nb_bytes_header+mod(((th-1)*nbproc+rank),everyall)&
                        *nb_bytes_model ! position in the file where the model is written

                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    ! Write the current model
                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                    call MPI_FILE_WRITE_AT(filehandle,position_file,nptfinal,1,MPI_INTEGER, &
                    status_write_mpi,ierror)
                    position_file=position_file+nb_bytes_integer

                    call MPI_FILE_WRITE_AT(filehandle,position_file,nic,1,MPI_INTEGER, &
                    status_write_mpi,ierror)
                    position_file=position_file+nb_bytes_integer

                    call MPI_FILE_WRITE_AT(filehandle,position_file,noc,1,MPI_INTEGER, &
                    status_write_mpi,ierror)
                    position_file=position_file+nb_bytes_integer

                    call MPI_FILE_WRITE_AT(filehandle,position_file,npt,1,MPI_INTEGER, &
                    status_write_mpi,ierror)
                    position_file=position_file+nb_bytes_integer

                    call MPI_FILE_WRITE_AT(filehandle,position_file,npt_ani,1,MPI_INTEGER, &
                    status_write_mpi,ierror)
                    position_file=position_file+nb_bytes_integer

                    call MPI_FILE_WRITE_AT(filehandle,position_file,Ad_R,1,MPI_DOUBLE_PRECISION, &
                    status_write_mpi,ierror)
                    position_file=position_file+nb_bytes_double

                    call MPI_FILE_WRITE_AT(filehandle,position_file,Ad_L,1,MPI_DOUBLE_PRECISION, &
                    status_write_mpi,ierror)
                    position_file=position_file+nb_bytes_double

                    call MPI_FILE_WRITE_AT(filehandle,position_file,r,mk,MPI_REAL, &
                    status_write_mpi,ierror)
                    position_file=position_file+mk*nb_bytes_real

                    call MPI_FILE_WRITE_AT(filehandle,position_file,vsv,mk,MPI_REAL, &
                    status_write_mpi,ierror)
                    position_file=position_file+mk*nb_bytes_real

                    call MPI_FILE_WRITE_AT(filehandle,position_file,xi,mk,MPI_REAL, &
                    status_write_mpi,ierror)
                    position_file=position_file+mk*nb_bytes_real

                    call MPI_FILE_WRITE_AT(filehandle,position_file,vpv,mk,MPI_REAL, &
                    status_write_mpi,ierror)
                    position_file=position_file+mk*nb_bytes_real

                    ! if add Ad_R/widening and Ad_L/widening to like_w, change postprocess_binary_outputs.f90
                    call MPI_FILE_WRITE_AT(filehandle,position_file,like_w,1,MPI_REAL, &
                    status_write_mpi,ierror)
                    position_file=position_file+nb_bytes_real

                    call MPI_FILE_WRITE_AT(filehandle,position_file,ndatad_R,1,MPI_INTEGER, &
                    status_write_mpi,ierror)
                    position_file=position_file+nb_bytes_integer

                    call MPI_FILE_WRITE_AT(filehandle,position_file,d_cR,ndatadmax,MPI_REAL, &
                    status_write_mpi,ierror)
                    position_file=position_file+ndatadmax*nb_bytes_real

                    call MPI_FILE_WRITE_AT(filehandle,position_file,ndatad_L,1,MPI_INTEGER, &
                    status_write_mpi,ierror)
                    position_file=position_file+nb_bytes_integer

                    call MPI_FILE_WRITE_AT(filehandle,position_file,d_cL,ndatadmax,MPI_REAL, &
                    status_write_mpi,ierror)
                    position_file=position_file+ndatadmax*nb_bytes_real

                    !CALL cpu_time(t2)
                    !write(*,*)'writing time (once every 20 iterations)',t2-t1

                endif
                !get convergence
                convBa(ount)=100*AcBa/PrBa ! acceptance rates in percent
                convB(ount)=100*AcB/PrB
                convDa(ount)=100*AcDa/PrDa
                convD(ount)=100*AcD/PrD
                convvp(ount)=100*Ac_vp/Pr_vp
                convvs(ount)=100*Acv/Prv
                convdp(ount)=100*Acp/Prp
                convxi(ount)=100*Acxi/Prxi
                convd_R(ount)=lsd_R
                convd_L(ount)=lsd_L
                ncell(ount)=npt
                convAd_R(ount)=Ad_R
                convAd_L(ount)=Ad_L
            endif
                !**********************************************************************

                !       Display what is going on every "Display" samples

                !**********************************************************************

            IF ((mod(ount,display).EQ.0).and.(mod(ran,10).EQ.0)) THEN
   
                write(*,*)'processor number',ran+1,'/',nbproc
                write(*,*)'widening',widening
                if (mod(ount,burn_in_widening+nsample_widening)<burn_in_widening) then
                    write(*,*)'transition burn-in'
                else
                    write(*,*)'iterating'
                endif
                write(*,*)'sample:',ount,'/',burn_in_widening+nsample_widening
                write(*,*)'number of cells:',npt
                write(*,*)'Ad_R',Ad_R,'Ad_L',Ad_L
                write(*,*)'Acceptance rates'
                write(*,*)'AR_move',100*AcP/PrP
                write(*,*)'AR_value',100*AcV/PrV
                write(*,*)'AR_Birth',100*AcB/PrB,'AR_Death',100*AcD/PrD,'sigmav',sigmav
                write(*,*)'AR_Birtha',100*AcBa/PrBa,'AR_Deatha',100*AcDa/PrDa
                write(*,*)'AR_xi',100*Acxi/Prxi,'pxi',pxi
                write(*,*)'AR__vp',100*Ac_vp/Pr_vp,'p_vp',p_vp
                write(*,*)'AR_Ad_R',100*Acnd_R/Prnd_R,'pAd_R',pAd_R
                write(*,*)'AR_Ad_L',100*Acnd_L/Prnd_L,'pAd_L',pAd_L
                write(*,*)'npt_iso',npt_iso,'npt_ani',npt_ani
                !write(*,*)Ar_birth_old,sigmav_old,sigmav
                write(*,*)'-----------------------------------------'
                write(*,*)
                write(*,*)

            endif

            !write(*,*)ount

            !endif
        end do !End Markov chain

        k=0
        if (th.ne.0) then !normalize averages
            inorout(ran+1)=1
            k=k+1
        else
            write(*,*)
            write(*,*)'rank=',ran,'th=',th

        endif



        do i=1,nbproc
            call MPI_REDUCE(inorout,inorouts,21000,MPI_Integer,MPI_Sum,i-1,MPI_COMM_WORLD,ierror)
        enddo

        write(*,*)'rank=',ran,'th=',th

        j=0
        k=0
        do i=1,nbproc
            j=j+inorouts(i)
            if (inorouts(i).ne.0) then
                k=k+1
                members(k)=i-1
            endif
        enddo

        IF (ran==0) write(*,*) 'k',k,'nbproc',nbproc
        !***************************************************************************

        ! Collect information from all the chains and average everything

        !***************************************************************************
        flag=0
        do i=1,j
            if (ran==members(i)) flag=1
        enddo


        call MPI_Group_incl(group_world, j, members, good_group, ierror)
        !IF (ran==0) write(*,*)'1'
        !write(*,*)ran,members(1:9)
        call MPI_Comm_create(MPI_COMM_WORLD, good_group, MPI_COMM_small, ierror)
        !IF (ran==0) write(*,*)'2'

        if (flag==1) then

            call MPI_REDUCE(convBa,convBas,(burn_in_widening+nsample_widening),&
            MPI_Real,MPI_Sum,0,MPI_COMM_WORLD,ierror)
            call MPI_REDUCE(convB,convBs,(burn_in_widening+nsample_widening),&
            MPI_Real,MPI_Sum,0,MPI_COMM_WORLD,ierror)
            call MPI_REDUCE(convDa,convDas,(burn_in_widening+nsample_widening),&
            MPI_Real,MPI_Sum,0,MPI_COMM_WORLD,ierror)
            call MPI_REDUCE(convD,convDs,(burn_in_widening+nsample_widening),&
            MPI_Real,MPI_Sum,0,MPI_COMM_WORLD,ierror)
            call MPI_REDUCE(convvp,convvps,(burn_in_widening+nsample_widening),&
            MPI_Real,MPI_Sum,0,MPI_COMM_WORLD,ierror)
            call MPI_REDUCE(convvs,convvss,(burn_in_widening+nsample_widening),&
            MPI_Real,MPI_Sum,0,MPI_COMM_WORLD,ierror)
            call MPI_REDUCE(convdp,convdps,(burn_in_widening+nsample_widening),&
            MPI_Real,MPI_Sum,0,MPI_COMM_WORLD,ierror)
            call MPI_REDUCE(convxi,convxis,(burn_in_widening+nsample_widening),&
            MPI_Real,MPI_Sum,0,MPI_COMM_WORLD,ierror)
            call MPI_REDUCE(convd_R,convd_Rs,(burn_in_widening+nsample_widening),&
            MPI_Real,MPI_Sum,0,MPI_COMM_WORLD,ierror)
            call MPI_REDUCE(convd_L,convd_Ls,(burn_in_widening+nsample_widening),&
            MPI_Real,MPI_Sum,0,MPI_COMM_WORLD,ierror)
            call MPI_REDUCE(convAd_R,convAd_Rs,(burn_in_widening+nsample_widening),&
            MPI_Real,MPI_Sum,0,MPI_COMM_WORLD,ierror)
            call MPI_REDUCE(convAd_L,convAd_Ls,(burn_in_widening+nsample_widening),&
            MPI_Real,MPI_Sum,0,MPI_COMM_WORLD,ierror)
            call MPI_REDUCE(ncell,ncells,(burn_in_widening+nsample_widening),&
            MPI_Real,MPI_Sum,0,MPI_COMM_WORLD,ierror)
            call MPI_REDUCE(convsigmav,convsigmavs,(burn_in_widening+nsample_widening),&
            MPI_Real,MPI_Sum,0,MPI_COMM_WORLD,ierror)

            convPs=convPs/nbproc
            convBs=convBs/nbproc
            convBas=convBas/nbproc
            convDs=convDs/nbproc
            convDas=convDas/nbproc
            convvps=convvps/nbproc
            convvss=convvss/nbproc
            convdps=convdps/nbproc
            convxis=convxis/nbproc
            convd_Rs=convd_Rs/nbproc
            convd_Ls=convd_Ls/nbproc
            convAd_Rs=convAd_Rs/nbproc
            convAd_Ls=convAd_Ls/nbproc
            convsigmavs=convsigmavs/nbproc
            ncells=ncells/nbproc

        endif

        !***********************************************************************

        !                      Write the results

        !***********************************************************************

        IF (ran==members(1)) THEN
            write(filenametmp,"('/Convergence_Birth_',f5.2'_',I3.3,'.out')")widening,rank
            open(54,file=dirname//filenametmp,status='replace')
            write(54,*)burn_in,n_w,burn_in_widening,nsample_widening
            do i=1,(burn_in_widening+nsample_widening)
                write(54,*)convB(i),convBa(i)
            enddo
            close(54)

            write(filenametmp,"('/Convergence_Death_',f5.2'_',I3.3,'.out')")widening,rank
            open(54,file=dirname//filenametmp,status='replace')
            write(54,*)burn_in,n_w,burn_in_widening,nsample_widening
            do i=1,(burn_in_widening+nsample_widening)
                write(54,*)convD(i),convDa(i)
            enddo
            close(54)

            write(filenametmp,"('/Convergence_vp_',f5.2'_',I3.3,'.out')")widening,rank
            open(54,file=dirname//filenametmp,status='replace')
            write(54,*)burn_in,n_w,burn_in_widening,nsample_widening
            do i=1,(burn_in_widening+nsample_widening)
                write(54,*)convvp(i)
            enddo
            close(54)

            write(filenametmp,"('/Convergence_vs_',f5.2'_',I3.3,'.out')")widening,rank
            open(54,file=dirname//filenametmp,status='replace')
            write(54,*)burn_in,n_w,burn_in_widening,nsample_widening
            do i=1,(burn_in_widening+nsample_widening)
                write(54,*)convvs(i)
            enddo
            close(54)

            write(filenametmp,"('/Convergence_dp_',f5.2'_',I3.3,'.out')")widening,rank
            open(54,file=dirname//filenametmp,status='replace')
            write(54,*)burn_in,n_w,burn_in_widening,nsample_widening
            do i=1,(burn_in_widening+nsample_widening)
                write(54,*)convdp(i)
            enddo
            close(54)

            write(filenametmp,"('/Convergence_xi_',f5.2'_',I3.3,'.out')")widening,rank
            open(54,file=dirname//filenametmp,status='replace')
            write(54,*)burn_in,n_w,burn_in_widening,nsample_widening
            do i=1,(burn_in_widening+nsample_widening)
                write(54,*)convxi(i)
            enddo
            close(54)

            write(filenametmp,"('/Convergence_misfit_',f5.2'_',I3.3,'.out')")widening,rank
            open(54,file=dirname//filenametmp,status='replace')
            write(54,*)burn_in,n_w,burn_in_widening,nsample_widening
            do i=1,(burn_in_widening+nsample_widening)
                write(54,*)convd_R(i),convd_L(i)
            enddo
            close(54)

            write(filenametmp,"('/Convergence_nb_',f5.2'_',I3.3,'.out')")widening,rank
            open(53,file=dirname//filenametmp,status='replace')
            write(53,*)burn_in,n_w,burn_in_widening,nsample_widening
            do i=1,(burn_in_widening+nsample_widening)
                write(53,*)ncell(i)
            enddo
            close(53)

            write(filenametmp,"('/Convergence_sigma_R_',f5.2'_',I3.3,'.out')")widening,rank
            open(52,file=dirname//filenametmp,status='replace')
            write(52,*)burn_in,n_w,burn_in_widening,nsample_widening
            do i=1,(burn_in_widening+nsample_widening)
                write(52,*)convAd_R(i)
            enddo
            close(52)

            write(filenametmp,"('/Convergence_sigma_L_',f5.2'_',I3.3,'.out')")widening,rank
            open(53,file=dirname//filenametmp,status='replace')
            write(53,*)burn_in,n_w,burn_in_widening,nsample_widening
            do i=1,(burn_in_widening+nsample_widening)
                write(53,*)convAd_L(i)
            enddo
            close(53)

            write(filenametmp,"('/Convergence_sigmav_',f5.2'_',I3.3,'.out')")widening,rank
            open(53,file=dirname//filenametmp,status='replace')
            write(53,*)burn_in,n_w,burn_in_widening,nsample_widening
            do i=1,(burn_in_widening+nsample_widening)
                write(53,*)convsigmav(i)
            enddo
            close(53)

            write(filenametmp,"('/Convergence_Birth_',f5.2'_ref.out')")widening
            open(54,file=dirname//filenametmp,status='replace')
            write(54,*)burn_in,n_w,burn_in_widening,nsample_widening
            do i=1,(burn_in_widening+nsample_widening)
                write(54,*)convBs(i),convBas(i)
            enddo
            close(54)

            write(filenametmp,"('/Convergence_Death_',f5.2'_ref.out')")widening
            open(54,file=dirname//filenametmp,status='replace')
            write(54,*)burn_in,n_w,burn_in_widening,nsample_widening
            do i=1,(burn_in_widening+nsample_widening)
                write(54,*)convDs(i),convDas(i)
            enddo
            close(54)

            write(filenametmp,"('/Convergence_vp_',f5.2'_ref.out')")widening
            open(54,file=dirname//filenametmp,status='replace')
            write(54,*)burn_in,n_w,burn_in_widening,nsample_widening
            do i=1,(burn_in_widening+nsample_widening)
                write(54,*)convvps(i)
            enddo
            close(54)

            write(filenametmp,"('/Convergence_vs_',f5.2'_ref.out')")widening
            open(54,file=dirname//filenametmp,status='replace')
            write(54,*)burn_in,n_w,burn_in_widening,nsample_widening
            do i=1,(burn_in_widening+nsample_widening)
                write(54,*)convvss(i)
            enddo
            close(54)

            write(filenametmp,"('/Convergence_dp_',f5.2'_ref.out')")widening
            open(54,file=dirname//filenametmp,status='replace')
            write(54,*)burn_in,n_w,burn_in_widening,nsample_widening
            do i=1,(burn_in_widening+nsample_widening)
                write(54,*)convdps(i)
            enddo
            close(54)

            write(filenametmp,"('/Convergence_xi_',f5.2'_ref.out')")widening
            open(54,file=dirname//filenametmp,status='replace')
            write(54,*)burn_in,n_w,burn_in_widening,nsample_widening
            do i=1,(burn_in_widening+nsample_widening)
                write(54,*)convxis(i)
            enddo
            close(54)

            write(filenametmp,"('/Convergence_misfit_',f5.2'_ref.out')")widening
            open(54,file=dirname//filenametmp,status='replace')
            write(54,*)burn_in,n_w,burn_in_widening,nsample_widening
            do i=1,(burn_in_widening+nsample_widening)
                write(54,*)convd_Rs(i),convd_Ls(i)
            enddo
            close(54)

            write(filenametmp,"('/Convergence_nb_',f5.2'_ref.out')")widening
            open(53,file=dirname//filenametmp,status='replace')
            write(53,*)burn_in,n_w,burn_in_widening,nsample_widening
            do i=1,(burn_in_widening+nsample_widening)
                write(53,*)ncells(i)
            enddo
            close(53)

            write(filenametmp,"('/Convergence_sigma_R_',f5.2'_ref.out')")widening
            open(52,file=dirname//filenametmp,status='replace')
            write(52,*)burn_in,n_w,burn_in_widening,nsample_widening
            do i=1,(burn_in_widening+nsample_widening)
                write(52,*)convAd_Rs(i)
            enddo
            close(52)

            write(filenametmp,"('/Convergence_sigma_L_',f5.2'_ref.out')")widening
            open(53,file=dirname//filenametmp,status='replace')
            write(53,*)burn_in,n_w,burn_in_widening,nsample_widening
            do i=1,(burn_in_widening+nsample_widening)
                write(53,*)convAd_Ls(i)
            enddo
            close(53)

            write(filenametmp,"('/Convergence_sigmav_',f5.2'_ref.out')")widening
            open(53,file=dirname//filenametmp,status='replace')
            write(53,*)burn_in,n_w,burn_in_widening,nsample_widening
            do i=1,(burn_in_widening+nsample_widening)
                write(53,*)convsigmavs(i)
            enddo
            close(53)

        else
            write(filenametmp,"('/Convergence_Birth_',f5.2'_',I3.3,'.out')")widening,rank
            open(54,file=dirname//filenametmp,status='replace')
            write(54,*)burn_in,n_w,burn_in_widening,nsample_widening
            do i=1,(burn_in_widening+nsample_widening)
                write(54,*)convB(i),convBa(i)
            enddo
            close(54)

            write(filenametmp,"('/Convergence_Death_',f5.2'_',I3.3,'.out')")widening,rank
            open(54,file=dirname//filenametmp,status='replace')
            write(54,*)burn_in,n_w,burn_in_widening,nsample_widening
            do i=1,(burn_in_widening+nsample_widening)
                write(54,*)convD(i),convDa(i)
            enddo
            close(54)

            write(filenametmp,"('/Convergence_vp_',f5.2'_',I3.3,'.out')")widening,rank
            open(54,file=dirname//filenametmp,status='replace')
            write(54,*)burn_in,n_w,burn_in_widening,nsample_widening
            do i=1,(burn_in_widening+nsample_widening)
                write(54,*)convvp(i)
            enddo
            close(54)

            write(filenametmp,"('/Convergence_vs_',f5.2'_',I3.3,'.out')")widening,rank
            open(54,file=dirname//filenametmp,status='replace')
            write(54,*)burn_in,n_w,burn_in_widening,nsample_widening
            do i=1,(burn_in_widening+nsample_widening)
                write(54,*)convvs(i)
            enddo
            close(54)

            write(filenametmp,"('/Convergence_dp_',f5.2'_',I3.3,'.out')")widening,rank
            open(54,file=dirname//filenametmp,status='replace')
            write(54,*)burn_in,n_w,burn_in_widening,nsample_widening
            do i=1,(burn_in_widening+nsample_widening)
                write(54,*)convdp(i)
            enddo
            close(54)

            write(filenametmp,"('/Convergence_xi_',f5.2'_',I3.3,'.out')")widening,rank
            open(54,file=dirname//filenametmp,status='replace')
            write(54,*)burn_in,n_w,burn_in_widening,nsample_widening
            do i=1,(burn_in_widening+nsample_widening)
                write(54,*)convxi(i)
            enddo
            close(54)

            write(filenametmp,"('/Convergence_misfit_',f5.2'_',I3.3,'.out')")widening,rank
            open(54,file=dirname//filenametmp,status='replace')
            write(54,*)burn_in,n_w,burn_in_widening,nsample_widening
            do i=1,(burn_in_widening+nsample_widening)
                write(54,*)convd_R(i),convd_L(i)
            enddo
            close(54)

            write(filenametmp,"('/Convergence_nb_',f5.2'_',I3.3,'.out')")widening,rank
            open(53,file=dirname//filenametmp,status='replace')
            write(53,*)burn_in,n_w,burn_in_widening,nsample_widening
            do i=1,(burn_in_widening+nsample_widening)
                write(53,*)ncell(i)
            enddo
            close(53)

            write(filenametmp,"('/Convergence_sigma_R_',f5.2'_',I3.3,'.out')")widening,rank
            open(52,file=dirname//filenametmp,status='replace')
            write(52,*)burn_in,n_w,burn_in_widening,nsample_widening
            do i=1,(burn_in_widening+nsample_widening)
                write(52,*)convAd_R(i)
            enddo
            close(52)

            write(filenametmp,"('/Convergence_sigma_L_',f5.2'_',I3.3,'.out')")widening,rank
            open(53,file=dirname//filenametmp,status='replace')
            write(53,*)burn_in,n_w,burn_in_widening,nsample_widening
            do i=1,(burn_in_widening+nsample_widening)
                write(53,*)convAd_L(i)
            enddo
            close(53)

            write(filenametmp,"('/Convergence_sigmav_',f5.2'_',I3.3,'.out')")widening,rank
            open(53,file=dirname//filenametmp,status='replace')
            write(53,*)burn_in,n_w,burn_in_widening,nsample_widening
            do i=1,(burn_in_widening+nsample_widening)
                write(53,*)convsigmav(i)
            enddo
            close(53)

        ENDIF

        call MPI_Group_free(good_group, ierror)
        call MPI_Comm_free(MPI_COMM_small, ierror)
        call MPI_BARRIER(MPI_COMM_WORLD, ierror)

        widening=widening+widening_step

    enddo
    CALL cpu_time(t2)
    if (ran==0) write(*,*)'time taken by the code was',t2-t1,'seconds'
    if (ran==0) write(*,*)'time taken by the code was',(t2-t1)/3600,'hours'


    call MPI_FINALIZE(ierror)! Terminate the parralelization


end! end the main program




!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



!----------------------------------------------------------------------

!               FUNCTIONS USED BY THE MAIN CODE

!---------------------------------------------------------------------


!-------------------------------------------------------------------
!
!    Numerical Recipes random number generator for
!       a Gaussian distribution
!
! ----------------------------------------------------------------------------



FUNCTION GASDEV(idum)

    !     ..Arguments..
    integer          idum
    real GASDEV

    !     ..Local..
    real v1,v2,r,fac
    real ran3

    if (idum.lt.0) iset=0
10  v1=2*ran3(idum)-1
    v2=2*ran3(idum)-1
    r=v1**2+v2**2
    if(r.ge.1.or.r.eq.0) GOTO 10
    fac=sqrt(-2*log(r)/r)
    GASDEV=v2*fac

    RETURN
END


!-------------------------------------------------------------------
!
!    Numerical Recipes random number generator
!
! ----------------------------------------------------------------------------

FUNCTION ran3(idum)
    INTEGER idum
    INTEGER MBIG,MSEED,MZ
    !     REAL MBIG,MSEED,MZ
    REAL ran3,FAC
    PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1./MBIG)
    !     PARAMETER (MBIG=4000000.,MSEED=1618033.,MZ=0.,FAC=1./MBIG)
    INTEGER i,iff,ii,inext,inextp,k
    INTEGER mj,mk,ma(55)
    !     REAL mj,mk,ma(55)
    SAVE iff,inext,inextp,ma
    DATA iff /0/
    ! write(*,*)' idum ',idum
    if(idum.lt.0.or.iff.eq.0)then
        iff=1
        mj=MSEED-iabs(idum)
        mj=mod(mj,MBIG)
        ma(55)=mj
        mk=1
        do 11 i=1,54
            ii=mod(21*i,55)
            ma(ii)=mk
            mk=mj-mk
            if(mk.lt.MZ)mk=mk+MBIG
            mj=ma(ii)
        !  write(*,*)' idum av',idum
11      continue
        do 13 k=1,4
            do 12 i=1,55
                ma(i)=ma(i)-ma(1+mod(i+30,55))
                if(ma(i).lt.MZ)ma(i)=ma(i)+MBIG
12          continue
13      continue
        ! write(*,*)' idum ap',idum
        inext=0
        inextp=31
        idum=1
    endif
    ! write(*,*)' idum app ',idum
    inext=inext+1
    if(inext.eq.56)inext=1
    inextp=inextp+1
    if(inextp.eq.56)inextp=1
    mj=ma(inext)-ma(inextp)
    if(mj.lt.MZ)mj=mj+MBIG
    ma(inext)=mj
    ran3=mj*FAC
    !  write(*,*)' idum ',idum

    return
END

!-------------------------------------------------------------------
!
!    create new mpi operation for alpha
!
! ----------------------------------------------------------------------------



!the code is finished.
