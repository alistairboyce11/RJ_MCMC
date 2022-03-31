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
    real th
    integer nlims_cur,nlims_cur_diff
    logical ier, error_flag ! error flags for dispersion curves and minos_bran
    integer npt,npt_prop,npt_iso,npt_ani ! numbers of points, isotropic, anisotropic
    logical accept,tes,birth,birtha,death,deatha,move,value,noisd_R,noisd_L,ani,change_vp !which parameter did we change?
    logical testing
    integer ind2,ind_vp,ind_vsv,ind_xi,j !more indices, mostly for posterior. j sometimes is the number of processors, careful at the end!
    real d !for calculating posterior, mostly depth
    real histo(maxlay),histos(maxlay),histoch(disd),histochs(disd) ! histo of layers, changepoints 
    real avvs(disd),avvp(disd),avxi(disd),avxis(disd),avvss(disd),avvps(disd) ! average of posterior
    real peri_R(ndatadmax),peri_L(ndatadmax) !periods of data
    real peri_R_tmp(ndatadmax),peri_L_tmp(ndatadmax) !periods of data
    real postvp(disd,disv),postvs(disd,disv),postxi(disd,disv),postvps(disd,disv),postvss(disd,disv),postxis(disd,disv) ! posteriors
    real probani(disd),probanis(disd) !anisotropy probabilities
    integer n_R(ndatadmax),n_L(ndatadmax),nmax_R,nmax_L,nmin_R,nmin_L !harmonic number of modes
    integer n_R_tmp(ndatadmax),n_L_tmp(ndatadmax)
    integer nlims_R(2,nharmo_max),nlims_L(2,nharmo_max),numharm_count,numharm_R(nharmo_max),numharm_L(nharmo_max)

    real PrB,AcB,PrD,AcD,Prnd_R,Acnd_R,Prnd_L,Acnd_L,&
        Acba,Prba,Acda,Prda,out,Prxi,Acxi,Pr_vp,&
        Ac_vp,PrP,PrV,AcP,AcV ! to adjust acceptance rates for all different parameters
    real lsd_L,lsd_L_prop,lsd_L_min,lsd_R,lsd_R_prop,lsd_R_min !logprob without uncertainties
    real logprob_vsv,logprob_vp !for reversible jump
    real convAd_R(nsample+burn_in),convAd_Rs(nsample+burn_in),convAd_L(nsample+burn_in),convAd_Ls(nsample+burn_in) !variations of uncertainty parameters
    real voro(malay,4),vsref_min,vsref_max,vpref_min,vpref_max !the model with its bounds
    real voro_prop(malay,4),ncell(nsample+burn_in),convd_R(nsample+burn_in),convd_L(nsample+burn_in),&
        convd_Rs(nsample+burn_in),convd_Ls(nsample+burn_in) ! the proposed model, history of lsd
    real ncells(nsample+burn_in),t1,t2 !timers, history of ls
    real like,like_prop,u,& !log-likelyhoods
        liked_R_prop,liked_R,liked_L_prop,liked_L
    real ML_Ad_R(disA),ML_Ad_Rs(disA),TRA(malay,malay+1),& !histogram of uncertainty parameters
        TRAs(malay,malay+1),ML_Ad_L(disA),ML_Ad_Ls(disA)
    double precision logrsig,Ad_R,Ad_R_prop,Ad_L,Ad_L_prop !uncertainty parameters
    real sigmav,sigmav_old,sigmav_new,AR_birth_old   !proposal on velocity when Birth move  
    real d_obsdcR(ndatadmax),d_obsdcL(ndatadmax),d_obsdcRe(ndatadmax),d_obsdcLe(ndatadmax) !observed data
    integer inorout(21000),inorouts(21000),members(21000),io ! mpi related variables
    integer ndatad_R,ndatad_L,ndatad_R_tmp,ndatad_L_tmp !number of observed data points
    real pxi,p_vp,pd,pv,pAd_R,pAd_L ! related to observed data for azimuthal anisotropy
    ! Geometry parameters
    ! Traces
    integer nptref,malayv ! number of points in reference model, number of layers in voro
    real convBs(nsample+burn_in),convB(nsample+burn_in),convPs(nsample+burn_in),convP(nsample+burn_in) ! birth rates, vp change rates
    real convvs(nsample+burn_in),convvss(nsample+burn_in)
    real convdp(nsample+burn_in),convdps(nsample+burn_in)
    real convvp(nsample+burn_in),convvps(nsample+burn_in),convxis(nsample+burn_in),convxi(nsample+burn_in)
    real convBas(nsample+burn_in),convBa(nsample+burn_in) ! anisotropic birth rates
    real convDs(nsample+burn_in),convD(nsample+burn_in),convDas(nsample+burn_in),convDa(nsample+burn_in)

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
    character filename*13, number*4
    real likemax ! log-likelihood of best model
    integer nptmax ! number of layers of the best model
    real voromax(malay,4) ! best model
    character filenamemax*30  !filename for writing the best model
    real d_cRmoy(ndatadmax), d_cRmoys(ndatadmax), d_cLmoy(ndatadmax), d_cLmoys(ndatadmax) ! average dispersion curve
    real d_cRdelta(ndatadmax), d_cRdeltas(ndatadmax), d_cLdelta(ndatadmax), d_cLdeltas(ndatadmax) ! standart deviation of dispersion curve
    real model_ref(mk,9) ! reference model, all models are deviations from it
    real,dimension(mk) :: r,rho,vpv,vph,vsv,vsh,qkappa,qshear,eta,xi,vpvsv_data ! input for forward modelling
    integer nptfinal,nic,noc,nic_ref,noc_ref,jcom ! more inputs
    integer nmodes,n_mode(nmodes_max),l_mode(nmodes_max) ! outputs of forward modelling
    real c_ph(nmodes_max),period(nmodes_max),raylquo(nmodes_max),tref ! more outputs, rayquo: error measurement, should be of the order of eps
    real wmin_R(nharmo_max),wmax_R(nharmo_max),wmin_L(nharmo_max),wmax_L(nharmo_max) ! more inputs f

    integer recalculated !counts the number of times we need to improve eps
    logical stuck
    integer nharm_R,nharm_L,iharm
    integer idis, numdis,numdis_tmp,i_al ! additional stuff for joint inversion
    real like_alt,like_alt_L,like_alt_R
    real logcratio(numdis_max),logalpharef,logalpha(numdis_max),alphasum(numdis_max),alpharefsum,logcratioref
    real alpha(numdis_max),alpharef
    real postvp_alt(disd,disv,numdis_max),postvs_alt(disd,disv,numdis_max),postxi_alt(disd,disv,numdis_max)
    real postvps_alt(disd,disv,numdis_max),postvss_alt(disd,disv,numdis_max),postxis_alt(disd,disv,numdis_max) 
    
    real avvs_alt(disd,numdis_max),avvp_alt(disd,numdis_max),avxi_alt(disd,numdis_max)
    real avxis_alt(disd,numdis_max),avvss_alt(disd,numdis_max),avvps_alt(disd,numdis_max)
    real d_obsdcR_alt(ndatadmax,numdis_max),d_obsdcL_alt(ndatadmax,numdis_max)
    real d_obsdcRe_alt(ndatadmax,numdis_max),d_obsdcLe_alt(ndatadmax,numdis_max) !observed data
    real d_cRmoy_alt(ndatadmax,numdis_max),d_cLmoy_alt(ndatadmax,numdis_max)
    real d_cRdelta_alt(ndatadmax,numdis_max),d_cLdelta_alt(ndatadmax,numdis_max) 
    real d_cRmoy_alts(ndatadmax,numdis_max),d_cLmoy_alts(ndatadmax,numdis_max)
    real d_cRdelta_alts(ndatadmax,numdis_max),d_cLdelta_alts(ndatadmax,numdis_max) 
    real histo_alt(maxlay,numdis_max),probani_alt(disd,numdis_max)
    real histo_alts(maxlay,numdis_max),probanis_alt(disd,numdis_max)
    real logalhists(num_logalpha,numdis_max),logalhist(num_logalpha,numdis_max)
    real logalpharefhist(num_logalpha),logalpharefhists(num_logalpha)
    real histoch_alt(disd,numdis_max),histoch_alts(disd,numdis_max)
    real ML_Ad_R_alt(disA,numdis_max),ML_Ad_R_alts(disA,numdis_max),& !histogram of uncertainty parameters
        ML_Ad_L_alt(disA,numdis_max),ML_Ad_L_alts(disA,numdis_max)
    real like_w,like_prop_w,TRA_alt(malay,malay+1,numdis_max),TRA_alts(malay,malay+1,numdis_max)
    
    ! widened results
    real th_wide
    real histo_wide(maxlay),probani_wide(disd)
    real histo_wides(maxlay),probanis_wide(disd)
    real histoch_wide(disd),histoch_wides(disd)
    real ML_Ad_R_wide(disA),ML_Ad_R_wides(disA),& !histogram of uncertainty parameters
        ML_Ad_L_wide(disA),ML_Ad_L_wides(disA)
    real avvs_wide(disd),avvp_wide(disd),avxi_wide(disd)
    real avxis_wide(disd),avvss_wide(disd),avvps_wide(disd)
    real postvp_wide(disd,disv),postvs_wide(disd,disv),postxi_wide(disd,disv)
    real postvp_wides(disd,disv),postvs_wides(disd,disv),postxi_wides(disd,disv) ! posteriors
    real d_cRmoy_wide(ndatadmax),d_cLmoy_wide(ndatadmax)
    real d_cRdelta_wide(ndatadmax),d_cLdelta_wide(ndatadmax) 
    real d_cRmoy_wides(ndatadmax),d_cLmoy_wides(ndatadmax)
    real d_cRdelta_wides(ndatadmax),d_cLdelta_wides(ndatadmax) 
    real TRA_wide(malay,malay+1),TRA_wides(malay,malay+1)
    
    integer i_w
    real widening_prop,widening
    real mean_best(numdis_max),mean_prop(numdis_max),mean_props(numdis_max),meandiff_hist(n_w,numdis_max)
    real alphahist(num_logalpha,numdis_max,n_w),alphahists(num_logalpha,numdis_max,n_w)
    real alphamax(numdis_max),alpharefmax
    real alphamax_prop(numdis_max),alpharefmax_prop
    real alphamax_props(numdis_max),alpharefmax_props
    real alphaall(nsample/thin,numdis_max)
    real, dimension(:), allocatable :: lat, lon
    
    character filebycore*15

1000 format(I4)


    !***********************************************************************

    CALL cpu_time(t1)  !Tic. start counting time 
    
    call MPI_INIT(ierror)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nbproc, ierror) ! MPI_COMM_WORLD: communicator (https://www.open-mpi.org/doc/v4.0/man3/MPI_Comm_size.3.php)
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror) ! https://www.open-mpi.org/doc/v4.0/man3/MPI_Comm_rank.3.php
    call MPI_Comm_group(MPI_COMM_WORLD,group_world,ierror) ! https://www.open-mpi.org/doc/v3.0/man3/MPI_Comm_group.3.php
    
    !Start Parralelization of the code. From now on, the code is run on each
    !processor independently, with ra = the number of the proc.

    !-!
    
    
    getting_old=.true.
    !testing=.true.
    !if (testing) write(*,*)'testing with synthetic model'
    
    ra=0
    
    !************************************************************
    !                READ PREM 
    !************************************************************
    open(7,file="Model_deep_REF.in",status='old',form='formatted')
    !  110 format(20a4)
    read(7,*) nptref,nic_ref,noc_ref
    read(7,'(f8.0, 3f9.2, 2f9.1, 2f9.2, f9.5)') (model_ref(i,1),&
        model_ref(i,2),model_ref(i,3),model_ref(i,4),model_ref(i,5),&
        model_ref(i,6),model_ref(i,7),model_ref(i,8),model_ref(i,9),&
        i=1,nptref)
    close(7)
    j=1
    do while (model_ref(j,1)<rearth-d_max*1000.)
        j=j+1
    end do
    j=j-1
    
    vsref_min=minval(model_ref(j:nptref,4)*(1-width)) ! setting min/max velocities for writing later
    vsref_max=maxval(model_ref(j:nptref,4)*(1+width))
    
    vpref_min=minval(model_ref(j:nptref,4)*vpvs*(1+vpvsv_min/100.)*(1-width)) 
    vpref_max=maxval(model_ref(j:nptref,4)*vpvs*(1+vpvsv_max/100.)*(1+width))
    
    ! GET SWD DATA ---------------------------------------------------------------- 
!     nlims_cur=1
!     open(65,file=dirname//'/dispersion.in',status='old')! 65: name of the opened file in memory (unit identifier)
!     read(65,*,IOSTAT=io)numdis
!     read(65,*,IOSTAT=io)lat(:numdis)
!     read(65,*,IOSTAT=io)lon(:numdis)
!     read(65,*,IOSTAT=io)ndatad_R
!     read(65,*,IOSTAT=io)nharm_R
!     do k=1,nharm_R
!         read(65,*,IOSTAT=io)numharm_R(k)
!         read(65,*,IOSTAT=io)nlims_cur_diff
!         nlims_R(1,k)=nlims_cur
!         nlims_R(2,k)=nlims_R(1,k)+nlims_cur_diff
!         !read(65,*,IOSTAT=io)wmin_R(k),wmax_R(k)
!         do i=nlims_cur,nlims_cur+nlims_cur_diff
!             read(65,*,IOSTAT=io)n_R(i),peri_R(i),d_obsdcR(i),d_obsdCRe(i)
!             do j=1,numdis
!                 read(65,*,IOSTAT=io)d_obsdcR_alt(i,j),d_obsdCRe_alt(i,j)
!             enddo
!         enddo
!         nlims_cur=nlims_R(2,k)+1
!     enddo
!     
!     
!     nlims_cur=1
!     read(65,*,IOSTAT=io)ndatad_L
!     read(65,*,IOSTAT=io)nharm_L
!     do k=1,nharm_L
!         read(65,*,IOSTAT=io)numharm_L(k)
!         read(65,*,IOSTAT=io)nlims_cur_diff
!         nlims_L(1,k)=nlims_cur
!         nlims_L(2,k)=nlims_L(1,k)+nlims_cur_diff
!         !read(65,*,IOSTAT=io)wmin_L(k),wmax_L(k)
!         do i=nlims_cur,nlims_cur+nlims_cur_diff
!             read(65,*,IOSTAT=io)n_L(i),peri_L(i),d_obsdcL(i),d_obsdCLe(i)
!             do j=1,numdis
!                 read(65,*,IOSTAT=io)d_obsdcL_alt(i,j),d_obsdCLe_alt(i,j)
!             enddo
!         enddo
!         nlims_cur=nlims_L(2,k)+1
!     enddo
!     close(65)
!     
!     do k=1,nharm_R
!         wmin_R(k)=1000./(maxval(peri_R(nlims_R(1,k):nlims_R(2,k)))+20)
!         wmax_R(k)=1000./(minval(peri_R(nlims_R(1,k):nlims_R(2,k)))-2)
!     enddo
!     
!     do k=1,nharm_L
!         wmin_L(k)=1000./(maxval(peri_L(nlims_L(1,k):nlims_L(2,k)))+20)
!         wmax_L(k)=1000./(minval(peri_L(nlims_L(1,k):nlims_L(2,k)))-2)
!     enddo
!     
!     nmin_R=minval(n_R(:ndatad_R))
!     nmax_R=maxval(n_R(:ndatad_R))
!     
!     nmin_L=minval(n_L(:ndatad_L))
!     nmax_L=maxval(n_L(:ndatad_L))
!     
!     if (ndatad_R>=ndatad_L) tref=sum(peri_R(:ndatad_R))/ndatad_R
!     if (ndatad_L>ndatad_R) tref=sum(peri_L(:ndatad_L))/ndatad_L
    nlims_cur=1
    open(65,file=dirname//'/dispersion.in',status='old')! 65: name of the opened file in memory (unit identifier)
    read(65,*,IOSTAT=io)numdis
    allocate(lat(numdis))
    read(65,*,IOSTAT=io)lat
    allocate(lon(numdis))
    read(65,*,IOSTAT=io)lon
    read(65,*,IOSTAT=io)ndatad_R
    read(65,*,IOSTAT=io)nharm_R
    do k=1,nharm_R
        read(65,*,IOSTAT=io)numharm_R(k)
        read(65,*,IOSTAT=io)nlims_cur_diff
        nlims_R(1,k)=nlims_cur
        nlims_R(2,k)=nlims_R(1,k)+nlims_cur_diff
        !read(65,*,IOSTAT=io)wmin_R(k),wmax_R(k)
        do i=nlims_cur,nlims_cur+nlims_cur_diff
            read(65,*,IOSTAT=io)n_R(i),peri_R(i),d_obsdcR(i),d_obsdCRe(i)
            do j=1,numdis
                read(65,*,IOSTAT=io)d_obsdcR_alt(i,j),d_obsdCRe_alt(i,j)
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
        !read(65,*,IOSTAT=io)wmin_L(k),wmax_L(k)
        do i=nlims_cur,nlims_cur+nlims_cur_diff
            read(65,*,IOSTAT=io)n_L(i),peri_L(i),d_obsdcL(i),d_obsdCLe(i)
            do j=1,numdis
                read(65,*,IOSTAT=io)d_obsdcL_alt(i,j),d_obsdCLe_alt(i,j)
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
    
    ! get correct widening from previous test
    open(51,file=dirname//'/goodalpha.inout',status='old')
    read(51,*,IOSTAT=io)widening,alpharefmax
    read(51,*,IOSTAT=io)numdis_tmp
    if (numdis/=numdis_tmp) then
        stop 'numdis from dispersion different from numdis from widening test'
    endif
    do i=1,numdis
        read(51,*,IOSTAT=io)alphamax(i)
    enddo
    close(51)

!**************************************************************
    
    !    True loop
    
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
    postxi=0
    postxis=0
    postvp=0
    postvps=0
    postvs=0
    postvss=0
    avvs=0
    avvp=0
    avxi=0
    avvss=0
    avvps=0
    avxis=0
 
    lsd_R=0
    lsd_L=0
    histoch=0
    histochs=0
    histoch_alt=0
    histoch_alts=0
    histoch_wide=0
    histoch_wides=0
    ML_Ad_R=0
    ML_Ad_Rs=0
    ML_Ad_L=0
    ML_Ad_Ls=0
    ML_Ad_R_alt=0
    ML_Ad_R_alts=0
    ML_Ad_L_alt=0
    ML_Ad_L_alts=0
    ML_Ad_R_wide=0
    ML_Ad_R_wides=0
    ML_Ad_L_wide=0
    ML_Ad_L_wides=0
    probani=0
    probanis=0
    
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
    histo=0
    histos=0
    convAd_R=0
    convAd_Rs=0
    convAd_L=0
    convAd_Ls=0
    TRA=0
    TRAs=0
    TRA_alt=0
    TRA_alts=0
    TRA_wide=0
    TRA_wides=0
    d_cRmoy=0
    d_cLmoy=0
    d_cRmoys=0
    d_cLmoys=0
    d_cRdelta=0
    d_cLdelta=0
    d_cRdeltas=0
    d_cLdeltas=0 
    d_cRmoy_alt=0
    d_cLmoy_alt=0
    d_cRmoy_alts=0
    d_cLmoy_alts=0
    d_cRdelta_alt=0
    d_cLdelta_alt=0
    d_cRdelta_alts=0
    d_cLdelta_alts=0 
    d_cRmoy_wide=0
    d_cLmoy_wide=0
    d_cRmoy_wides=0
    d_cLmoy_wides=0
    d_cRdelta_wide=0
    d_cLdelta_wide=0
    d_cRdelta_wides=0
    d_cLdelta_wides=0 
    ier=.false.
    
    avvs_alt=0
    avvp_alt=0
    avxi_alt=0
    avxis_alt=0
    avvss_alt=0
    avvps_alt=0
    histo_alt=0
    histo_alts=0
    probani_alt=0
    probanis_alt=0
    logalhists=0
    logalhist=0
    alphasum=0
    avvs_wide=0
    avvp_wide=0
    avxi_wide=0
    avxis_wide=0
    avvss_wide=0
    avvps_wide=0
    histo_wide=0
    histo_wides=0
    probani_wide=0
    probanis_wide=0
    
    !**************************************************************

    !    Draw the first model randomly from the prior distribution

    !**************************************************************
    
    
    ! Initilalise sigma (noise parameter)
    Ad_R = Ad_R_max-pAd_R! Ad_min+ran3(ra)*(Ad_R_max-Ad_min)
    Ad_L = Ad_L_max-pAd_L

    ! Initial number of cells
    !------------------------------------


    j=0
    tes=.false.
    do while(.not.tes)
50      tes=.true.
        
        if (getting_old) then
            write(filenamemax,"('/last_model_',I3.3,'.inout')") rank    
            open(65,file=dirname//filenamemax,status='old')
            read(65,*,IOSTAT=io)Ad_R
            read(65,*,IOSTAT=io)Ad_L
            read(65,*,IOSTAT=io)pxi
            read(65,*,IOSTAT=io)p_vp
            read(65,*,IOSTAT=io)pAd_R
            read(65,*,IOSTAT=io)pAd_L
            read(65,*,IOSTAT=io)pv
            read(65,*,IOSTAT=io)pd
            read(65,*,IOSTAT=io)sigmav
            read(65,*,IOSTAT=io)npt
            do i=1,npt
                read(65,*,IOSTAT=io)voro(i,1),voro(i,2),voro(i,3),voro(i,4)
            enddo
            close(65)
            j=j+1
        else
        
            ! create a starting model randomly
            npt = milay+ran3(ra)*(malay-milay)!(12-milay)!(maxlay-minlay)
            !npt=11
            if ((npt>malay).or.(npt<milay)) goto 50
            !write(*,*)npt
            j=j+1
            do i=1,npt
    
                voro(i,1)= d_min+ran3(ra)*(d_max-d_min) 
                voro(i,2)= (-width+2*width*ran3(ra))
                if (ran3(ra)<0.5) then
                    voro(i,3)= xi_min+(xi_max-xi_min)*ran3(ra)
                else
                    voro(i,3)= -1
                endif
                !voro(i,3)=-1
                voro(i,4)= vpvsv_min+(vpvsv_max-vpvsv_min)*ran3(ra)
            enddo
            
    !         open(65,file=dirname//'/stuck_039.out',status='old')
    !         read(65,*,IOSTAT=io)npt
    !         do i=1,npt
    !             read(65,*,IOSTAT=io)voro(i,1),voro(i,2),voro(i,3),voro(i,4)
    !         enddo
    !         close(65)

        endif

        call combine(model_ref,nptref,nic_ref,noc_ref,voro,npt,d_max,&
            r,rho,vpv,vph,vsv,vsh,qkappa,qshear,eta,nptfinal,nic,noc,xi,vpvsv_data)
        write(*,*)"Combined"
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
                else
                    write(*,*)"Minos_bran SUCCESS for RAYLEIGH"
                end if
                peri_R_tmp=peri_R(nlims_R(1,iharm):nlims_R(2,iharm))
                n_R_tmp=n_R(nlims_R(1,iharm):nlims_R(2,iharm))
                ndatad_R_tmp=nlims_R(2,iharm)-nlims_R(1,iharm)+1 ! fortran slices take the first and the last element
                call dispersion_minos(nmodes_max,nmodes,n_mode,c_ph,period,raylquo,&
                    peri_R_tmp,n_R_tmp,d_cR_tmp,rq_R,ndatad_R_tmp,ier) ! extract phase velocities from minos output (pretty ugly)
                write(*,*)"dispersion_minos for RAYLEIGH"
                d_cR(nlims_R(1,iharm):nlims_R(2,iharm))=d_cR_tmp
                if (ier) tes=.false.
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
                else
                    write(*,*)"Minos_bran SUCCESS for LOVE"
                end if
                peri_L_tmp=peri_L(nlims_L(1,iharm):nlims_L(2,iharm))
                n_L_tmp=n_L(nlims_L(1,iharm):nlims_L(2,iharm))
                ndatad_L_tmp=nlims_L(2,iharm)-nlims_L(1,iharm)+1 ! fortran slices take the first and the last element
                call dispersion_minos(nmodes_max,nmodes,n_mode,c_ph,period,raylquo,&
                peri_L_tmp,n_L_tmp,d_cL_tmp,rq_L,ndatad_L_tmp,ier) ! extract phase velocities from minos output (pretty ugly)
                write(*,*)"dispersion_minos for LOVE"
                d_cL(nlims_L(1,iharm):nlims_L(2,iharm))=d_cL_tmp
                if (ier) tes=.false.
                if (maxval(abs(rq_L(:ndatad_L_tmp)))>maxrq*eps) tes=.false.
            enddo
        
        endif
        
        if (j>250) stop "CAN'T INITIALIZE MODEL" ! if it doesn't work after 250 tries, give up
        
        if (.not.tes) getting_old=.false.

    enddo

    write(*,*)'DONE INITIALIZING'
    
    ! isoflag says if a layer is isotropic
    do i=1,npt
        if (voro(i,3)==-1) then
            isoflag(i)=.true.
        else
            isoflag(i)=.false.
        endif
    enddo
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Check and double-check
    write(*,*)'getting the ratio of the uncertainties'
    
    write(*,*)'getting initial likelihood'
    
    
    !***********************************************************

    !                 Get initial likelihood

    !***********************************************************

    lsd_R=0
    lsd_L=0
    liked_R=0
    liked_L=0
    
    do i=1,ndatad_R ! calculate misfit -> log-likelihood of initial model
        lsd_R=lsd_R+(d_obsdCR(i)-d_cR(i))**2
        liked_R=liked_R+(d_obsdCR(i)-d_cR(i))**2/(2*(Ad_R*d_obsdCRe(i))**2) ! gaussian errors
    enddo
    do i=1,ndatad_L
        lsd_L=lsd_L+(d_obsdCL(i)-d_cL(i))**2
        liked_L=liked_L+(d_obsdCL(i)-d_cL(i))**2/(2*(Ad_L*d_obsdCLe(i))**2) 
    enddo
    lsd_R_min=lsd_R
    lsd_L_min=lsd_L
    likemax=lsd_R+lsd_L
    
    like= (liked_R + liked_L)
    like_w=like/widening
    
    write(*,*)like,like_w
    
    stuck=.false.
    
    sample=0
    th=0
    if (getting_old) then
        ount=burn_in-burn_in_widening
    else
        ount=0
    endif
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
    
    recalculated=0
    
    write(*,*)'all initialized' 

    do while (sample<nsample) ! main loop, sample: number of sample post burn-in
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
            
!             if ((Ac_vp/(Pr_vp+1))<0.01) then 
!                 write(filenamemax,"('/stuck_',I3.3,'.out')") rank    
!                 open(56,file=dirname//filenamemax,status='replace')
!                 write(56,*) npt
!                 do i=1,npt
!                     write(56,*)voro(i,1),voro(i,2),voro(i,3),voro(i,4)
!                 enddo
!                 close(56)
!             endif
            
            !-----------------------------------------------

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
        logprob_vp=0
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
                if ((voro_prop(ind,4)<=vpvsv_min).or.(voro_prop(ind,4)>=vpvsv_max)) out=0
                
            endif     
        elseif (u<0.3) then !change position--------------------------------------------
            move=.true.
            ind=1+ceiling(ran3(ra)*(npt-1))
            if (ind==1) ind=2
            !write(*,*)ind
 
            !if (ount.GT.burn_in) then 
            PrP=PrP+1
            voro_prop(ind,1)=voro(ind,1)+gasdev(ra)*pd
        
            if ((voro_prop(ind,1)<=d_min).or.(voro_prop(ind,1)>=d_max)) then
                out=0
            endif
            if ((voro_prop(ind,2)<=-width).or.(voro_prop(ind,2)>=width)) then
                out=0
            endif
            !if (out==1) write(*,*)voro_prop(ind,1)
 
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
                voro_prop(npt_prop,4) = vpvsv_min+(vpvsv_max-vpvsv_min)*ran3(ra)
                isoflag_prop(npt_prop) = .true. 
                
                logprob_vsv=log(1/(sigmav*sqrt(2*pi)))-((voro(ind,2)-voro_prop(npt_prop,2))**2)/(2*sigmav**2) ! correct acceptance rates because transdimensional
                !logprob_vp=log(1/(p_vp*sqrt(2*pi)))-((voro(ind,4)-voro_prop(npt_prop,4))**2)/(2*p_vp**2)
                
                !Check bounds                    
                if ((voro_prop(npt_prop,2)<=-width).or.(voro_prop(npt_prop,2)>=width)) then
                    out=0
                end if
                if ((voro_prop(npt_prop,4)<=vpvsv_min).or.(voro_prop(npt_prop,4)>=vpvsv_max)) then
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
            
            if ((npt_prop<milay)) then!.or.(ind==1)) then
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

            call combine(model_ref,nptref,nic_ref,noc_ref,voro_prop,npt_prop,d_max,&
                r,rho,vpv,vph,vsv,vsh,qkappa,qshear,eta,nptfinal,nic,noc,xi,vpvsv_data)
            
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
                    
                    peri_R_tmp=peri_R(nlims_R(1,iharm):nlims_R(2,iharm))
                    n_R_tmp=n_R(nlims_R(1,iharm):nlims_R(2,iharm))
                    ndatad_R_tmp=nlims_R(2,iharm)-nlims_R(1,iharm)+1 ! fortran slices take the first and the last element
                    call dispersion_minos(nmodes_max,nmodes,n_mode,c_ph,period,raylquo,&
                        peri_R_tmp,n_R_tmp,d_cR_prop_tmp,rq_R,ndatad_R_tmp,ier) ! extract phase velocities from minos output (pretty ugly)
                    !write(*,*)"dispersion_minos for RAYLEIGH"
                    d_cR_prop(nlims_R(1,iharm):nlims_R(2,iharm))=d_cR_prop_tmp
                    
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
                    
                    peri_L_tmp=peri_L(nlims_L(1,iharm):nlims_L(2,iharm))
                    n_L_tmp=n_L(nlims_L(1,iharm):nlims_L(2,iharm))
                    ndatad_L_tmp=nlims_L(2,iharm)-nlims_L(1,iharm)+1 ! fortran slices take the first and the last element
                    call dispersion_minos(nmodes_max,nmodes,n_mode,c_ph,period,raylquo,&
                        peri_L_tmp,n_L_tmp,d_cL_prop_tmp,rq_L,ndatad_L_tmp,ier) ! extract phase velocities from minos output (pretty ugly)
                    !write(*,*)"dispersion_minos for LOVE"
                    d_cL_prop(nlims_L(1,iharm):nlims_L(2,iharm))=d_cL_prop_tmp
                    
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
                liked_R_prop=liked_R_prop+(d_obsdCR(i)-d_cR_prop(i))**2/(2*(Ad_R_prop*d_obsdCRe(i))**2) ! prendre en compte erreurs mesurées
            end do
            do i=1,ndatad_L
                lsd_L_prop=lsd_L_prop+(d_obsdCL(i)-d_cL_prop(i))**2
                liked_L_prop=liked_L_prop+(d_obsdCL(i)-d_cL_prop(i))**2/(2*(Ad_L_prop*d_obsdCLe(i))**2) ! prendre en compte erreurs mesurées
            end do
            
        endif 
        
        like_prop=(liked_R_prop+liked_L_prop) !log-likelihood of the proposed model
        like_prop_w=like_prop/widening
        
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
!         if (ind<3) then
!         write(*,*)accept,ind,move,value,noisd_R
!         endif
        
        if (accept) then 
            
!             if (ind<3) then
!             write(*,*)ind,move,value,noisd_R
!             endif
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
                    
            
            if (lsd_R<lsd_R_min)then
                lsd_R_min = lsd_R            
            endif
            if (lsd_L<lsd_L_min)then
                lsd_L_min = lsd_L            
            endif

            
            
            !**********************************************************************
                
            !       Save best model
            
            !**********************************************************************
            
            if ((lsd_L+lsd_R).lt.likemax) then
                        
                
                voromax=voro
                nptmax=npt
                likemax=(lsd_L+lsd_R)
            endif
            
            
            logalpharef=ndatad_R*(1/widening-1)*log(Ad_R)&
            +ndatad_L*(1/widening-1)*log(Ad_L)&
            +like/widening-like&
            -alpharefmax
            logalpharef=min(logalpharef,0.)
            !+logcratioref
            !-(1/widening-1)*(ndatad_R+ndatad_L)/2&
            !
            
            do idis=1,numdis
                like_alt_R=0
                like_alt_L=0
                do i=1,ndatad_R
                    like_alt_R=like_alt_R+(d_obsdCR_alt(i,idis)-d_cR(i))**2/(2*(Ad_R*d_obsdCRe_alt(i,idis))**2) ! prendre en compte erreurs mesurées
                end do
                do i=1,ndatad_L
                    like_alt_L=like_alt_L+(d_obsdCL_alt(i,idis)-d_cL(i))**2/(2*(Ad_L*d_obsdCLe_alt(i,idis))**2) ! prendre en compte erreurs mesurées
                end do
                like_alt=like_alt_R+like_alt_L
                
                
                
!                 do i=1,ndatad_R
!                     liked_R_prop=liked_R_prop+(d_obsdCR(i)-d_cR_prop(i))**2/(2*(Ad_R_prop*d_obsdCRe(i))**2) ! prendre en compte erreurs mesurées
!                 end do
!                 do i=1,ndatad_L
!                     liked_L_prop=liked_L_prop+(d_obsdCL(i)-d_cL_prop(i))**2/(2*(Ad_L_prop*d_obsdCLe(i))**2) ! prendre en compte erreurs mesurées
!                 end do
                
                
                
                
                
                logalpha(idis)=ndatad_R*(1/widening-1)*log(Ad_R)&
                +ndatad_L*(1/widening-1)*log(Ad_L)&
                +like/widening-like_alt& 
                -alphamax(idis)
                logalpha(idis)=min(logalpha(idis),0.)

            enddo
            
    
        endif

        npt_ani=0
        npt_iso=0
    
        do i=1,npt
            if (isoflag(i).eqv..true.) npt_iso=npt_iso+1
            if (isoflag(i).eqv..false.) npt_ani=npt_ani+1
        enddo
        if (npt_iso+npt_ani.ne.npt) stop "Error here"
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

        IF (ount.GT.burn_in) THEN
            sample=sample+1
            
            IF (mod(ount,thin)==0) THEN
                
                alpharef=exp(logalpharef)
                th = th + alpharef
                th_wide = th_wide + 1 
                histo(npt)=histo(npt)+alpharef !histogram of number of layers
                histo_wide(npt)=histo_wide(npt)+1
                
                alpha=exp(logalpha)
                alphasum=alphasum+alpha
                do idis=1,numdis
                    histo_alt(npt,idis)=histo_alt(npt,idis)+alpha(idis)
                enddo
                ! average dispersion curve and variance
                d_cRmoy=d_cRmoy+d_cR*alpharef
                d_cLmoy=d_cLmoy+d_cL*alpharef
                d_cRdelta=d_cRdelta+(d_cR**2)*alpharef
                d_cLdelta=d_cLdelta+(d_cL**2)*alpharef
                
                d_cRmoy_wide=d_cRmoy_wide+d_cR
                d_cLmoy_wide=d_cLmoy_wide+d_cL
                d_cRdelta_wide=d_cRdelta_wide+(d_cR**2)
                d_cLdelta_wide=d_cLdelta_wide+(d_cL**2)
                
                do idis=1,numdis
                    d_cRmoy_alt(:,idis)=d_cRmoy_alt(:,idis)+d_cR*alpha(idis)
                    d_cLmoy_alt(:,idis)=d_cLmoy_alt(:,idis)+d_cL*alpha(idis)
                    d_cRdelta_alt(:,idis)=d_cRdelta_alt(:,idis)+(d_cR**2)*alpha(idis)
                    d_cLdelta_alt(:,idis)=d_cLdelta_alt(:,idis)+(d_cL**2)*alpha(idis)
                enddo
                
!                 write(*,*)'logcratio',logcratio
!                 write(*,*)'like_alt',like_alt
!                 write(*,*)'logalpha',logalpha
!                 write(*,*)'alpha',alpha
!                 write(*,*)'logcratioref',logcratioref
!                 write(*,*)'like',like
!                 write(*,*)'logalpharef',logalpharef
!                 write(*,*)'alpharef',alpharef
                
                i_al=ceiling((logalpharef-logalpha_min)/(logalpha_max-logalpha_min)*num_logalpha)
                logalpharefhist(i_al)=logalpharefhist(i_al)+1
                
                do idis=1,numdis
                    i_al=ceiling((logalpha(idis)-logalpha_min)/(logalpha_max-logalpha_min)*num_logalpha)
                    if (logalpha(idis)<logalpha_min) i_al=1
                    if (logalpha(idis)>logalpha_max) i_al=num_logalpha
                    logalhist(i_al,idis)=logalhist(i_al,idis)+1
                enddo
                
                call combine(model_ref,nptref,nic_ref,noc_ref,voro,npt,d_max,&
                    r,rho,vpv,vph,vsv,vsh,qkappa,qshear,eta,nptfinal,nic,noc,xi,vpvsv_data)   
                j=1
                do i=disd,1,-1 ! average model
                    d=rearth-((i-1)*prof/real(disd-1))*1000.
                    
                    do while (r(j)<d)
                        j=j+1
                    end do
                    avvs(i)=avvs(i)+interp(r(j-1),r(j),vsv(j-1),vsv(j),d)*alpharef
                    avvs_wide(i)=avvs_wide(i)+interp(r(j-1),r(j),vsv(j-1),vsv(j),d)
                    avxi(i)=avxi(i)+interp(r(j-1),r(j),xi(j-1),xi(j),d)*alpharef
                    avxi_wide(i)=avxi_wide(i)+interp(r(j-1),r(j),xi(j-1),xi(j),d)
                    avvp(i)=avvp(i)+interp(r(j-1),r(j),vpvsv_data(j-1),vpvsv_data(j),d)*alpharef
                    avvp_wide(i)=avvp_wide(i)+interp(r(j-1),r(j),vpvsv_data(j-1),vpvsv_data(j),d)
                    
                    if (abs(interp(r(j-1),r(j),xi(j-1),xi(j),d)-1.)>0.01) then
                        probani(i)=probani(i)+alpharef
                        probani_wide(i)=probani_wide(i)+1
                    endif
                    
                    avvs_alt(i,:)=avvs_alt(i,:)+interp(r(j-1),r(j),vsv(j-1),vsv(j),d)*alpha
                    avxi_alt(i,:)=avxi_alt(i,:)+interp(r(j-1),r(j),xi(j-1),xi(j),d)*alpha
                    avvp_alt(i,:)=avvp_alt(i,:)+interp(r(j-1),r(j),vpvsv_data(j-1),vpvsv_data(j),d)*alpha
                        
                    if (abs(interp(r(j-1),r(j),xi(j-1),xi(j),d)-1.)>0.01) then
                        probani_alt(i,:)=probani_alt(i,:)+alpha
                    endif
                    
                end do
                
                
                ! Write Posterior matrix!-----------------------------------------------------------------
                j=1
                do i=disd,1,-1
                    d=rearth-((i-1)*prof/real(disd-1))*1000.
                    
                    do while (r(j)<d)
                        j=j+1
                    end do
                    ind_vsv=ceiling((interp(r(j-1),r(j),vsv(j-1),vsv(j),d)-vsref_min)*disv/(vsref_max-vsref_min))
                    if (ind_vsv==0) ind_vsv=1
                    ind_xi=ceiling(((interp(r(j-1),r(j),xi(j-1),&
                        xi(j),d))-xi_min)*disv/(xi_max-xi_min))
                    if (ind_xi==0) ind_xi=1
                    ind_vp=ceiling((interp(r(j-1),r(j),vpvsv_data(j-1),&
                        vpvsv_data(j),d)-vpvsv_min)*disv/(vpvsv_max-vpvsv_min))
                    if (ind_vp==0) ind_vp=1
                    if (ind_vsv<1)then
                        write(*,*)'ICI ICI ICI'
                        write(*,*)ind_vsv,i,interp(r(j-1),r(j),vsv(j-1),vsv(j),d),vsref_min,vsref_max,d
                        stop "vsv<disv"
                    endif
                    if (ind_vsv>disv) then
                        write(*,*)"vsv>disv"
                        write(*,*)interp(r(j-1),r(j),vsv(j-1),vsv(j),d),vsref_min,vsref_max
                        write(*,*)r(j-1),r(j),vsv(j-1),vsv(j)
                        write(*,*)npt
                        write(*,*)voro(:,2)
                        
                        stop "vsv>disv"
                    endif
                    if (ind_vp>disv) then
                        write(*,*)"vp>disv"
                        write(*,*)interp(r(j-1),r(j),vpvsv_data(j-1),vpvsv_data(j),d),vpvsv_min,vpvsv_max
                        write(*,*)voro(:,4)
                        write(*,*)r(j-1),r(j),vpv(j-1),vpv(j)
                        write(*,*)vpvsv_data
                        write(*,*)j,i
                        
                        stop "vp>disv"
                    endif
                    
                    if (ind_xi>disv) then
                        write(*,*)"xi>disv"
                        write(*,*)(interp(r(j-1),r(j),xi(j-1),xi(j),d)),xi_min,xi_max
                        write(*,*)voro(:,3)
                        write(*,*)xi
                        write(*,*)j,i
                        
                        stop "xi>disv"
                    endif
                    
                    postvs(i,ind_vsv)=postvs(i,ind_vsv)+alpharef
                    postvs_wide(i,ind_vsv)=postvs_wide(i,ind_vsv)+1
                    postvs_alt(i,ind_vsv,:)=postvs_alt(i,ind_vsv,:)+alpha
                    

                    if (ind_xi>0) then 
                        postxi(i,ind_xi)=postxi(i,ind_xi)+alpharef  
                        postxi_wide(i,ind_xi)=postxi_wide(i,ind_xi)+1  
                        postxi_alt(i,ind_xi,:)=postxi_alt(i,ind_xi,:)+alpha
                    endif  
                    if (ind_vp>0) then 
                        postvp(i,ind_vp)=postvp(i,ind_vp)+alpharef
                        postvp_wide(i,ind_vp)=postvp_wide(i,ind_vp)+1
                        postvp_alt(i,ind_vp,:)=postvp_alt(i,ind_vp,:)+alpha
                    endif    
                enddo

                i=ceiling((Ad_R-Ad_R_min)*disA/(Ad_R_max-Ad_R_min)) !histogram of Ad_R
                ML_Ad_R(i) = ML_Ad_R(i)+alpharef
                ML_Ad_R_wide(i) = ML_Ad_R_wide(i)+1
                ML_Ad_R_alt(i,:) = ML_Ad_R_alt(i,:)+alpha
        
                i=ceiling((Ad_L-Ad_L_min)*disA/(Ad_L_max-Ad_L_min)) !histogram of Ad_L
                ML_Ad_L(i) = ML_Ad_L(i)+alpharef
                ML_Ad_L_wide(i) = ML_Ad_L_wide(i)+1
                ML_Ad_L_alt(i,:) = ML_Ad_L_alt(i,:)+alpha
    
                !Get distribution on changepoint locations.
                
                do i=2,npt
                    j=ceiling((voro(i,1))*disd/(prof))
                    if (j.le.disd) then
                        histoch(j)=histoch(j)+alpharef
                        histoch_wide(j)=histoch_wide(j)+1
                        histoch_alt(j,:)=histoch_alt(j,:)+alpha
                    endif
                enddo    
                
                
                
                TRA(npt,npt_ani+1)=TRA(npt,npt_ani+1)+alpharef
                TRA_wide(npt,npt_ani+1)=TRA_wide(npt,npt_ani+1)+1
                TRA_alt(npt,npt_ani+1,:)=TRA_alt(npt,npt_ani+1,:)+alpha
                if(npt_ani>npt) stop "npt_ani>npt"
                if(npt>malayv) stop "npt>malay" ! MODIF 
                if(npt_ani>malayv) stop "npt_ani>malay" ! MODIF
                
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
 
 
            do i=1,disd
                if (histoch(i)<0) write(*,*)'hahahahahahahahahh' !lol
            enddo

 
            !**********************************************************************
    
            !       Display what is going on every "Display" samples

            !**********************************************************************

        endif
        IF ((mod(ount,display).EQ.0).and.(mod(ran,10).EQ.0)) THEN

            write(*,*)'processor number',ran+1,'/',nbproc
            write(*,*)'sample:',ount,'/',burn_in+nsample
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
            write(*,*)'recalculated',recalculated
            !write(*,*)Ar_birth_old,sigmav_old,sigmav
            write(*,*)'-----------------------------------------'
            write(*,*)
            write(*,*)
            
            
            recalculated=0
            
        END IF

        !write(*,*)ount

        !endif
    end do !End Markov chain
    
!     ! best fitting model
!     write(filenamemax,"('/bestmodel_',I3.3,'_.out')") rank    
!     open(56,file=dirname//filenamemax,status='replace')
!     write(56,*) nptmax,likemax
!     do i=1,npt
!         write(56,*)voromax(i,1),voromax(i,2),voromax(i,3),voromax(i,4)
!     enddo
!     close(56)
    
    k=0
    if (th.ne.0) then !normalize averages
        avvs=avvs/th
        avvp=avvp/th
        avxi=avxi/th
        probani=100*probani/th
        
        d_cRmoy=d_cRmoy/th
        d_cLmoy=d_cLmoy/th
        d_cRdelta=d_cRdelta/th
        d_cLdelta=d_cLdelta/th
        
        avvs_wide=avvs_wide/th_wide
        avvp_wide=avvp_wide/th_wide
        avxi_wide=avxi_wide/th_wide
        probani_wide=100*probani_wide/th_wide
        
        d_cRmoy_wide=d_cRmoy_wide/th_wide
        d_cLmoy_wide=d_cLmoy_wide/th_wide
        d_cRdelta_wide=d_cRdelta_wide/th_wide
        d_cLdelta_wide=d_cLdelta_wide/th_wide
        
        inorout(ran+1)=1
        k=k+1
        do idis=1,numdis
            avvs_alt(:,idis)=avvs_alt(:,idis)/alphasum(idis)
            avvp_alt(:,idis)=avvp_alt(:,idis)/alphasum(idis)
            avxi_alt(:,idis)=avxi_alt(:,idis)/alphasum(idis)
            probani_alt(:,idis)=100*probani_alt(:,idis)/alphasum(idis)
            
            d_cRmoy_alt(:,idis)=d_cRmoy_alt(:,idis)/alphasum(idis)
            d_cLmoy_alt(:,idis)=d_cLmoy_alt(:,idis)/alphasum(idis)
            d_cRdelta_alt(:,idis)=d_cRdelta_alt(:,idis)/alphasum(idis)
            d_cLdelta_alt(:,idis)=d_cLdelta_alt(:,idis)/alphasum(idis)
        enddo
    else
        write(*,*)
        write(*,*)'rank=',ran,'th=',th
        do k=325,350
            write(*,*)ran,probani(k),'%'
        enddo

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
    
    

    do i=1,disd
        if (histochs(i)<0) write(*,*)'haaaaaaaaaaaaaaaaaaaaaaaa'
    enddo
    !-!
    if (flag==1) then
        call MPI_REDUCE(convBa,convBas,nsample+burn_in,MPI_Real,MPI_Sum,0,MPI_COMM_small,ierror)
        call MPI_REDUCE(convB,convBs,nsample+burn_in,MPI_Real,MPI_Sum,0,MPI_COMM_small,ierror)
        call MPI_REDUCE(convDa,convDas,nsample+burn_in,MPI_Real,MPI_Sum,0,MPI_COMM_small,ierror)
        call MPI_REDUCE(convD,convDs,nsample+burn_in,MPI_Real,MPI_Sum,0,MPI_COMM_small,ierror)
        call MPI_REDUCE(convvp,convvps,nsample+burn_in,MPI_Real,MPI_Sum,0,MPI_COMM_small,ierror)
        call MPI_REDUCE(convvs,convvss,nsample+burn_in,MPI_Real,MPI_Sum,0,MPI_COMM_small,ierror)
        call MPI_REDUCE(convdp,convdps,nsample+burn_in,MPI_Real,MPI_Sum,0,MPI_COMM_small,ierror)
        call MPI_REDUCE(convxi,convxis,nsample+burn_in,MPI_Real,MPI_Sum,0,MPI_COMM_small,ierror)
        call MPI_REDUCE(histoch,histochs,disd,MPI_Real,MPI_Sum,0,MPI_COMM_small,ierror)
        call MPI_REDUCE(histoch_alt,histoch_alts,disd*numdis_max,MPI_Real,MPI_Sum,0,MPI_COMM_small,ierror)
        call MPI_REDUCE(histoch_wide,histoch_wides,disd,MPI_Real,MPI_Sum,0,MPI_COMM_small,ierror)
        call MPI_REDUCE(ML_Ad_R,ML_Ad_Rs,disA,MPI_Real,MPI_Sum,0,MPI_COMM_small,ierror)
        call MPI_REDUCE(ML_Ad_L,ML_Ad_Ls,disA,MPI_Real,MPI_Sum,0,MPI_COMM_small,ierror)
        call MPI_REDUCE(ML_Ad_R_alt,ML_Ad_R_alts,disA*numdis_max,MPI_Real,MPI_Sum,0,MPI_COMM_small,ierror)
        call MPI_REDUCE(ML_Ad_L_alt,ML_Ad_L_alts,disA*numdis_max,MPI_Real,MPI_Sum,0,MPI_COMM_small,ierror)
        call MPI_REDUCE(ML_Ad_R_wide,ML_Ad_R_wides,disA,MPI_Real,MPI_Sum,0,MPI_COMM_small,ierror)
        call MPI_REDUCE(ML_Ad_L_wide,ML_Ad_L_wides,disA,MPI_Real,MPI_Sum,0,MPI_COMM_small,ierror)
        call MPI_REDUCE(avvs,avvss,disd,MPI_Real,MPI_Sum,0,MPI_COMM_small,ierror)
        call MPI_REDUCE(avvp,avvps,disd,MPI_Real,MPI_Sum,0,MPI_COMM_small,ierror)
        call MPI_REDUCE(avxi,avxis,disd,MPI_Real,MPI_Sum,0,MPI_COMM_small,ierror)
        call MPI_REDUCE(probani,probanis,disd,MPI_Real,MPI_Sum,0,MPI_COMM_small,ierror)
        call MPI_REDUCE(postvs,postvss,disd*disv,MPI_Real,MPI_Sum,0,MPI_COMM_small,ierror)
        call MPI_REDUCE(postxi,postxis,disd*disv,MPI_Real,MPI_Sum,0,MPI_COMM_small,ierror)
        call MPI_REDUCE(postvp,postvps,disd*disv,MPI_Real,MPI_Sum,0,MPI_COMM_small,ierror)
        call MPI_REDUCE(ncell,ncells,nsample+burn_in,MPI_Real,MPI_Sum,0,MPI_COMM_small,ierror)
        call MPI_REDUCE(convd_R,convd_Rs,nsample+burn_in,MPI_Real,MPI_Sum,0,MPI_COMM_small,ierror)
        call MPI_REDUCE(convd_L,convd_Ls,nsample+burn_in,MPI_Real,MPI_Sum,0,MPI_COMM_small,ierror)
        call MPI_REDUCE(histo,histos,malay,MPI_Real,MPI_Sum,0,MPI_COMM_small,ierror)
        call MPI_REDUCE(histo_alt,histo_alts,malay*numdis_max,MPI_Real,MPI_Sum,0,MPI_COMM_small,ierror)
        call MPI_REDUCE(histo_wide,histo_wides,malay,MPI_Real,MPI_Sum,0,MPI_COMM_small,ierror)
        call MPI_REDUCE(TRA,TRAs,malay*(malay+1),MPI_Real,MPI_Sum,0,MPI_COMM_small,ierror)
        call MPI_REDUCE(TRA_alt,TRA_alts,malay*(malay+1)*numdis_max,MPI_Real,MPI_Sum,0,MPI_COMM_small,ierror)
        call MPI_REDUCE(TRA_wide,TRA_wides,malay*(malay+1),MPI_Real,MPI_Sum,0,MPI_COMM_small,ierror)
        call MPI_REDUCE(convAd_R,convAd_Rs,nsample+burn_in,MPI_Real,MPI_Sum,0,MPI_COMM_small,ierror)
        call MPI_REDUCE(convAd_L,convAd_Ls,nsample+burn_in,MPI_Real,MPI_Sum,0,MPI_COMM_small,ierror)
        ! 
        call MPI_REDUCE(d_cRmoy,d_cRmoys,ndatadmax,MPI_Real,MPI_Sum,0,MPI_COMM_small,ierror)
        call MPI_REDUCE(d_cLmoy,d_cLmoys,ndatadmax,MPI_Real,MPI_Sum,0,MPI_COMM_small,ierror)
        call MPI_REDUCE(d_cRdelta,d_cRdeltas,ndatadmax,MPI_Real,MPI_Sum,0,MPI_COMM_small,ierror)
        call MPI_REDUCE(d_cLdelta,d_cLdeltas,ndatadmax,MPI_Real,MPI_Sum,0,MPI_COMM_small,ierror)
        
        call MPI_REDUCE(d_cRmoy_alt,d_cRmoy_alts,ndatadmax*numdis_max,MPI_Real,MPI_Sum,0,MPI_COMM_small,ierror)
        call MPI_REDUCE(d_cLmoy_alt,d_cLmoy_alts,ndatadmax*numdis_max,MPI_Real,MPI_Sum,0,MPI_COMM_small,ierror)
        call MPI_REDUCE(d_cRdelta_alt,d_cRdelta_alts,ndatadmax*numdis_max,MPI_Real,MPI_Sum,0,MPI_COMM_small,ierror)
        call MPI_REDUCE(d_cLdelta_alt,d_cLdelta_alts,ndatadmax*numdis_max,MPI_Real,MPI_Sum,0,MPI_COMM_small,ierror)
        
        call MPI_REDUCE(d_cRmoy_wide,d_cRmoy_wides,ndatadmax,MPI_Real,MPI_Sum,0,MPI_COMM_small,ierror)
        call MPI_REDUCE(d_cLmoy_wide,d_cLmoy_wides,ndatadmax,MPI_Real,MPI_Sum,0,MPI_COMM_small,ierror)
        call MPI_REDUCE(d_cRdelta_wide,d_cRdelta_wides,ndatadmax,MPI_Real,MPI_Sum,0,MPI_COMM_small,ierror)
        call MPI_REDUCE(d_cLdelta_wide,d_cLdelta_wides,ndatadmax,MPI_Real,MPI_Sum,0,MPI_COMM_small,ierror)
        
        call MPI_REDUCE(avvs_alt,avvss_alt,disd*numdis_max,MPI_Real,MPI_Sum,0,MPI_COMM_small,ierror)
        call MPI_REDUCE(avvp_alt,avvps_alt,disd*numdis_max,MPI_Real,MPI_Sum,0,MPI_COMM_small,ierror)
        call MPI_REDUCE(avxi_alt,avxis_alt,disd*numdis_max,MPI_Real,MPI_Sum,0,MPI_COMM_small,ierror)
        call MPI_REDUCE(probani_alt,probanis_alt,disd*numdis_max,MPI_Real,MPI_Sum,0,MPI_COMM_small,ierror)
        call MPI_REDUCE(postvs_alt,postvss_alt,disd*disv*numdis_max,MPI_Real,MPI_Sum,0,MPI_COMM_small,ierror)
        call MPI_REDUCE(postxi_alt,postxis_alt,disd*disv*numdis_max,MPI_Real,MPI_Sum,0,MPI_COMM_small,ierror)
        call MPI_REDUCE(postvp_alt,postvps_alt,disd*disv*numdis_max,MPI_Real,MPI_Sum,0,MPI_COMM_small,ierror)
        call MPI_REDUCE(logalhist,logalhists,num_logalpha*numdis_max,MPI_Real,MPI_Sum,0,MPI_COMM_small,ierror)
        call MPI_REDUCE(logalpharefhist,logalpharefhists,num_logalpha,MPI_Real,MPI_Sum,0,MPI_COMM_small,ierror)
        
        call MPI_REDUCE(avvs_wide,avvss_wide,disd,MPI_Real,MPI_Sum,0,MPI_COMM_small,ierror)
        call MPI_REDUCE(avvp_wide,avvps_wide,disd,MPI_Real,MPI_Sum,0,MPI_COMM_small,ierror)
        call MPI_REDUCE(avxi_wide,avxis_wide,disd,MPI_Real,MPI_Sum,0,MPI_COMM_small,ierror)
        call MPI_REDUCE(probani_wide,probanis_wide,disd,MPI_Real,MPI_Sum,0,MPI_COMM_small,ierror)
        call MPI_REDUCE(postvs_wide,postvs_wides,disd*disv,MPI_Real,MPI_Sum,0,MPI_COMM_small,ierror)
        call MPI_REDUCE(postxi_wide,postxi_wides,disd*disv,MPI_Real,MPI_Sum,0,MPI_COMM_small,ierror)
        call MPI_REDUCE(postvp_wide,postvp_wides,disd*disv,MPI_Real,MPI_Sum,0,MPI_COMM_small,ierror)

        do i=1,disd
            if (histochs(i)<0) write(*,*)'hiiiiiiiiiiiiiiiii'
        enddo
        
        write(*,*)'here'

        d_cRmoys=d_cRmoys/j
        d_cLmoys=d_cLmoys/j
        d_cRdeltas=sqrt(d_cRdeltas/j-d_cRmoys**2)
        d_cLdeltas=sqrt(d_cLdeltas/j-d_cLmoys**2)
        
        d_cRmoy_alts=d_cRmoy_alts/j
        d_cLmoy_alts=d_cLmoy_alts/j
        d_cRdelta_alts=sqrt(d_cRdelta_alts/j-d_cRmoy_alts**2)
        d_cLdelta_alts=sqrt(d_cLdelta_alts/j-d_cLmoy_alts**2)
        
        d_cRmoy_wides=d_cRmoy_wides/j
        d_cLmoy_wides=d_cLmoy_wides/j
        d_cRdelta_wides=sqrt(d_cRdelta_wides/j-d_cRmoy_wides**2)
        d_cLdelta_wides=sqrt(d_cLdelta_wides/j-d_cLmoy_wides**2)
        
        convPs=convPs/j
        convBs=convBs/j
        convBas=convBas/j
        convDs=convDs/j
        convDas=convDas/j
        convvps=convvps/j
        convvss=convvss/j
        convdps=convdps/j
        convxis=convxis/j
 
        probanis=probanis/j
        probanis_alt=probanis_alt/j

        write(*,*)'rank=',ran,'th=',th

        avvss=avvss/j
        avvps=avvps/j
        avxis=avxis/j
        convd_Rs=convd_Rs/j
        convd_Ls=convd_Ls/j
        convAd_Rs=convAd_Rs/j
        convAd_Ls=convAd_Ls/j
        ncells=ncells/j
        
        avvss_wide=avvss_wide/j
        avvps_wide=avvps_wide/j
        avxis_wide=avxis_wide/j
        
        avvss_alt=avvss_alt/j
        avvps_alt=avvps_alt/j
        avxis_alt=avxis_alt/j

    endif
    
    !***********************************************************************

    !                      Write the results

    !***********************************************************************

    IF (ran==members(1)) THEN
        write(*,*)'writing'
        write(*,*)dirname//'/Change_points.out'
        open(65,file=dirname//'/Change_points.out',status='replace')
        do i=1,disd
            d=d_min+(i-0.5)*prof/real(disd)
            write(65,*)d,histochs(i)
        enddo
        close(65)
        
        open(65,file=dirname//'/Change_points_wide.out',status='replace')
        do i=1,disd
            d=d_min+(i-0.5)*prof/real(disd)
            write(65,*)d,histoch_wides(i)
        enddo
        close(65)
        
        open(65,file=dirname//'/Change_points_alt.out',status='replace')
        write(65,*)numdis
        do i=1,disd
            d=d_min+(i-0.5)*prof/real(disd)
            write(65,*)d,histoch_alts(i,:numdis)
        enddo
        close(65)

        open(95,file=dirname//'/Tradeoff.out',status='replace')
        write(95,*)malay
        do i=1,malay
            do j=1,malay+1
      
                write(95,*)TRAs(i,j)
            enddo
        enddo
        close(95)
        
        open(95,file=dirname//'/Tradeoff_wide.out',status='replace')
        write(95,*)malay
        do i=1,malay
            do j=1,malay+1
      
                write(95,*)TRA_wides(i,j)
            enddo
        enddo
        close(95)
        
        open(95,file=dirname//'/Tradeoff_alt.out',status='replace')
        write(95,*)malay,numdis
        do i=1,malay
            do j=1,malay+1
      
                write(95,*)TRA_alts(i,j,:numdis)
            enddo
        enddo
        close(95)

        open(56,file=dirname//'/Average.out',status='replace')
        do i=1,disd
            d=(i-1)*prof/real(disd-1)
            write(56,*)d,avvss(i),avxis(i),avvps(i),probanis(i)
        enddo
        close(56)
        
        open(56,file=dirname//'/Average_wide.out',status='replace')
        do i=1,disd
            d=(i-1)*prof/real(disd-1)
            write(56,*)d,avvss_wide(i),avxis_wide(i),avvps_wide(i),probanis_wide(i)
        enddo
        close(56)
        
        open(56,file=dirname//'/Average_alt.out',status='replace')
        write(56,*)numdis
        do i=1,disd
            d=(i-1)*prof/real(disd-1)
            do idis=1,numdis
                write(56,*)d,avvss_alt(i,idis),avxis_alt(i,idis),&
                avvps_alt(i,idis),probanis_alt(i,idis)
            enddo
        enddo
        close(56)

        open(46,file=dirname//'/Sigmad_R.out',status='replace')
        write(46,*)Ad_R_min,Ad_R_max,disa
        do i=1,disA
            d=Ad_R_min+(i-0.5)*(Ad_R_max-Ad_R_min)/real(disA)
            write(46,*)d,ML_Ad_Rs(i)
        enddo
        close(46)

        open(47,file=dirname//'/Sigmad_L.out',status='replace')
        write(47,*)Ad_L_min,Ad_L_max,disa
        do i=1,disA
            d=Ad_L_min+(i-0.5)*(Ad_L_max-Ad_L_min)/real(disA)
            write(47,*)d,ML_Ad_Ls(i)
        enddo
        close(47)
        
        open(46,file=dirname//'/Sigmad_R_alt.out',status='replace')
        write(46,*)numdis
        write(46,*)Ad_R_min,Ad_R_max,disa
        do i=1,disA
            d=Ad_R_min+(i-0.5)*(Ad_R_max-Ad_R_min)/real(disA)
            write(46,*)d,ML_Ad_R_alts(i,:numdis)
        enddo
        close(46)

        open(47,file=dirname//'/Sigmad_L_alt.out',status='replace')
        write(47,*)numdis
        write(47,*)Ad_L_min,Ad_L_max,disa
        do i=1,disA
            d=Ad_L_min+(i-0.5)*(Ad_L_max-Ad_L_min)/real(disA)
            write(47,*)d,ML_Ad_L_alts(i,:numdis)
        enddo
        close(47)
        
        open(46,file=dirname//'/Sigmad_R_wide.out',status='replace')
        write(46,*)Ad_R_min,Ad_R_max,disa
        do i=1,disA
            d=Ad_R_min+(i-0.5)*(Ad_R_max-Ad_R_min)/real(disA)
            write(46,*)d,ML_Ad_R_wides(i)
        enddo
        close(46)

        open(47,file=dirname//'/Sigmad_L_wide.out',status='replace')
        write(47,*)Ad_L_min,Ad_L_max,disa
        do i=1,disA
            d=Ad_L_min+(i-0.5)*(Ad_L_max-Ad_L_min)/real(disA)
            write(47,*)d,ML_Ad_L_wides(i)
        enddo
        close(47)

        open(71,file=dirname//'/Dispersion_mean.out',status='replace')
        write(71,*)ndatad_R,ndatad_L
        do i=1,ndatad_R
            write(71,*)peri_R(i),n_R(i),d_cRmoys(i),d_cRdeltas(i)
        enddo
        do i=1,ndatad_L
            write(71,*)peri_L(i),n_L(i),d_cLmoys(i),d_cLdeltas(i)
        enddo
        close(71)
        
        open(71,file=dirname//'/Dispersion_obs.out',status='replace')
        write(71,*)ndatad_R,ndatad_L
        do i=1,ndatad_R
            write(71,*)peri_R(i),n_R(i),d_obsdcR(i),d_obsdcRe(i)
        enddo
        do i=1,ndatad_L
            write(71,*)peri_L(i),n_L(i),d_obsdcL(i),d_obsdcLe(i)
        enddo
        close(71)
        
        
        open(71,file=dirname//'/Dispersion_mean_wide.out',status='replace')
        write(71,*)ndatad_R,ndatad_L
        do i=1,ndatad_R
            write(71,*)peri_R(i),n_R(i),d_cRmoy_wides(i),d_cRdelta_wides(i)
        enddo
        do i=1,ndatad_L
            write(71,*)peri_L(i),n_L(i),d_cLmoy_wides(i),d_cLdelta_wides(i)
        enddo
        close(71)
        
        open(71,file=dirname//'/Dispersion_mean_alt.out',status='replace')
        write(71,*)numdis
        write(71,*)ndatad_R,ndatad_L
        do j=1,numdis
            do i=1,ndatad_R
                write(71,*)peri_R(i),n_R(i),d_cRmoy_alts(i,j),d_cRdelta_alts(i,j)
            enddo
            do i=1,ndatad_L
                write(71,*)peri_L(i),n_L(i),d_cLmoy_alts(i,j),d_cLdelta_alts(i,j)
            enddo
        enddo
        close(71)
        
        open(71,file=dirname//'/Dispersion_obs_alt.out',status='replace')
        write(71,*)numdis
        write(71,*)ndatad_R,ndatad_L
        do j=1,numdis
            do i=1,ndatad_R
                write(71,*)peri_R(i),n_R(i),d_obsdcR_alt(i,j),d_obsdcRe_alt(i,j)
            enddo
            do i=1,ndatad_L
                write(71,*)peri_L(i),n_L(i),d_obsdcL_alt(i,j),d_obsdcLe_alt(i,j)
            enddo
        enddo
        close(71)

        open(71,file=dirname//'/Posterior.out',status='replace')
        write(71,*)prof,disd,d_max
        write(71,*)vsref_min,vsref_max,disv,width,xi_min,xi_max,vpvsv_min,vpvsv_max
        do i=1,disd
            do j=1,disv
                write(71,*)postvss(i,j),postxis(i,j),postvps(i,j)
            enddo
        enddo
        close(71)
        
        open(71,file=dirname//'/Posterior_wide.out',status='replace')
        write(71,*)prof,disd,d_max
        write(71,*)vsref_min,vsref_max,disv,width,xi_min,xi_max,vpvsv_min,vpvsv_max
        do i=1,disd
            do j=1,disv
                write(71,*)postvs_wides(i,j),postxi_wides(i,j),postvp_wides(i,j)
            enddo
        enddo
        close(71)
        
        open(71,file=dirname//'/Posterior_alt.out',status='replace')
        write(71,*)numdis
        write(71,*)prof,disd,d_max
        write(71,*)vsref_min,vsref_max,disv,width,xi_min,xi_max,vpvsv_min,&
        vpvsv_max
        do i=1,disd
            do j=1,disv
                do idis=1,numdis
                    write(71,*)postvss_alt(i,j,idis),postxis_alt(i,j,idis),&
                    postvps_alt(i,j,idis)
                enddo
            enddo
        enddo
        close(71)

        open(54,file=dirname//'/Convergence_misfit.out',status='replace')
        write(54,*)burn_in,nsample,burn_in,nsample
        do i=1,nsample+burn_in
            write(54,*)convd_R(i),convd_Rs(i),convd_L(i),convd_Ls(i)
        enddo
        close(54)

        open(53,file=dirname//'/Convergence_nb_layers.out',status='replace')
        write(53,*)burn_in,nsample
        do i=1,nsample+burn_in
            write(53,*)ncell(i),ncells(i)
        enddo
        close(53)

        open(52,file=dirname//'/Convergence_sigma_R.out',status='replace')
        write(52,*)burn_in,nsample
        do i=1,nsample+burn_in
            write(52,*)convAd_R(i),convAd_Rs(i)
        enddo
        close(52)

        open(53,file=dirname//'/Convergence_sigma_L.out',status='replace')
        write(53,*)burn_in,nsample
        do i=1,nsample+burn_in
            write(53,*)convAd_L(i),convAd_Ls(i)
        enddo
        close(53)


        open(45,file=dirname//'/NB_layers.out',status='replace')
        write(45,*)malay
        do i=1,malay
            write(45,*)histos(i)
        enddo
        close(45)
        
        open(45,file=dirname//'/NB_layers_wide.out',status='replace')
        write(45,*)malay
        do i=1,malay
            write(45,*)histo_wides(i)
        enddo
        close(45)
        
        open(45,file=dirname//'/NB_layers_alt.out',status='replace')
        write(45,*)malay,numdis
        do i=1,malay
            write(45,*)histo_alts(i,:numdis)
        enddo
        close(45)

        open(54,file=dirname//'/Convergence_Birth.out',status='replace')
        write(54,*)burn_in,nsample,burn_in,nsample
        do i=1,nsample+burn_in
            write(54,*)convB(i),convBs(i),convBa(i),convBas(i)
        enddo
        close(54)
        
        open(54,file=dirname//'/Convergence_Death.out',status='replace')
        write(54,*)burn_in,nsample,burn_in,nsample
        do i=1,nsample+burn_in
            write(54,*)convD(i),convDs(i),convDa(i),convDas(i)
        enddo
        close(54)
        
        open(54,file=dirname//'/Convergence_vp.out',status='replace')
        write(54,*)burn_in,nsample,burn_in,nsample
        do i=1,nsample+burn_in
            write(54,*)convvp(i),convvps(i)
        enddo
        close(54)
        
        open(54,file=dirname//'/Convergence_vs.out',status='replace')
        write(54,*)burn_in,nsample,burn_in,nsample
        do i=1,nsample+burn_in
            write(54,*)convvs(i),convvss(i)
        enddo
        close(54)
        
        open(54,file=dirname//'/Convergence_dp.out',status='replace')
        write(54,*)burn_in,nsample,burn_in,nsample
        do i=1,nsample+burn_in
            write(54,*)convdp(i),convdps(i)
        enddo
        close(54)
        
        open(54,file=dirname//'/Convergence_xi.out',status='replace')
        write(54,*)burn_in,nsample,burn_in,nsample
        do i=1,nsample+burn_in
            write(54,*)convxi(i),convxis(i)
        enddo
        close(54)
        
        open(56,file=dirname//'/alpha_hist.out',status='replace')
        write(56,*)numdis
        write(56,*)num_logalpha,logalpha_min,logalpha_max
        do i=1,num_logalpha
            do idis=1,numdis
                write(56,*)logalhists(i,idis)
            enddo
        enddo
        close(56)
        
        open(56,file=dirname//'/alpharef_hist.out',status='replace')
        write(56,*)num_logalpha,logalpha_min,logalpha_max
        do i=1,num_logalpha
            write(56,*)logalpharefhists(i)
        enddo
        close(56)
        
        open(56,file=dirname//'/parameters.out',status='replace')
        write(56,*)'burn-in',burn_in
        write(56,*)'nsample',nsample
        write(56,*)'nsample_widening',nsample_widening
        write(56,*)'burn-in_widening',burn_in_widening
        close(56)


    endif
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
