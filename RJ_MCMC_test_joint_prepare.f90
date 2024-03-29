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
    integer th
    integer nlims_cur,nlims_cur_diff
    logical ier, error_flag ! error flags for dispersion curves and minos_bran
    integer npt,npt_prop,npt_iso,npt_ani ! numbers of points, isotropic, anisotropic
    logical accept,tes,birth,birtha,death,deatha,move,value,noisd_R,noisd_L,ani,change_vp !which parameter did we change?
    logical testing
    integer ind2,ind_vp,ind_vsv,ind_xi,j !more indices, mostly for posterior. j sometimes is the number of processors, careful at the end!
    real d !for calculating posterior, mostly depth
    real peri_R(ndatadmax),peri_L(ndatadmax) !periods of data
    real peri_R_tmp(ndatadmax),peri_L_tmp(ndatadmax) !periods of data
    
    integer n_R(ndatadmax),n_L(ndatadmax),nmax_R,nmax_L,nmin_R,nmin_L !harmonic number of modes
    integer n_R_tmp(ndatadmax),n_L_tmp(ndatadmax)
    integer nlims_R(2,nharmo_max),nlims_L(2,nharmo_max),numharm_count,numharm_R(nharmo_max),numharm_L(nharmo_max)

    real PrB,AcB,PrD,AcD,Prnd_R,Acnd_R,Prnd_L,Acnd_L,&
        Acba,Prba,Acda,Prda,out,Prxi,Acxi,Pr_vp,&
        Ac_vp,PrP,PrV,AcP,AcV ! to adjust acceptance rates for all different parameters
    real lsd_L,lsd_L_prop,lsd_L_min,lsd_R,lsd_R_prop,lsd_R_min !logprob without uncertainties
    real logprob_vsv,logprob_vp !for reversible jump
    real voro(malay,4),vsref_min,vsref_max,vpref_min,vpref_max !the model with its bounds
    real voro_prop(malay,4)
    real t1,t2 !timers
    real like,like_prop,u,& !log-likelyhoods
        liked_R_prop,liked_R,liked_L_prop,liked_L
    double precision logrsig,Ad_R,Ad_R_prop,Ad_L,Ad_L_prop !uncertainty parameters
    real sigmav,sigmav_old,sigmav_new,AR_birth_old   !proposal on velocity when Birth move  
    real d_obsdcR(ndatadmax),d_obsdcL(ndatadmax),d_obsdcRe(ndatadmax),d_obsdcLe(ndatadmax) !observed data
    integer inorout(21000),inorouts(21000),members(21000),io ! mpi related variables
    integer ndatad_R,ndatad_L,ndatad_R_tmp,ndatad_L_tmp !number of observed data points
    real pxi,p_vp,pd,pv,pAd_R,pAd_L ! related to observed data for azimuthal anisotropy
    ! Geometry parameters
    ! Traces
    integer nptref,malayv ! number of points in reference model, number of layers in voro

    logical isoflag(malay),isoflag_prop(malay) ! is this layer isotropic?
    real d_cR(ndatadmax),d_cL(ndatadmax),d_cR_tmp(ndatadmax),d_cL_tmp(ndatadmax) ! phase velocity as simulated by forward modelling
    real rq_R(ndatadmax),rq_L(ndatadmax) ! errors from modelling
    real d_cR_prop(ndatadmax),d_cL_prop(ndatadmax),d_cR_prop_tmp(ndatadmax),d_cL_prop_tmp(ndatadmax) ! phase velocity as simulated by forward modelling
    real rq_R_prop(ndatadmax),rq_L_prop(ndatadmax)
    !for MPI
    integer ra,ran,rank, nbproc, ierror ,status(MPI_STATUS_SIZE),group_world,good_group,MPI_COMM_small,flag
    ! ierror: MPI-related error status, nbproc: MPI-related, number of processes, rank, group_world: MPI-related
    character filename*13, number*4
    real likemax ! log-likelihood of best model
    integer nptmax ! number of layers of the best model
    character filenamemax*30  !filename for writing the best model
    real model_ref(mk,9) ! reference model, all models are deviations from it
    real,dimension(mk) :: r,rho,vpv,vph,vsv,vsh,qkappa,qshear,eta,xi,vpvsv_data ! input for forward modelling
    integer nptfinal,nic,noc,nic_ref,noc_ref,jcom ! more inputs
    integer nmodes,n_mode(nmodes_max),l_mode(nmodes_max) ! outputs of forward modelling
    real c_ph(nmodes_max),period(nmodes_max),raylquo(nmodes_max),tref ! more outputs, rayquo: error measurement, should be of the order of eps

    real wmin_R(nharmo_max),wmax_R(nharmo_max),wmin_L(nharmo_max),wmax_L(nharmo_max) ! more inputs f
    logical stuck
    integer nharm_R,nharm_L,iharm
    
    integer idis, numdis,i_al ! additional stuff for joint inversion
    real like_alt,like_alt_L,like_alt_R
    real logcratio(numdis_max),logalpharef,logalpha(numdis_max),alphasum(numdis_max),alpharefsum,logcratioref
    real alpha(numdis_max),alpharef
    real postvp_alt(disd,disv,numdis_max),postvs_alt(disd,disv,numdis_max),postxi_alt(disd,disv,numdis_max)
    real postvps_alt(disd,disv,numdis_max),postvss_alt(disd,disv,numdis_max),postxis_alt(disd,disv,numdis_max) 
    
    real avvs_alt(disd,numdis_max),avvp_alt(disd,numdis_max),avxi_alt(disd,numdis_max)
    real avxis_alt(disd,numdis_max),avvss_alt(disd,numdis_max),avvps_alt(disd,numdis_max)
    real d_obsdcR_alt(ndatadmax,numdis_max),d_obsdcL_alt(ndatadmax,numdis_max)
    real d_obsdcRe_alt(ndatadmax,numdis_max),d_obsdcLe_alt(ndatadmax,numdis_max) !observed data
    real like_w,like_prop_w
    
    integer i_w
    real widening_prop,widening
    real mean_best(numdis_max),mean_prop(numdis_max),mean_props(numdis_max),meandiff_hist(n_w,numdis_max)
    real alphahist(num_logalpha,numdis_max,n_w),alphahists(num_logalpha,numdis_max,n_w)
    real alpharefhist(num_logalpha,n_w),alpharefhists(num_logalpha,n_w)
    real alphamax_hist(n_w,numdis_max),alphamaxref_hist(n_w)
    real alphamax(numdis_max),alpharefmax
    real alphamax_true(n_w,numdis_max),alpharefmax_true(n_w)
    real alphamax_prop(numdis_max),alpharefmax_prop
    real alphamax_props(numdis_max),alpharefmax_props
    real alphaall(nsample_widening/thin,numdis_max)
    real alpharefall(nsample_widening/thin)
    real, dimension(:), allocatable :: lat, lon
    !real lat(numdis_max),lon(numdis_max)
    integer isum
    real totsum,sum_tmp
    
    character filebycore*15
    integer outputfile
    
    ! to save the current state
    real best_current_Ad_R,best_current_Ad_L,best_current_pxi,best_current_p_vp
    real best_current_pAd_R,best_current_pAd_L,best_current_pv
    real best_current_pd,best_current_sigmav
    integer best_current_npt
    real best_current_voro(malay,4)
    logical burnin_in_progress
    
    real ncell(n_w*(burn_in_widening+nsample_widening)+burn_in)
    real ncells(n_w*(burn_in_widening+nsample_widening)+burn_in)
    real convBs(n_w*(burn_in_widening+nsample_widening)+burn_in)
    real convB(n_w*(burn_in_widening+nsample_widening)+burn_in)
    real convPs(n_w*(burn_in_widening+nsample_widening)+burn_in)
    real convP(n_w*(burn_in_widening+nsample_widening)+burn_in) ! birth rates, vp change rates
    real convvs(n_w*(burn_in_widening+nsample_widening)+burn_in)
    real convvss(n_w*(burn_in_widening+nsample_widening)+burn_in)
    real convdp(n_w*(burn_in_widening+nsample_widening)+burn_in)
    real convdps(n_w*(burn_in_widening+nsample_widening)+burn_in)
    real convvp(n_w*(burn_in_widening+nsample_widening)+burn_in)
    real convvps(n_w*(burn_in_widening+nsample_widening)+burn_in)
    real convxis(n_w*(burn_in_widening+nsample_widening)+burn_in)
    real convxi(n_w*(burn_in_widening+nsample_widening)+burn_in)
    real convBas(n_w*(burn_in_widening+nsample_widening)+burn_in)
    real convBa(n_w*(burn_in_widening+nsample_widening)+burn_in) ! anisotropic birth rates
    real convDs(n_w*(burn_in_widening+nsample_widening)+burn_in)
    real convD(n_w*(burn_in_widening+nsample_widening)+burn_in)
    real convDas(n_w*(burn_in_widening+nsample_widening)+burn_in)
    real convDa(n_w*(burn_in_widening+nsample_widening)+burn_in)
    real convd_R(n_w*(burn_in_widening+nsample_widening)+burn_in)
    real convd_L(n_w*(burn_in_widening+nsample_widening)+burn_in)
    real convd_Rs(n_w*(burn_in_widening+nsample_widening)+burn_in)
    real convd_Ls(n_w*(burn_in_widening+nsample_widening)+burn_in) ! the proposed model, history of lsd
    real convAd_R(n_w*(burn_in_widening+nsample_widening)+burn_in)
    real convAd_Rs(n_w*(burn_in_widening+nsample_widening)+burn_in)
    real convAd_L(n_w*(burn_in_widening+nsample_widening)+burn_in)
    real convAd_Ls(n_w*(burn_in_widening+nsample_widening)+burn_in) !variations of uncertainty parameters
    
    integer i_opt,num_prop(numdis_max),num_props(numdis_max)
    
    ! todo: implement a test with raylquo
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
    
    
    
    testing=.false.
    if (testing) write(*,*)'testing with synthetic model'
    
    ra=rank
    ran=rank
    
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
    
    
    if (testing) then !!!!!!!testing: create synthetic model
        ndatad_R=0
        j=1
        numharm_count=1
        ! careful: periods listed by decreasing period, increasing harmonic order
        nlims_R(1,numharm_count)=j
        do i=200,40,-5 !synthetic periods, harmonics, uncertainties
            peri_R(j)=real(i)
            d_obsdcRe(j)=0.01
            n_R(j)=0
            j=j+1
        end do
        numharm_R(numharm_count)=0
        nlims_R(2,numharm_count)=j-1
        wmin_R(numharm_count)=1000./(maxval(peri_R(nlims_R(1,numharm_count):nlims_R(2,numharm_count)))+20)
        wmax_R(numharm_count)=1000./(minval(peri_R(nlims_R(1,numharm_count):nlims_R(2,numharm_count)))-2)
        
!         numharm_count=numharm_count+1
!         nlims_R(1,numharm_count)=j
!         do i=240,40,-10
!             
!             peri_R(j)=real(i)
!             d_obsdcRe(j)=0.01
!             n_R(j)=1
!             j=j+1
!         end do
!         numharm_R(numharm_count)=1
!         nlims_R(2,numharm_count)=j-1
!         wmin_R(numharm_count)=1000./(maxval(peri_R(nlims_R(1,numharm_count):nlims_R(2,numharm_count)))+20)
!         wmax_R(numharm_count)=1000./(minval(peri_R(nlims_R(1,numharm_count):nlims_R(2,numharm_count)))-2)

!         numharm_count=numharm_count+1
!         nlims_R(1,numharm_count)=j
!         do i=160,40,-10
!             
!             peri_R(j)=real(i)
!             d_obsdcRe(j)=0.01
!             n_R(j)=2
!             j=j+1
!         end do
!         numharm_R(numharm_count)=2
!         nlims_R(2,numharm_count)=j-1
!         wmin_R(numharm_count)=1000./(maxval(peri_R(nlims_R(1,numharm_count):nlims_R(2,numharm_count)))+20)
!         wmax_R(numharm_count)=1000./(minval(peri_R(nlims_R(1,numharm_count):nlims_R(2,numharm_count)))-2)

!         numharm_count=numharm_count+1
!         nlims_R(1,numharm_count)=j
!         do i=90,40,-10
!             
!             peri_R(j)=real(i)
!             d_obsdcRe(j)=0.01
!             n_R(j)=3
!             j=j+1
!         end do
!         numharm_R(numharm_count)=3
!         nlims_R(2,numharm_count)=j-1
!         wmin_R(numharm_count)=1000./(maxval(peri_R(nlims_R(1,numharm_count):nlims_R(2,numharm_count)))+20)
!         wmax_R(numharm_count)=1000./(minval(peri_R(nlims_R(1,numharm_count):nlims_R(2,numharm_count)))-2)

!         numharm_count=numharm_count+1
!         nlims_R(1,numharm_count)=j
!         do i=50,40,-5
!             
!             peri_R(j)=real(i)
!             d_obsdcRe(j)=0.01
!             n_R(j)=4
!             j=j+1
!         end do
!         numharm_R(numharm_count)=4
!         nlims_R(2,numharm_count)=j-1
!         wmin_R(numharm_count)=1000./(maxval(peri_R(nlims_R(1,numharm_count):nlims_R(2,numharm_count)))+20)
!         wmax_R(numharm_count)=1000./(minval(peri_R(nlims_R(1,numharm_count):nlims_R(2,numharm_count)))-2)

!         numharm_count=numharm_count+1
!         nlims_R(1,numharm_count)=j
!         do i=50,40,-5
!             
!             peri_R(j)=real(i)
!             d_obsdcRe(j)=0.01
!             n_R(j)=5
!             j=j+1
!         end do
!         numharm_R(numharm_count)=5
!         nlims_R(2,numharm_count)=j-1
!         wmin_R(numharm_count)=1000./(maxval(peri_R(nlims_R(1,numharm_count):nlims_R(2,numharm_count)))+20)
!         wmax_R(numharm_count)=1000./(minval(peri_R(nlims_R(1,numharm_count):nlims_R(2,numharm_count)))-2)

        ndatad_R=j-1
        nharm_R=numharm_count
!         wmin_R=1000./(maxval(peri_R(:ndatad_R))+20)
!         wmax_R=1000./(minval(peri_R(:ndatad_R))-2)
        nmin_R=minval(n_R(:ndatad_R))
        nmax_R=maxval(n_R(:ndatad_R))
        
        ndatad_L=0
        j=1
        numharm_count=1
        nlims_L(1,numharm_count)=j
        do i=200,40,-5
            peri_L(j)=real(i)
            d_obsdcLe(j)=0.02
            n_L(j)=0
            j=j+1
        end do
        numharm_L(numharm_count)=0
        nlims_L(2,numharm_count)=j-1
        wmin_L(numharm_count)=1000./(maxval(peri_L(nlims_L(1,numharm_count):nlims_L(2,numharm_count)))+20)
        wmax_L(numharm_count)=1000./(minval(peri_L(nlims_L(1,numharm_count):nlims_L(2,numharm_count)))-2)
        
!         numharm_count=numharm_count+1
!         nlims_L(1,numharm_count)=j
!         do i=240,40,-10
!             peri_L(j)=real(i)
!             d_obsdcLe(j)=0.02
!             n_L(j)=1
!             j=j+1
!         end do
!         numharm_L(numharm_count)=1
!         nlims_L(2,numharm_count)=j-1
!         wmin_L(numharm_count)=1000./(maxval(peri_L(nlims_L(1,numharm_count):nlims_L(2,numharm_count)))+20)
!         wmax_L(numharm_count)=1000./(minval(peri_L(nlims_L(1,numharm_count):nlims_L(2,numharm_count)))-2)

!         numharm_count=numharm_count+1
!         nlims_L(1,numharm_count)=j
!         do i=160,40,-10
!             peri_L(j)=real(i)
!             d_obsdcLe(j)=0.02
!             n_L(j)=2
!             j=j+1
!         end do
!         numharm_L(numharm_count)=2
!         nlims_L(2,numharm_count)=j-1
!         wmin_L(numharm_count)=1000./(maxval(peri_L(nlims_L(1,numharm_count):nlims_L(2,numharm_count)))+20)
!         wmax_L(numharm_count)=1000./(minval(peri_L(nlims_L(1,numharm_count):nlims_L(2,numharm_count)))-2)

!         numharm_count=numharm_count+1
!         nlims_L(1,numharm_count)=j
!         do i=90,40,-10
!             peri_L(j)=real(i)
!             d_obsdcLe(j)=0.02
!             n_L(j)=3
!             j=j+1
!         end do
!         numharm_L(numharm_count)=3
!         nlims_L(2,numharm_count)=j-1
!         wmin_L(numharm_count)=1000./(maxval(peri_L(nlims_L(1,numharm_count):nlims_L(2,numharm_count)))+20)
!         wmax_L(numharm_count)=1000./(minval(peri_L(nlims_L(1,numharm_count):nlims_L(2,numharm_count)))-2)

!         numharm_count=numharm_count+1
!         nlims_L(1,numharm_count)=j
!         do i=50,40,-5
!             peri_L(j)=real(i)
!             d_obsdcLe(j)=0.02
!             n_L(j)=4
!             j=j+1
!         end do
!         numharm_L(numharm_count)=4
!         nlims_L(2,numharm_count)=j-1
!         wmin_L(numharm_count)=1000./(maxval(peri_L(nlims_L(1,numharm_count):nlims_L(2,numharm_count)))+20)
!         wmax_L(numharm_count)=1000./(minval(peri_L(nlims_L(1,numharm_count):nlims_L(2,numharm_count)))-2)

!         numharm_count=numharm_count+1
!         nlims_L(1,numharm_count)=j
!         do i=50,40,-5
!             peri_L(j)=real(i)
!             d_obsdcLe(j)=0.02
!             n_L(j)=5
!             j=j+1
!         end do
!         numharm_L(numharm_count)=5
!         nlims_L(2,numharm_count)=j-1
!         wmin_L(numharm_count)=1000./(maxval(peri_L(nlims_L(1,numharm_count):nlims_L(2,numharm_count)))+20)
!         wmax_L(numharm_count)=1000./(minval(peri_L(nlims_L(1,numharm_count):nlims_L(2,numharm_count)))-2)
        ndatad_L=j-1
        nharm_L=numharm_count
        wmin_L=1000./(maxval(peri_L(:ndatad_L))+20)
        wmax_L=1000./(minval(peri_L(:ndatad_L))-2)
        nmin_L=minval(n_L(:ndatad_L))
        nmax_L=maxval(n_L(:ndatad_L))
        
        tref=sum(peri_R(:ndatad_R))/ndatad_R ! average period for minos
        !tref=sum(peri_L(:ndatad_L))/ndatad_L
        
        ! create synthetic model
        npt=0
        
        voro(1,1)=1 !depth of interface
        voro(1,2)=0.0  !vsv=vsv_prem*(1+voro(i,2))
        voro(1,3)=-1!0.7+0.5/33.*i !xi=voro(i,3), -1 for isotropic layer
        voro(1,4)=0!0.3-0.5/33.*i !vpv=vsv*vpvs*(1+voro(i,4)), vpvs set to 1.73 in
        npt=1
        do i=2,11
            voro(i,1)=30*(i-1) !depth of interface
            voro(i,2)=0.0  !vsv=vsv_prem*(1+voro(i,2))
            voro(i,3)=-1!0.7+0.5/33.*i !xi=voro(i,3), -1 for isotropic layer
            voro(i,4)=0!0.3-0.5/33.*i !vpv=vsv*vpvs*(1+voro(i,4)), vpvs set to 1.73 in params.h
            npt=npt+1
        end do
        
!         voro(10,4)=0.3
!         voro(11,4)=0.3
!         voro(12,4)=0.3
!         voro(13,4)=0.3
!         voro(14,4)=0.3
!         voro(15,4)=-0.3
!         voro(16,4)=-0.3
!         voro(17,4)=-0.3
!         voro(18,4)=-0.3
!         voro(19,4)=-0.3
        
        ! take voro, combine it with prem into a format suitable for minos
        call combine_linear(model_ref,nptref,nic_ref,noc_ref,voro,npt,d_max,&
            r,rho,vpv,vph,vsv,vsh,qkappa,qshear,eta,nptfinal,nic,noc,xi,vpvsv_data)
        
        !calculate synthetic dispersion curves
        if (ndatad_R>0) then
            jcom=3 !rayleigh waves
            
            do iharm=1,nharm_R
                nmodes=0
                call minos_bran(1,tref,nptfinal,nic,noc,r,rho,vpv,vph,vsv,vsh,&
                    qkappa,qshear,eta,eps,wgrav,jcom,lmin,lmax,&
                    wmin_R(iharm),wmax_R(iharm),numharm_R(iharm),numharm_R(iharm),&
                    nmodes_max,nmodes,n_mode,l_mode,c_ph,period,raylquo,error_flag) ! calculate normal modes
                if (error_flag) stop "INVALID INITIAL MODEL - RAYLEIGH - minos_bran.f FAIL 001"
                peri_R_tmp=peri_R(nlims_R(1,iharm):nlims_R(2,iharm))
                ndatad_R_tmp=nlims_R(2,iharm)-nlims_R(1,iharm)+1 ! fortran slices take the first and the last element
                n_R_tmp(:ndatad_R_tmp)=n_R(nlims_R(1,iharm):nlims_R(2,iharm))
                
                call dispersion_minos(nmodes_max,nmodes,n_mode,c_ph,period,raylquo,&
                    peri_R_tmp,n_R_tmp,d_cR_tmp,rq_R,ndatad_R_tmp,ier) ! extract phase velocities from minos output (pretty ugly)
                d_cR(nlims_R(1,iharm):nlims_R(2,iharm))=d_cR_tmp
                if (ier) stop "INVALID INITIAL MODEL - RAYLEIGH - CHANGE PERIODS or MODEL"
                if (maxval(abs(rq_R(:ndatad_R_tmp)))>maxrq*eps) stop "INVALID INITIAL MODEL - RAYLEIGH - CHANGE PERIODS or MODEL"
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
                if (error_flag) stop "INVALID INITIAL MODEL - LOVE - minos_bran.f FAIL 002"
                peri_L_tmp=peri_L(nlims_L(1,iharm):nlims_L(2,iharm))
                ndatad_L_tmp=nlims_L(2,iharm)-nlims_L(1,iharm)+1 ! fortran slices take the first and the last element
                n_L_tmp(:ndatad_L_tmp)=n_L(nlims_L(1,iharm):nlims_L(2,iharm))
                
                call dispersion_minos(nmodes_max,nmodes,n_mode,c_ph,period,raylquo,&
                    peri_L_tmp,n_L_tmp,d_cL_tmp,rq_L,ndatad_L_tmp,ier) ! extract phase velocities from minos output (pretty ugly)
                d_cL(nlims_L(1,iharm):nlims_L(2,iharm))=d_cL_tmp
                if (ier) stop "ier INVALID INITIAL MODEL - LOVE - CHANGE PERIODS or MODEL"
                if (maxval(abs(rq_L(:ndatad_L_tmp)))>maxrq*eps) stop "rq INVALID INITIAL MODEL - LOVE - CHANGE PERIODS or MODEL"
            enddo
        endif

        d_obsdcR(:ndatad_R)=d_cR(:ndatad_R)        
        d_obsdcL(:ndatad_L)=d_cL(:ndatad_L)
        
        !add errors
        
        IF (nbproc.gt.1) THEN
            IF (ran==0) THEN
                do i=1,ndatad_R
                    d_obsdcR(i)=d_obsdcR(i)+gasdev(ra)*0.1
                end do
                do i=1,ndatad_L
                    d_obsdcL(i)=d_obsdcL(i)+gasdev(ra)*0.1
                end do
                do i=2,nbproc
                    call MPI_SEND(d_obsdcR, ndatadmax, MPI_Real, i-1, 1, MPI_COMM_WORLD, ierror)
                    call MPI_SEND(d_obsdcL, ndatadmax, MPI_Real, i-1, 1, MPI_COMM_WORLD, ierror)
                enddo
            else
                call MPI_RECV(d_obsdcR, ndatadmax, MPI_Real, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierror)
                call MPI_RECV(d_obsdcL, ndatadmax, MPI_Real, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierror)
            endif
        else
            do i=1,ndatad_R
                d_obsdcR(i)=d_obsdcR(i)+gasdev(ra)*0.1
            end do
            do i=1,ndatad_L
                d_obsdcL(i)=d_obsdcL(i)+gasdev(ra)*0.1
            end do
        endif
        
        ! write synthetic model into a file
        open(65,file=dirname//'/true_model.out',status='replace')
        do i=1,nptfinal
            write(65,*)(rearth-r(i))/1000,vsv(i),xi(i),vpvsv_data(i)
        enddo
        close(65)
        
        allocate(lat(numdis_max))
        allocate(lon(numdis_max))
        ! make small deviations for joint inversion
        d_obsdcRe_alt(:,1)=d_obsdcRe
        d_obsdcLe_alt(:,1)=d_obsdcLe
        
        ! create synthetic model
        lat(1)=45
        lon(1)=45
        npt=0
        
        voro(1,1)=1 !depth of interface
        voro(1,2)=0.0  !vsv=vsv_prem*(1+voro(i,2))
        voro(1,3)=-1!0.7+0.5/33.*i !xi=voro(i,3), -1 for isotropic layer
        voro(1,4)=0!0.3-0.5/33.*i !vpv=vsv*vpvs*(1+voro(i,4)), vpvs set to 1.73 in
        do i=2,11
            voro(i,1)=30*(i-1) !depth of interface
            voro(i,2)=0.0  !vsv=vsv_prem*(1+voro(i,2))
            voro(i,3)=-1!0.7+0.5/33.*i !xi=voro(i,3), -1 for isotropic layer
            voro(i,4)=0!0.3-0.5/33.*i !vpv=vsv*vpvs*(1+voro(i,4)), vpvs set to 1.73 in params.h
            npt=npt+1
        end do
        
        voro(5,2)=voro(5,2)+0.05
        voro(6,2)=voro(6,2)-0.05
!         do i=5,10
!             voro(i,2)=voro(i,2)+0.05
!         enddo
!         voro(10,4)=0.3
!         voro(11,4)=0.3
!         voro(12,4)=0.3
!         voro(13,4)=0.3
!         voro(14,4)=0.3
!         voro(15,4)=-0.3
!         voro(16,4)=-0.3
!         voro(17,4)=-0.3
!         voro(18,4)=-0.3
!         voro(19,4)=-0.3
        
        
        ! take voro, combine it with prem into a format suitable for minos
        call combine_linear(model_ref,nptref,nic_ref,noc_ref,voro,npt,d_max,&
            r,rho,vpv,vph,vsv,vsh,qkappa,qshear,eta,nptfinal,nic,noc,xi,vpvsv_data)
        
        !calculate synthetic dispersion curves
        if (ndatad_R>0) then
            jcom=3 !rayleigh waves
            
            do iharm=1,nharm_R
                nmodes=0
                call minos_bran(1,tref,nptfinal,nic,noc,r,rho,vpv,vph,vsv,vsh,&
                    qkappa,qshear,eta,eps,wgrav,jcom,lmin,lmax,&
                    wmin_R(iharm),wmax_R(iharm),numharm_R(iharm),numharm_R(iharm),&
                    nmodes_max,nmodes,n_mode,l_mode,c_ph,period,raylquo,error_flag) ! calculate normal modes
                if (error_flag) stop "INVALID INITIAL MODEL - RAYLEIGH - minos_bran.f FAIL 001"
                peri_R_tmp=peri_R(nlims_R(1,iharm):nlims_R(2,iharm))
                ndatad_R_tmp=nlims_R(2,iharm)-nlims_R(1,iharm)+1 ! fortran slices take the first and the last element
                n_R_tmp(:ndatad_R_tmp)=n_R(nlims_R(1,iharm):nlims_R(2,iharm))
                
                call dispersion_minos(nmodes_max,nmodes,n_mode,c_ph,period,raylquo,&
                    peri_R_tmp,n_R_tmp,d_cR_tmp,rq_R,ndatad_R_tmp,ier) ! extract phase velocities from minos output (pretty ugly)
                d_cR(nlims_R(1,iharm):nlims_R(2,iharm))=d_cR_tmp
                if (ier) stop "INVALID INITIAL MODEL - RAYLEIGH - CHANGE PERIODS or MODEL"
                if (maxval(abs(rq_R(:ndatad_R_tmp)))>maxrq*eps) stop "INVALID INITIAL MODEL - RAYLEIGH - CHANGE PERIODS or MODEL"
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
                if (error_flag) stop "INVALID INITIAL MODEL - LOVE - minos_bran.f FAIL 002"
                peri_L_tmp=peri_L(nlims_L(1,iharm):nlims_L(2,iharm))
                ndatad_L_tmp=nlims_L(2,iharm)-nlims_L(1,iharm)+1 ! fortran slices take the first and the last elementwrite(*,*)'now here',rank
                n_L_tmp(:ndatad_L_tmp)=n_L(nlims_L(1,iharm):nlims_L(2,iharm))
                
                call dispersion_minos(nmodes_max,nmodes,n_mode,c_ph,period,raylquo,&
                    peri_L_tmp,n_L_tmp,d_cL_tmp,rq_L,ndatad_L_tmp,ier) ! extract phase velocities from minos output (pretty ugly)
                d_cL(nlims_L(1,iharm):nlims_L(2,iharm))=d_cL_tmp
                if (ier) stop "ier INVALID INITIAL MODEL - LOVE - CHANGE PERIODS or MODEL"
                if (maxval(abs(rq_L(:ndatad_L_tmp)))>maxrq*eps) stop "rq INVALID INITIAL MODEL - LOVE - CHANGE PERIODS or MODEL"
            enddo
        endif
        
        d_obsdcR_alt(:ndatad_R,1)=d_cR(:ndatad_R)
        d_obsdcL_alt(:ndatad_L,1)=d_cL(:ndatad_L)
        
        !add errors
        
        IF (nbproc.gt.1) THEN
            IF (ran==0) THEN
                do i=1,ndatad_R
                    d_obsdcR_alt(i,1)=d_obsdcR_alt(i,1)+gasdev(ra)*0.1
                end do
                do i=1,ndatad_L
                    d_obsdcL_alt(i,1)=d_obsdcL_alt(i,1)+gasdev(ra)*0.1
                end do
!                 do i=2,nbproc
!                     call MPI_SEND(d_obsdcR_alt, ndatadmax*numdis_max, MPI_Real, i-1, 2, MPI_COMM_WORLD, ierror)
!                     call MPI_SEND(d_obsdcL_alt, ndatadmax*numdis_max, MPI_Real, i-1, 2, MPI_COMM_WORLD, ierror)
!                 enddo
!             else
!                 call MPI_RECV(d_obsdcR_alt, ndatadmax*numdis_max, MPI_Real, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierror)
!                 call MPI_RECV(d_obsdcL_alt, ndatadmax*numdis_max, MPI_Real, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierror)
            endif
        else
            do i=1,ndatad_R
                d_obsdcR_alt(i,1)=d_obsdcR_alt(i,1)+gasdev(ra)*0.1
            end do
            do i=1,ndatad_L
                d_obsdcL_alt(i,1)=d_obsdcL_alt(i,1)+gasdev(ra)*0.1
            end do
        endif
        
        open(65,file=dirname//'/true_model_1.out',status='replace')
        do i=1,nptfinal
            write(65,*)(rearth-r(i))/1000,vsv(i),xi(i),vpvsv_data(i)
        enddo
        close(65)
        
        d_obsdcRe_alt(:,2)=d_obsdcRe
        d_obsdcLe_alt(:,2)=d_obsdcLe
        
        ! create synthetic model
        lat(2)=0
        lon(2)=0
        npt=0
        
        voro(1,1)=1 !depth of interface
        voro(1,2)=0.0  !vsv=vsv_prem*(1+voro(i,2))
        voro(1,3)=-1!0.7+0.5/33.*i !xi=voro(i,3), -1 for isotropic layer
        voro(1,4)=0!0.3-0.5/33.*i !vpv=vsv*vpvs*(1+voro(i,4)), vpvs set to 1.73 in
        do i=2,11
            voro(i,1)=30*(i-1) !depth of interface
            voro(i,2)=0.0  !vsv=vsv_prem*(1+voro(i,2))
            voro(i,3)=-1!0.7+0.5/33.*i !xi=voro(i,3), -1 for isotropic layer
            voro(i,4)=0!0.3-0.5/33.*i !vpv=vsv*vpvs*(1+voro(i,4)), vpvs set to 1.73 in params.h
            npt=npt+1
        end do
        
        voro(5,2)=voro(5,2)-0.05
        voro(6,2)=voro(6,2)+0.05
!         do i=5,10
!             voro(i,2)=voro(i,2)-0.05
!         enddo
!         voro(10,4)=0.3
!         voro(11,4)=0.3
!         voro(12,4)=0.3
!         voro(13,4)=0.3
!         voro(14,4)=0.3
!         voro(15,4)=-0.3
!         voro(16,4)=-0.3
!         voro(17,4)=-0.3
!         voro(18,4)=-0.3
!         voro(19,4)=-0.3
        
        
        ! take voro, combine it with prem into a format suitable for minos
        call combine_linear(model_ref,nptref,nic_ref,noc_ref,voro,npt,d_max,&
            r,rho,vpv,vph,vsv,vsh,qkappa,qshear,eta,nptfinal,nic,noc,xi,vpvsv_data)
        
        !calculate synthetic dispersion curves
        if (ndatad_R>0) then
            jcom=3 !rayleigh waves
            
            do iharm=1,nharm_R
                nmodes=0
                call minos_bran(1,tref,nptfinal,nic,noc,r,rho,vpv,vph,vsv,vsh,&
                    qkappa,qshear,eta,eps,wgrav,jcom,lmin,lmax,&
                    wmin_R(iharm),wmax_R(iharm),numharm_R(iharm),numharm_R(iharm),&
                    nmodes_max,nmodes,n_mode,l_mode,c_ph,period,raylquo,error_flag) ! calculate normal modes
                if (error_flag) stop "INVALID INITIAL MODEL - RAYLEIGH - minos_bran.f FAIL 001"
                peri_R_tmp=peri_R(nlims_R(1,iharm):nlims_R(2,iharm))
                ndatad_R_tmp=nlims_R(2,iharm)-nlims_R(1,iharm)+1 ! fortran slices take the first and the last element
                n_R_tmp(:ndatad_R_tmp)=n_R(nlims_R(1,iharm):nlims_R(2,iharm))
                
                call dispersion_minos(nmodes_max,nmodes,n_mode,c_ph,period,raylquo,&
                    peri_R_tmp,n_R_tmp,d_cR_tmp,rq_R,ndatad_R_tmp,ier) ! extract phase velocities from minos output (pretty ugly)
                d_cR(nlims_R(1,iharm):nlims_R(2,iharm))=d_cR_tmp
                if (ier) stop "INVALID INITIAL MODEL - RAYLEIGH - CHANGE PERIODS or MODEL"
                if (maxval(abs(rq_R(:ndatad_R_tmp)))>maxrq*eps) stop "INVALID INITIAL MODEL - RAYLEIGH - CHANGE PERIODS or MODEL"
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
                if (error_flag) stop "INVALID INITIAL MODEL - LOVE - minos_bran.f FAIL 002"
                peri_L_tmp=peri_L(nlims_L(1,iharm):nlims_L(2,iharm))
                ndatad_L_tmp=nlims_L(2,iharm)-nlims_L(1,iharm)+1 ! fortran slices take the first and the last element
                n_L_tmp(:ndatad_L_tmp)=n_L(nlims_L(1,iharm):nlims_L(2,iharm))
                
                call dispersion_minos(nmodes_max,nmodes,n_mode,c_ph,period,raylquo,&
                    peri_L_tmp,n_L_tmp,d_cL_tmp,rq_L,ndatad_L_tmp,ier) ! extract phase velocities from minos output (pretty ugly)
                d_cL(nlims_L(1,iharm):nlims_L(2,iharm))=d_cL_tmp
                if (ier) stop "ier INVALID INITIAL MODEL - LOVE - CHANGE PERIODS or MODEL"
                if (maxval(abs(rq_L(:ndatad_L_tmp)))>maxrq*eps) stop "rq INVALID INITIAL MODEL - LOVE - CHANGE PERIODS or MODEL"
            enddo
        endif

        d_obsdcR_alt(:ndatad_R,2)=d_cR(:ndatad_R)        
        d_obsdcL_alt(:ndatad_L,2)=d_cL(:ndatad_L)
        
        !add errors
        IF (nbproc.gt.1) THEN
            IF (ran==0) THEN
                do i=1,ndatad_R
                    d_obsdcR_alt(i,2)=d_obsdcR_alt(i,2)+gasdev(ra)*0.1
                end do
                do i=1,ndatad_L
                    d_obsdcL_alt(i,2)=d_obsdcL_alt(i,2)+gasdev(ra)*0.1
                end do
                do i=2,nbproc
                    call MPI_SEND(d_obsdcR_alt, ndatadmax*numdis_max, MPI_Real, i-1, 3, MPI_COMM_WORLD, ierror)
                    call MPI_SEND(d_obsdcL_alt, ndatadmax*numdis_max, MPI_Real, i-1, 3, MPI_COMM_WORLD, ierror)
                enddo
            else
                call MPI_RECV(d_obsdcR_alt, ndatadmax*numdis_max, MPI_Real, 0, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierror)
                call MPI_RECV(d_obsdcL_alt, ndatadmax*numdis_max, MPI_Real, 0, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierror)
            endif
        else
            do i=1,ndatad_R
                d_obsdcR_alt(i,2)=d_obsdcR_alt(i,2)+gasdev(ra)*0.1
            end do
            do i=1,ndatad_L
                d_obsdcL_alt(i,2)=d_obsdcL_alt(i,2)+gasdev(ra)*0.1
            end do
        endif
        
        numdis=2
        
        open(65,file=dirname//'/true_model_2.out',status='replace')
        do i=1,nptfinal
            write(65,*)(rearth-r(i))/1000,vsv(i),xi(i),vpvsv_data(i)
        enddo
        close(65)
        
        open(65,file=dirname//'/dispersion.in',status='replace')
        write(65,*)numdis
        write(65,*)lat(:numdis)
        write(65,*)lon(:numdis)
        write(65,*)ndatad_R
        write(65,*)nharm_R
        do k=1,nharm_R
            write(65,*)numharm_R(k)
            write(65,*)nlims_R(2,k)-nlims_R(1,k)
            !write(65,*)wmin_R(k),wmax_R(k)
            do i=nlims_R(1,k),nlims_R(2,k)
                write(65,*)n_R(i),peri_R(i),d_obsdcR(i),d_obsdCRe(i)
                do j=1,numdis
                    write(65,*)d_obsdcR_alt(i,j),d_obsdCRe_alt(i,j)
                enddo
            enddo
        enddo
        
        write(65,*)ndatad_L
        write(65,*)nharm_L
        do k=1,nharm_L
            write(65,*)numharm_L(k)
            write(65,*)nlims_L(2,k)-nlims_L(1,k)
            !write(65,*)wmin_L(k),wmax_L(k)
            do i=nlims_L(1,k),nlims_L(2,k)
                write(65,*)n_L(i),peri_L(i),d_obsdcL(i),d_obsdCLe(i)
                do j=1,numdis
                    write(65,*)d_obsdcL_alt(i,j),d_obsdCLe_alt(i,j)
                enddo
            enddo
        enddo
        close(65)
        
        write(*,*)'DONE INITIALIZING'
    
    else ! real data , untested , unedited, will probably need a little work to get working
        ! GET SWD DATA ---------------------------------------------------------------- 
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
    endif
    
    open(56,file=dirname//'/parameters.out',status='replace')
    write(56,*)'burn-in',burn_in
    write(56,*)'nsample',nsample
    write(56,*)'nsample_widening',nsample_widening
    write(56,*)'burn-in_widening',burn_in_widening
    close(56)
    !***************************************************
    !***************************************************
    
    !***************************************************
    
    !    First loop to find the good widening
    
    !***************************************************
    
    !**************************************************************

    !    Draw the first model randomly from the prior distribution

    !**************************************************************
    
    write(*,*)dirname
    write(1000*rank,*)dirname
    
    widening_prop=widening_start
    
    mean_best=-1000000000
    alphahists=0
    alphahist=0
    
    
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
    stuck=.false.
    
    
    ! to do: branch out on different cores
    
    ! Initilalise sigma (noise parameter)
    Ad_R = Ad_R_max-pAd_R! Ad_min+ran3(ra)*(Ad_R_max-Ad_min)
    Ad_L = Ad_L_max-pAd_L

    

    j=0
    tes=.false.
    do while(.not.tes)
150      tes=.true.
        
        ! create a starting model randomly
        npt = milay+ran3(ra)*(malay-milay)!(12-milay)!(maxlay-minlay)
        !npt=11
        if ((npt>malay).or.(npt<milay)) goto 150
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
        
        call combine_linear(model_ref,nptref,nic_ref,noc_ref,voro,npt,d_max,&
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
                    tes=.false.
                    write(*,*)"Minos_bran FAILED for RAYLEIGH 004"
                end if
                peri_R_tmp=peri_R(nlims_R(1,iharm):nlims_R(2,iharm))
                n_R_tmp=n_R(nlims_R(1,iharm):nlims_R(2,iharm))
                ndatad_R_tmp=nlims_R(2,iharm)-nlims_R(1,iharm)+1 ! fortran slices take the first and the last element
                call dispersion_minos(nmodes_max,nmodes,n_mode,c_ph,period,raylquo,&
                    peri_R_tmp,n_R_tmp,d_cR_tmp,rq_R,ndatad_R_tmp,ier) ! extract phase velocities from minos output (pretty ugly)
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
                end if
                peri_L_tmp=peri_L(nlims_L(1,iharm):nlims_L(2,iharm))
                n_L_tmp=n_L(nlims_L(1,iharm):nlims_L(2,iharm))
                ndatad_L_tmp=nlims_L(2,iharm)-nlims_L(1,iharm)+1 ! fortran slices take the first and the last element
                call dispersion_minos(nmodes_max,nmodes,n_mode,c_ph,period,raylquo,&
                peri_L_tmp,n_L_tmp,d_cL_tmp,rq_L,ndatad_L_tmp,ier) ! extract phase velocities from minos output (pretty ugly)
                d_cL(nlims_L(1,iharm):nlims_L(2,iharm))=d_cL_tmp
                if (ier) tes=.false.
                if (maxval(abs(rq_L(:ndatad_L_tmp)))>maxrq*eps) tes=.false.
            enddo
        
        endif
        
        if (j>250) stop "CAN'T INITIALIZE MODEL" ! if it doesn't work after 250 tries, give up

    enddo

    
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
    like_w=like/widening_prop
    
    if (ran==0) write(*,*)widening_prop
    if (ran==0) write(*,*)like,like_w
    
    
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
    
    
    burnin_in_progress=.true.
    do i_w=1,n_w
        
        
        
        alphamax_prop=-1000000
        alpharefmax_prop=-1000000
        alphamax_props=-1000000
        alpharefmax_props=-1000000
        mean_prop=0
        mean_props=0
        num_prop=0
        num_props=0
        
      
        
      
        ra=rank !seed for RNG
        ran=rank
        inorout=0
        inorouts=0
        
        ier=.false.
        
        call combine_linear(model_ref,nptref,nic_ref,noc_ref,voro,npt,d_max,&
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
                    tes=.false.
                    write(*,*)"Minos_bran FAILED for RAYLEIGH 004"
                end if
                peri_R_tmp=peri_R(nlims_R(1,iharm):nlims_R(2,iharm))
                n_R_tmp=n_R(nlims_R(1,iharm):nlims_R(2,iharm))
                ndatad_R_tmp=nlims_R(2,iharm)-nlims_R(1,iharm)+1 ! fortran slices take the first and the last element
                call dispersion_minos(nmodes_max,nmodes,n_mode,c_ph,period,raylquo,&
                    peri_R_tmp,n_R_tmp,d_cR_tmp,rq_R,ndatad_R_tmp,ier) ! extract phase velocities from minos output (pretty ugly)
                d_cR(nlims_R(1,iharm):nlims_R(2,iharm))=d_cR_tmp
                if (ier) then 
                    tes=.false.
                    write(*,*)'Disperdion_minos RAYLEIGH error'
                endif
                if (maxval(abs(rq_R(:ndatad_R_tmp)))>maxrq*eps) then
                    tes=.false.
                    write(*,*)'raylquo RAYLEIGH error'
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
                    tes=.false.
                    write(*,*)"Minos_bran FAILED for LOVE 004"
                end if
                peri_L_tmp=peri_L(nlims_L(1,iharm):nlims_L(2,iharm))
                n_L_tmp=n_L(nlims_L(1,iharm):nlims_L(2,iharm))
                ndatad_L_tmp=nlims_L(2,iharm)-nlims_L(1,iharm)+1 ! fortran slices take the first and the last element
                call dispersion_minos(nmodes_max,nmodes,n_mode,c_ph,period,raylquo,&
                peri_L_tmp,n_L_tmp,d_cL_tmp,rq_L,ndatad_L_tmp,ier) ! extract phase velocities from minos output (pretty ugly)
                d_cL(nlims_L(1,iharm):nlims_L(2,iharm))=d_cL_tmp
                if (ier) then 
                    tes=.false.
                    write(*,*)'Disperdion_minos LOVE error'
                endif
                if (maxval(abs(rq_L(:ndatad_L_tmp)))>maxrq*eps) then
                    tes=.false.
                    write(*,*)'raylquo LOVE error'
                endif
                
            enddo
        
        endif
        
        if (.not.tes) stop 'Transition model failed ?????'
        
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
        like_w=like/widening_prop
        
        if (ran==0) write(*,*)widening_prop
        if (ran==0) write(*,*)like,like_w
        
        
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
           
                ! UPDATE SIGMAV
                if ((abs((AcB/PrB)-Ar_birth_old)>2*abs((AcB/PrB)-(AcD/PrD))).and.&
                    (AcB.ne.0)) then !special treatement for adding/removing layers
    
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
                
!                 if ((Ac_vp/(Pr_vp+1))<0.001).and. then 
!                     write(filenamemax,"('/stuck_',I3.3,'.out')") rank    
!                     open(56,file=dirname//filenamemax,status='replace')
!                     write(56,*) npt
!                     do i=1,npt
!                         write(56,*)voro(i,1),voro(i,2),voro(i,3),voro(i,4)
!                     enddo
!                     close(56)
!                 endif
                
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
                Prxi=Prxi+1
                ani=.true.
                if (npt_ani.ne.0) then
                    ani=.true.
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
                    
                     !increase counter to calculate acceptance rates
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
     
                PrP=PrP+1
                voro_prop(ind,1)=voro(ind,1)+gasdev(ra)*pd
            
                if ((voro_prop(ind,1)<=d_min).or.(voro_prop(ind,1)>=d_max)) then
                    out=0
                endif
                if ((voro_prop(ind,2)<=-width).or.(voro_prop(ind,2)>=width)) then
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

                if (ind>npt) then
                    write(*,*)npt,ind
                    stop "763"
                endif
                
                PrV=PrV+1
                voro_prop(ind,2)=voro(ind,2)+gasdev(ra)*pv
                
                !voro_prop(ind,2) = -width+2*width*ran3(ra)
                
                !Check if oustide bounds of prior
                if ((voro_prop(ind,2)<=-width).or.(voro_prop(ind,2)>=width)) then
                    out=0
                endif
    
    
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
            elseif (u<1.1) then !death of an anisotropic layer!---------------------------------------    
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
                call combine_linear(model_ref,nptref,nic_ref,noc_ref,voro_prop,npt_prop,d_max,&
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
                            goto 11142
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
                            goto 11142
                        endif
                        
                        if (maxval(abs(rq_R(:ndatad_R_tmp)))>maxrq*eps) then
                            out=0
                            !write(*,*)"BAD UNDERTAINTIES - RAYLEIGH..."
                            goto 11142
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
                            goto 11142
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
                            goto 11142
                        endif
                        
                        if (maxval(abs(rq_L(:ndatad_L_tmp)))>maxrq*eps) then
                            out=0
                            !write(*,*)"BAD UNDERTAINTIES - LOVE..."
                            goto 11142
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
            like_prop_w=like_prop/widening_prop
            
       11142 Accept = .false.
            
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
                logrsig=ndatad_R*log(Ad_R/Ad_R_prop)/widening_prop ! ATTENTION avc ld 2
                if (log(ran3(ra))<logrsig+log(out)-like_prop_w+like_w) then ! hierarchical case
                    accept=.true.
                    Acnd_R=Acnd_R+1
                endif
            elseif (noisd_L) then !@@@@@@@@@@@@@@@@@@@@@@@@@@@@ logrsig  @@@@@@@@@@@@
                logrsig=ndatad_L*log(Ad_L/Ad_L_prop)/widening_prop ! ATTENTION avc ld 2
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
                    
                    logalpha(idis)=ndatad_R*(1/widening_prop-1)*log(Ad_R)&
                    +ndatad_L*(1/widening_prop-1)*log(Ad_L)&
                    +like_w-like_alt
                enddo
                !write(*,*)'logcratio',logcratio
                !write(*,*)'logalpha',logalpha
                
                logalpharef=ndatad_R*(1/widening_prop-1)*log(Ad_R)&
                +ndatad_L*(1/widening_prop-1)*log(Ad_L)&
                +like_w-like
                
        
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
                
                
                IF ((mod(ount,display).EQ.0)) THEN !.and.(mod(ran,10).EQ.0)

                    write(*,*)'processor number',ran+1,'/',nbproc
                    write(*,*)'widening step:','burn-in'
                    write(*,*)'sample:',ount,'/',burn_in
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
                    write(*,*)'widening',widening_prop
                    !write(*,*)Ar_birth_old,sigmav_old,sigmav
                    write(*,*)'-----------------------------------------'
                    write(*,*)
                    write(*,*)
                    
                    
                END IF
                
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
                
                if (ount>burn_in) then
                    burnin_in_progress=.false.
                    ount=0
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
    
            IF (ount.GT.burn_in_widening) THEN
                sample=sample+1
                
                do idis=1,numdis
                    if (logalpha(idis)>alphamax_prop(idis)) then
                        alphamax_prop(idis)=logalpha(idis)
                    endif
                enddo
                
                
                
                if (logalpharef>alpharefmax_prop) then
                    alpharefmax_prop=logalpharef
                endif
                
                IF (mod(ount,thin)==0) THEN
                    
                    !mean_prop=mean_prop+logalpha
                    th = th + 1
                    
                    alpha=exp(logalpha)
                    alphasum=alphasum+exp(logalpha)
                    
                    
                    
                    alphaall(th,:)=logalpha
                    alpharefall(th)=logalpharef
!                     do idis=1,numdis
!                         i_al=ceiling((logalpha(idis)-logalpha_min)/(logalpha_max-logalpha_min)*num_logalpha)
!                         if (i_al<logalpha_min) i_al=1
!                         if (i_al>logalpha_max) i_al=num_logalpha
!                         !alphahist(i_al,idis,i_w)=alphahist(i_al,idis,i_w)+1
!                     enddo

                endif
    
            endif
            
            
!             if (burnin_in_progress) then
!                 convBa(ount)=100*AcBa/PrBa ! acceptance rates in percent
!                 convB(ount)=100*AcB/PrB
!                 convDa(ount)=100*AcDa/PrDa
!                 convD(ount)=100*AcD/PrD
!                 convvp(ount)=100*Ac_vp/Pr_vp
!                 convvs1(ount)=100*Acv(1)/Prv(1)
!                 convvs2(ount)=100*Acv(2)/Prv(2)
!                 convdp1(ount)=100*Acp(1)/Prp(1)
!                 convdp2(ount)=100*Acp(2)/Prp(2)
!                 convxi(ount)=100*Acxi/Prxi
!                 convd_R(ount)=lsd_R
!                 convd_L(ount)=lsd_L
!                 ncell(ount)=npt
!                 convAd_R(ount)=Ad_R
!                 convAd_L(ount)=Ad_L
!             else
            convBa(burn_in+ount+(nsample_widening+burn_in_widening)*(i_w-1))=100*AcBa/PrBa ! acceptance rates in percent
            convB(burn_in+ount+(nsample_widening+burn_in_widening)*(i_w-1))=100*AcB/PrB
            convDa(burn_in+ount+(nsample_widening+burn_in_widening)*(i_w-1))=100*AcDa/PrDa
            convD(burn_in+ount+(nsample_widening+burn_in_widening)*(i_w-1))=100*AcD/PrD
            convvp(burn_in+ount+(nsample_widening+burn_in_widening)*(i_w-1))=100*Ac_vp/Pr_vp
            convvs(burn_in+ount+(nsample_widening+burn_in_widening)*(i_w-1))=100*Acv/Prv
            convdp(burn_in+ount+(nsample_widening+burn_in_widening)*(i_w-1))=100*Acp/Prp
            convxi(burn_in+ount+(nsample_widening+burn_in_widening)*(i_w-1))=100*Acxi/Prxi
            convd_R(burn_in+ount+(nsample_widening+burn_in_widening)*(i_w-1))=lsd_R
            convd_L(burn_in+ount+(nsample_widening+burn_in_widening)*(i_w-1))=lsd_L
            ncell(burn_in+ount+(nsample_widening+burn_in_widening)*(i_w-1))=npt
            convAd_R(burn_in+ount+(nsample_widening+burn_in_widening)*(i_w-1))=Ad_R
            convAd_L(burn_in+ount+(nsample_widening+burn_in_widening)*(i_w-1))=Ad_L
!             endif
            
            IF ((mod(ount,display).EQ.0)) THEN !.and.(mod(ran,10).EQ.0)

                write(*,*)'processor number',ran+1,'/',nbproc
                write(*,*)'widening step:',i_w,'/',n_w
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
                write(*,*)'widening',widening_prop
                !write(*,*)Ar_birth_old,sigmav_old,sigmav
                write(*,*)'-----------------------------------------'
                write(*,*)
                write(*,*)
                
                
            END IF
        end do !End Markov chain
        
        !write(*,*)'alphamax1',alphamax_prop(:numdis)
        
        !alphamax_prop=maxval(alphaall,1)
        
        !write(*,*)'alphamax2',alphamax_prop(:numdis)
        
        k=0
        if (th.ne.0) then !normalize averages
            
            !mean_prop=mean_prop/th
            
            inorout(ran+1)=1
            k=k+1

        endif



        do i=1,nbproc
            call MPI_REDUCE(inorout,inorouts,21000,MPI_Integer,MPI_Sum,i-1,MPI_COMM_WORLD,ierror)
        enddo
        

        !write(*,*)'rank=',ran,'th=',th

        j=0
        k=0
        do i=1,nbproc
            j=j+inorouts(i)
            if (inorouts(i).ne.0) then
                k=k+1
                members(k)=i-1
            endif
        enddo

        !IF (ran==0) write(*,*) 'k',k,'nbproc',nbproc
        !***************************************************************************

        ! Collect information from all the chains and average everything

        !***************************************************************************
        flag=0
        do i=1,j
            if (ran==members(i)) flag=1
        enddo

        call MPI_Group_incl(group_world, j, members, good_group, ierror)
        call MPI_Comm_create(MPI_COMM_WORLD, good_group, MPI_COMM_small, ierror)
        

        if (flag==1) then
            
            ! reduce alphamax to get the max over the cores, send it to each core
            
            call MPI_ALLREDUCE(alphamax_prop,alphamax_props,numdis_max,MPI_Real,MPI_MAX,MPI_COMM_small,ierror)
            call MPI_ALLREDUCE(alpharefmax_prop,alpharefmax_props,1,MPI_Real,MPI_MAX,MPI_COMM_small,ierror)
            
            ! make alpha histogram
            
            do idis=1,numdis
                
                do i=1,th
                    i_al=ceiling(((alphaall(i,idis)-alphamax_props(idis))-logalpha_min)/(logalpha_max-logalpha_min)*num_logalpha)
                    if (i_al<1) i_al=1
                    if (i_al>num_logalpha) i_al=num_logalpha
                    alphahist(i_al,idis,i_w)=alphahist(i_al,idis,i_w)+1
                    
                enddo
            enddo
            
            do i=1,th
                i_al=ceiling(((alpharefall(i)-alpharefmax_props)-logalpha_min)/(logalpha_max-logalpha_min)*num_logalpha)
                if (i_al<1) i_al=1
                if (i_al>num_logalpha) i_al=num_logalpha
                alpharefhist(i_al,i_w)=alpharefhist(i_al,i_w)+1
            enddo
            
            call MPI_ALLREDUCE(alphahist,alphahists,num_logalpha*numdis_max*n_w,MPI_Real,MPI_SUM,&
                 MPI_COMM_small,ierror)
            
            call MPI_ALLREDUCE(alpharefhist,alpharefhists,num_logalpha*n_w,MPI_Real,MPI_SUM,&
                 MPI_COMM_small,ierror)
            
            ! get new alphamax as to cap at 1% to whatever is over -15
            
            alphamax_true(i_w,:)=alphamax_props
            alpharefmax_true(i_w)=alpharefmax_props
            
            i_opt=int((-15-logalpha_min)/(logalpha_max-logalpha_min)*num_logalpha)
            
            do idis=1,numdis
                
                sum_tmp=0
                isum=num_logalpha+1
                
                !i_opt=isum-int((15)/(logalpha_max-logalpha_min)*num_logalpha)
                i_opt=1
                totsum=sum(alphahists(i_opt:,idis,i_w))
                do while ((sum_tmp<0.01*totsum).or.(totsum==0))
                    isum=isum-1
                    sum_tmp=sum_tmp+alphahists(isum,idis,i_w)
                    !i_opt=i_opt-1
                    !totsum=sum(alphahists(i_opt:,idis,i_w))
                enddo
                alphamax_props(idis)=alphamax_props(idis)+(logalpha_min+isum*((logalpha_max-logalpha_min)/num_logalpha))
            enddo
            
            sum_tmp=0
            isum=num_logalpha+1
            !i_opt=isum-int((15)/(logalpha_max-logalpha_min)*num_logalpha)
            i_opt=1
            totsum=sum(alpharefhists(i_opt:,i_w))
            do while ((sum_tmp<0.01*totsum).or.(totsum==0))
                isum=isum-1
                sum_tmp=sum_tmp+alpharefhists(isum,i_w)
                !i_opt=i_opt-1
                !totsum=sum(alpharefhists(i_opt:,i_w))
            enddo
            alpharefmax_props=alpharefmax_props+(logalpha_min+isum*((logalpha_max-logalpha_min)/num_logalpha))
            
            do idis=1,numdis
                
                do i=1,th
                    !if ((alphaall(i,idis)-alphamax_props(idis))>-15) then
                        mean_prop(idis)=mean_prop(idis)+min(0.,alphaall(i,idis)-alphamax_props(idis))
                        num_prop(idis)=num_prop(idis)+1
                    !endif
                enddo
            enddo
            
            call MPI_ALLREDUCE(mean_prop,mean_props,numdis_max,MPI_Real,MPI_sum,MPI_COMM_small,ierror)
            call MPI_ALLREDUCE(num_prop,num_props,numdis_max,MPI_Real,MPI_sum,MPI_COMM_small,ierror)
            
            call MPI_Group_free(good_group, ierror)
            call MPI_Comm_free(MPI_COMM_small, ierror)
            call MPI_BARRIER(MPI_COMM_WORLD, ierror)
            
            mean_props=mean_props/num_props!-alphamax_props
            meandiff_hist(i_w,:)=mean_props
            alphamax_hist(i_w,:)=alphamax_props
            alphamaxref_hist(i_w)=alpharefmax_props
            if (minval(mean_props(:numdis))>minval(mean_best(:numdis))) then
                mean_best=mean_props
                widening=widening_prop
                alphamax=alphamax_props
                alpharefmax=alpharefmax_props
                
                best_current_Ad_R=Ad_R
                best_current_Ad_L=Ad_L
                best_current_pxi=pxi
                best_current_p_vp=p_vp
                best_current_pAd_R=pAd_R
                best_current_pAd_L=pAd_L
                best_current_pv=pv
                best_current_pd=pd
                best_current_sigmav=sigmav
                best_current_npt=npt
                best_current_voro=voro
                
            endif
            IF (ran==members(1)) THEN
                write(*,*)'New widening tested: '
                write(*,*)i_w,widening_prop,widening,'best',mean_best(:numdis),'props',mean_props(:numdis)
                write(*,*)'alphamax',alphamax_props(:numdis)
            endif
            widening_prop=widening_prop+widening_step
            
            
        endif
    enddo
    
    IF (ran==members(1)) THEN
        write(*,*)'Best widening: '
        write(*,*)'widening',widening
        write(*,*)'mean',mean_best(:numdis)
        write(*,*)'alphamax',alphamax(:numdis)
    endif
    
    do i=1,nbproc
        call MPI_REDUCE(inorout,inorouts,21000,MPI_Integer,MPI_Sum,i-1,MPI_COMM_WORLD,ierror)
    enddo

    !write(*,*)'rank=',ran,'th=',th

    j=0
    k=0
    do i=1,nbproc
        j=j+inorouts(i)
        if (inorouts(i).ne.0) then
            k=k+1
            members(k)=i-1
        endif
    enddo

    !IF (ran==0) write(*,*) 'k',k,'nbproc',nbproc
    !***************************************************************************

    ! Collect information from all the chains and average everything

    !***************************************************************************
    flag=0
    do i=1,j
        if (ran==members(i)) flag=1
    enddo
    
    !write(*,*)'alphamax_prop',alphamax_prop


    call MPI_Group_incl(group_world, j, members, good_group, ierror)
    call MPI_Comm_create(MPI_COMM_WORLD, good_group, MPI_COMM_small, ierror)
    
    if (flag==1) then
        
        call MPI_REDUCE(convBa,convBas,n_w*(burn_in_widening+nsample_widening)+burn_in,&
        MPI_Real,MPI_Sum,0,MPI_COMM_small,ierror)
        call MPI_REDUCE(convB,convBs,n_w*(burn_in_widening+nsample_widening)+burn_in,&
        MPI_Real,MPI_Sum,0,MPI_COMM_small,ierror)
        call MPI_REDUCE(convDa,convDas,n_w*(burn_in_widening+nsample_widening)+burn_in,&
        MPI_Real,MPI_Sum,0,MPI_COMM_small,ierror)
        call MPI_REDUCE(convD,convDs,n_w*(burn_in_widening+nsample_widening)+burn_in,&
        MPI_Real,MPI_Sum,0,MPI_COMM_small,ierror)
        call MPI_REDUCE(convvp,convvps,n_w*(burn_in_widening+nsample_widening)+burn_in,&
        MPI_Real,MPI_Sum,0,MPI_COMM_small,ierror)
        call MPI_REDUCE(convvs,convvss,n_w*(burn_in_widening+nsample_widening)+burn_in,&
        MPI_Real,MPI_Sum,0,MPI_COMM_small,ierror)
        call MPI_REDUCE(convdp,convdps,n_w*(burn_in_widening+nsample_widening)+burn_in,&
        MPI_Real,MPI_Sum,0,MPI_COMM_small,ierror)
        call MPI_REDUCE(convxi,convxis,n_w*(burn_in_widening+nsample_widening)+burn_in,&
        MPI_Real,MPI_Sum,0,MPI_COMM_small,ierror)
        call MPI_REDUCE(convd_R,convd_Rs,n_w*(burn_in_widening+nsample_widening)+burn_in,&
        MPI_Real,MPI_Sum,0,MPI_COMM_small,ierror)
        call MPI_REDUCE(convd_L,convd_Ls,n_w*(burn_in_widening+nsample_widening)+burn_in,&
        MPI_Real,MPI_Sum,0,MPI_COMM_small,ierror)
        call MPI_REDUCE(convAd_R,convAd_Rs,n_w*(burn_in_widening+nsample_widening)+burn_in,&
        MPI_Real,MPI_Sum,0,MPI_COMM_small,ierror)
        call MPI_REDUCE(convAd_L,convAd_Ls,n_w*(burn_in_widening+nsample_widening)+burn_in,&
        MPI_Real,MPI_Sum,0,MPI_COMM_small,ierror)
        call MPI_REDUCE(ncell,ncells,n_w*(burn_in_widening+nsample_widening)+burn_in,&
        MPI_Real,MPI_Sum,0,MPI_COMM_small,ierror)

        
        convPs=convPs/j
        convBs=convBs/j
        convBas=convBas/j
        convDs=convDs/j
        convDas=convDas/j
        convvps=convvps/j
        convvss=convvss/j
        convdps=convdps/j
        convxis=convxis/j
        convd_Rs=convd_Rs/j
        convd_Ls=convd_Ls/j
        convAd_Rs=convAd_Rs/j
        convAd_Ls=convAd_Ls/j
        ncells=ncells/j
    endif
    
    call MPI_Group_free(good_group, ierror)
    call MPI_Comm_free(MPI_COMM_small, ierror)
    call MPI_BARRIER(MPI_COMM_WORLD, ierror)
    
    write(*,*)'writing outputs of widening tests'
    
    
    write(filenamemax,"('/last_model_',I3.3,'.inout')") rank    
    open(65,file=dirname//filenamemax,status='replace')
    write(65,*)best_current_Ad_R
    write(65,*)best_current_Ad_L
    write(65,*)best_current_pxi
    write(65,*)best_current_p_vp
    write(65,*)best_current_pAd_R
    write(65,*)best_current_pAd_L
    write(65,*)best_current_pv
    write(65,*)best_current_pd
    write(65,*)best_current_sigmav
    write(65,*)best_current_npt
    do i=1,best_current_npt
        write(65,*)best_current_voro(i,:)
    enddo
    close(65)
    
    IF (ran==members(1)) THEN

        open(65,file=dirname//'/mean_prop.out',status='replace')
        write(65,*)widening_start,widening_step,n_w,numdis
        do i=1,n_w
            write(65,*)widening_start+(i-1)*widening_step,alphamaxref_hist(i),&
            alpharefmax_true(i)
            write(65,*)meandiff_hist(i,:numdis)
            write(65,*)alphamax_hist(i,:numdis)
            write(65,*)alphamax_true(i,:numdis)
        enddo
        close(65)
        
        open(65,file=dirname//'/alphahist.out',status='replace')
        write(65,*)logalpha_min,logalpha_max,num_logalpha
        write(65,*)widening_start,widening_step,n_w,numdis
        do i=1,n_w
            do j=1,numdis
                write(65,*)alphahists(:,j,i)
            enddo
        enddo
        close(65)
        
        open(65,file=dirname//'/goodalpha.inout',status='replace')
        write(65,*)widening,alpharefmax
        write(65,*)numdis
        do i=1,numdis
            write(65,*)alphamax(i)
        enddo
        write(65,*)mean_best(:numdis)
        close(65)
        
        !write(*,*)'final alphamax',alpharefmax,alphamax(:numdis)
        
        open(54,file=dirname//'/Convergence_Birth_prep.out',status='replace')
        write(54,*)burn_in,n_w,burn_in_widening,nsample_widening
        do i=1,n_w*(burn_in_widening+nsample_widening)+burn_in
            write(54,*)convB(i),convBs(i),convBa(i),convBas(i)
        enddo
        close(54)
        
        open(54,file=dirname//'/Convergence_Death_prep.out',status='replace')
        write(54,*)burn_in,n_w,burn_in_widening,nsample_widening
        do i=1,n_w*(burn_in_widening+nsample_widening)+burn_in
            write(54,*)convD(i),convDs(i),convDa(i),convDas(i)
        enddo
        close(54)
        
        open(54,file=dirname//'/Convergence_vp_prep.out',status='replace')
        write(54,*)burn_in,n_w,burn_in_widening,nsample_widening
        do i=1,n_w*(burn_in_widening+nsample_widening)+burn_in
            write(54,*)convvp(i),convvps(i)
        enddo
        close(54)
        
        open(54,file=dirname//'/Convergence_vs_prep.out',status='replace')
        write(54,*)burn_in,n_w,burn_in_widening,nsample_widening
        do i=1,n_w*(burn_in_widening+nsample_widening)+burn_in
            write(54,*)convvs(i),convvss(i)
        enddo
        close(54)
        
        open(54,file=dirname//'/Convergence_dp_prep.out',status='replace')
        write(54,*)burn_in,n_w,burn_in_widening,nsample_widening
        do i=1,n_w*(burn_in_widening+nsample_widening)+burn_in
            write(54,*)convdp(i),convdps(i)
        enddo
        close(54)
        
        open(54,file=dirname//'/Convergence_xi_prep.out',status='replace')
        write(54,*)burn_in,n_w,burn_in_widening,nsample_widening
        do i=1,n_w*(burn_in_widening+nsample_widening)+burn_in
            write(54,*)convxi(i),convxis(i)
        enddo
        close(54)
        
        open(54,file=dirname//'/Convergence_misfit_prep.out',status='replace')
        write(54,*)burn_in,n_w,burn_in_widening,nsample_widening
        do i=1,n_w*(burn_in_widening+nsample_widening)+burn_in
            write(54,*)convd_R(i),convd_Rs(i),convd_L(i),convd_Ls(i)
        enddo
        close(54)

        open(53,file=dirname//'/Convergence_nb_layers_prep.out',status='replace')
        write(53,*)burn_in,n_w,burn_in_widening,nsample_widening
        do i=1,n_w*(burn_in_widening+nsample_widening)+burn_in
            write(53,*)ncell(i),ncells(i)
        enddo
        close(53)

        open(52,file=dirname//'/Convergence_sigma_R_prep.out',status='replace')
        write(52,*)burn_in,n_w,burn_in_widening,nsample_widening
        do i=1,n_w*(burn_in_widening+nsample_widening)+burn_in
            write(52,*)convAd_R(i),convAd_Rs(i)
        enddo
        close(52)

        open(53,file=dirname//'/Convergence_sigma_L_prep.out',status='replace')
        write(53,*)burn_in,n_w,burn_in_widening,nsample_widening
        do i=1,n_w*(burn_in_widening+nsample_widening)+burn_in
            write(53,*)convAd_L(i),convAd_Ls(i)
        enddo
        close(53)
    ENDIF
    
    CALL cpu_time(t2)
    if (ran==0) write(*,*)'time taken by the code was',t2-t1,'seconds'
    if (ran==0) write(*,*)'time taken by the code was',(t2-t1)/3600,'hours'
    
    close(outputfile)
    
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
