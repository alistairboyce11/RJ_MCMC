program postprocess_binary_outputs
    use mpi
    implicit none

    include 'params.h'
    include 'data_joint.h'

    real widening_prop
    real widening_prop2
    integer i_w,i,io_file,io
    character filenamemax*300,filenamemax2*300
    integer numsample,everyall_2,burn_in_2,thin_2
    real d_min_2,d_max_2,width_2,xi_min_2,xi_max_2,vp_min_2,vp_max_2
    DOUBLE PRECISION Ad_R_min_2,Ad_R_max_2,Ad_L_min_2,Ad_L_max_2
    integer milay_2,malay_2,mk_2,ndatadmax_2
    integer nptfinal,npt,npt_ani,nic,noc
    DOUBLE PRECISION Ad_R,Ad_L
    real, DIMENSION(mk) :: r,vsv,xi,vpv
    real like_w
    integer ndatad_R,ndatad_L
    real, DIMENSION(ndatadmax) :: d_cR,d_cL

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Read binary files containing models and converts them into ascii for python to read
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    widening_prop2=widening_start
    do i_w=1,n_w
        write(*,*)widening_prop2
        do i=0,50
            write(filenamemax,"('/All_models_prepare_',I4.4,'_',f5.2,'.out')")i,widening_prop2
            write(*,*)dirname//filenamemax
            open(100,file=dirname//filenamemax,status='old',form='unformatted',access='stream',iostat=io_file)

            if (io_file/=0) goto 500
            write(filenamemax2,"('/All_models_processed_prepare_',I4.4,'_',f5.2,'.out')")i,widening_prop2
            write(*,*)filenamemax2
            open(200,file=dirname//filenamemax2,status='replace')

            read(100,IOSTAT=io)numsample

            read(100,IOSTAT=io)everyall_2

            write(200,*)numsample,everyall_2
            write(*,*)numsample,everyall_2

            read(100,IOSTAT=io)burn_in_2
            read(100,IOSTAT=io)widening_prop
            read(100,IOSTAT=io)thin_2

            write(200,*)burn_in_2,widening_prop,thin_2

            read(100,IOSTAT=io)d_min_2
            read(100,IOSTAT=io)d_max_2

            write(200,*)d_min_2,d_max_2

            read(100,IOSTAT=io)width_2

            write(200,*)width_2

            read(100,IOSTAT=io)xi_min_2
            read(100,IOSTAT=io)xi_max_2

            write(200,*)xi_min_2,xi_max_2

            read(100,IOSTAT=io)vp_min_2
            read(100,IOSTAT=io)vp_max_2

            write(200,*)vp_min_2,vp_max_2

            read(100,IOSTAT=io)Ad_R_min_2
            read(100,IOSTAT=io)Ad_R_max_2

            write(200,*)Ad_R_min_2,Ad_R_max_2

            read(100,IOSTAT=io)Ad_L_min_2
            read(100,IOSTAT=io)Ad_L_max_2

            write(200,*)Ad_L_min_2,Ad_L_max_2

            read(100,IOSTAT=io)milay_2
            read(100,IOSTAT=io)malay_2

            write(200,*)milay_2,malay_2

            read(100,IOSTAT=io)mk_2
            read(100,IOSTAT=io)ndatadmax_2

            io=0

            do while (io==0)

                read(100,IOSTAT=io)nptfinal
                !write(*,*)nptfinal
                if (io/=0) goto 500
                if (nptfinal>mk) CONTINUE

                read(100,IOSTAT=io)nic
                if (io/=0) goto 500

                read(100,IOSTAT=io)noc
                if (io/=0) goto 500

                read(100,IOSTAT=io)npt
                if (io/=0) goto 500

                read(100,IOSTAT=io)npt_ani
                if (io/=0) goto 500

                write(200,*)nptfinal-noc,npt,npt_ani

                read(100,IOSTAT=io)Ad_R
                if (io/=0) goto 500

                read(100,IOSTAT=io)Ad_L
                if (io/=0) goto 500

                write(200,*)Ad_R,Ad_L

                read(100,IOSTAT=io)r
                if (io/=0) goto 500

                write(200,*)(rearth-r(nptfinal:noc+1:-1))/1000. ! only take what is needed for postprocessing

                read(100,IOSTAT=io)vsv
                if (io/=0) goto 500

                write(200,*)vsv(nptfinal:noc+1:-1)/1000.

                read(100,IOSTAT=io)xi
                if (io/=0) goto 500

                write(200,*)xi(nptfinal:noc+1:-1)

                read(100,IOSTAT=io)vpv ! careful! We store vpv, so we need to check it for the python postprocessing
                if (io/=0) goto 500

                write(200,*)vpv(nptfinal:noc+1:-1)/1000.

                read(100,IOSTAT=io)like_w
                if (io/=0) goto 500

                read(100,IOSTAT=io)ndatad_R
                if (io/=0) goto 500

                read(100,IOSTAT=io)d_cR
                if (io/=0) goto 500

                read(100,IOSTAT=io)ndatad_L
                if (io/=0) goto 500

                read(100,IOSTAT=io)d_cL
                if (io/=0) goto 500

                ! because it wasn't put in the initial like_w
                like_w=-like_w-ndatad_R*log(Ad_R)/widening_prop-ndatad_L*log(Ad_L)/widening_prop

                write(200,*)like_w

                write(200,*)ndatad_R

                write(200,*)d_cR(:ndatad_R)

                write(200,*)ndatad_L

                write(200,*)d_cL(:ndatad_L)
            enddo


        600 close(100)
        500 close(200)
        enddo
        widening_prop2=widening_prop2+widening_step
    enddo

end program postprocess_binary_outputs
