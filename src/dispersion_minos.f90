subroutine dispersion_minos(nmodes_max,nmodes,n_mode,c_ph,period,raylquo,peri,n,d_c,rq,ndatad,ier,idis,imin_idis)

    implicit none
    include '../params.h'
    integer,intent(in) :: nmodes_max,nmodes,ndatad
    integer,dimension(nmodes_max),intent(in) :: n_mode
    real,dimension(nmodes_max),intent(in) :: period,c_ph,raylquo
    real,dimension(ndatadmax),intent(in) :: peri
    integer,dimension(ndatadmax),intent(in) :: n
    real,dimension(ndatadmax),intent(out) :: d_c,rq
    integer :: j,i,imin,imin2
    real :: dmin,dmin_new
    logical,intent(out) :: ier,idis
    integer,intent(out) :: imin_idis
    logical :: allfine
    real, external :: interp
    
    imin2=2
    ier=.false.
    idis=.false.
    allfine=.true.
    
    do j=1,ndatad
        if (n(j)==n_mode(imin2-1)) then
            imin=imin2-1
            dmin=abs(peri(j)-period(imin))
        else
            imin=1
            dmin=1000000
        end if
        do i=imin,nmodes
            if (n(j)==n_mode(i)) then
                dmin_new=abs(peri(j)-period(i))
                if (dmin_new>dmin) then
                    if (peri(j)<period(i-1)) then
                        d_c(j)=interp(period(i),period(i-1),c_ph(i),c_ph(i-1),peri(j))
                        rq(j)=interp(period(i),period(i-1),raylquo(i),raylquo(i-1),peri(j))
                    else 
                        if ((i-2>0).and.(n(j)==n_mode(i-2))) then
                            d_c(j)=interp(period(i-2),period(i-1),c_ph(i-2),c_ph(i-1),peri(j))
                            rq(j)=interp(period(i-2),period(i-1),raylquo(i-2),raylquo(i-1),peri(j))
                        else
                            d_c(j)=c_ph(i-1)
                            rq(j)=raylquo(i-1)
                            if (dmin>maxdis) ier=.true.
                        endif
                    endif
        
                    imin2=i
                    exit
                elseif (i==nmodes) then
                    if (peri(j)<period(i)) then
                        d_c(j)=c_ph(i)    
                        rq(j)=raylquo(i)    
                        imin2=i
                        dmin=abs(peri(j)-period(i))
                        if (dmin>maxdis) ier=.true.
                        exit
                    else
                        d_c(j)=interp(period(i),period(i-1),c_ph(i),c_ph(i-1),peri(j))
                        rq(j)=interp(period(i),period(i-1),raylquo(i),raylquo(i-1),peri(j))
                        imin2=i
                        exit
                    endif
                else
                    dmin=dmin_new
                end if
            elseif (n(j)<n_mode(i)) then
                if ((i-2>0).and.(n(j)==n_mode(i-2)).and.(peri(j)>period(i-1))) then
                    d_c(j)=interp(period(i-2),period(i-1),c_ph(i-2),c_ph(i-1),peri(j))
                    rq(j)=interp(period(i-2),period(i-1),raylquo(i-2),raylquo(i-1),peri(j))
                else
                    d_c(j)=c_ph(i-1) 
                    rq(j)=raylquo(i-1) 
                endif                   
                imin2=i
                if (dmin>maxdis) ier=.true.
                exit
            endif
        end do
        
    end do
    
end subroutine dispersion_minos