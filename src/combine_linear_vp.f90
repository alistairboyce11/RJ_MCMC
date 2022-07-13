subroutine combine_linear_vp(model_ref,nptref,nic_ref,noc_ref,voro,npt,d_max,&
    r,rho,vpv,vph,vsv,vsh,qkappa,qshear,eta,nptfinal,nic,noc,vs_data,xi,vp_data)
    implicit none
    include '../params.h'
    
    real, external :: interp
    real,intent(in) :: d_max
    integer,intent(in) :: nic_ref,noc_ref,nptref,npt
    real,dimension(mk,9),intent(in) :: model_ref
    real,dimension(malay,4),intent(in) :: voro
    real,dimension(malay,4) :: voro2
    real,dimension(4) :: t
    real, parameter :: phi = 1 ! P-wave radial anisotropy variable: phi= (vpv/vph)^2
    
    real,dimension(malay) :: vshvsv,vpvsv
    real,dimension(malay+1) :: radius
    integer :: i,j,k,c
    
    real,dimension(mk),intent(out) :: r,rho,vpv,vph,vsv,vsh,qkappa,qshear,eta
    real,dimension(mk),intent(out) :: xi,vp_data,vs_data
    ! real,dimension(mk) :: vs_data
    integer,intent(out) :: nptfinal,nic,noc
    
    ! write(*,*) nptref, nic_ref, noc_ref
    ! do i = 1,nptref
    !      write(*,*) i, model_ref(i,1),&
    !                 model_ref(i,2),model_ref(i,3),model_ref(i,4),model_ref(i,5),&
    !                 model_ref(i,6),model_ref(i,7),model_ref(i,8),model_ref(i,9)
    ! end do


    voro2=voro
    !call quicksort(voro2,npt)
    
    do i=1,npt
        do j=1,npt-1 
            if (voro2(j,1).gt.voro2(j+1,1)) then
                t = voro2(j,:);  voro2(j,:) = voro2(j+1,:);  voro2(j+1,:) = t
            ENDIF
        ENDDO
    ENDDO

    !voro(i,:)=[thickness,vsv,(vsh/vsv)**2,vp/vsv]
    do i=1,npt
        ! radius(i)=rearth-float(int(voro2(i,1)*1000)) ! Round to nearest meter
        radius(i)=rearth-float(int(voro2(i,1))*1000) ! Round to nearest km
    end do
    radius(npt+1)=rearth-d_max*1000
    
    !vshvsv=voro2(:,3)
    do i=1,npt
        vshvsv(i)=voro2(i,3)
        if (vshvsv(i)==-1) vshvsv(i)=1
    end do
    
    ! do i=1,npt
    !     vpvsv(i)=vpvs*(1+voro2(i,4))
    ! end do
    
    i=nptref
    j=1
    k=nptref+2*npt
    nic=nic_ref
    noc=noc_ref
    nptfinal=nptref+2*npt
    c=0 ! Number of discont. in voro that match reference model.

    r(k)=model_ref(i,1)
    rho(k)=model_ref(i,2)
    qkappa(k)=model_ref(i,5)
    qshear(k)=model_ref(i,6)
    eta(k)=model_ref(i,9)
    
    vsv(k)=model_ref(i,4)*(1+voro2(j,2))
    ! vph(k)=vsv(k)*vpvsv(j)
    ! vpv(k)=vph(k)
    vsh(k)=vsv(k)*sqrt(vshvsv(j))

    !!!!!!!!!!!!!!!!! New !!!!!!!!!!!!!!!!!!!!!!!!
    vph(k)=model_ref(i,7)*(1+voro2(j,4))
    vpv(k)=vph(k)*sqrt(phi)


    xi(k)=vshvsv(j)
    vp_data(k)=voro2(j,4)
    vs_data(k)=voro2(j,2)
    
    j=j+1

    k=k-1
    i=i-1
    
    do while ((i>=1).or.(j<npt+2))
        ! write(*,*)j
        ! write(*,*)i
        ! write(*,*)'model_ref',model_ref(i,1)
        ! write(*,*)'radius',radius(j)
        ! write(*,*)j,npt,i,radius(j)
        
        if (j>npt) then
        
            if ((j>npt+1).and.(model_ref(i,1)<radius(npt+1))) then
                r(k)=model_ref(i,1)
                rho(k)=model_ref(i,2)
                vpv(k)=model_ref(i,3)
                vsv(k)=model_ref(i,4)
                qkappa(k)=model_ref(i,5)
                qshear(k)=model_ref(i,6)
                vph(k)=model_ref(i,7)
                vsh(k)=model_ref(i,8)
                eta(k)=model_ref(i,9)
                xi(k)=(vsh(k)/vsv(k))**2
                vp_data(k)=0.0
                vs_data(k)=0.0
                ! write(*,*)r(k), "Here 1"

                k=k-1
                i=i-1
            
            elseif ((j==npt+1).and.(model_ref(i,1)<=radius(j))) then
                r(k)=radius(j)
                rho(k)=interp(model_ref(i,1),model_ref(i+1,1),model_ref(i,2),model_ref(i+1,2),r(k))
                qkappa(k)=interp(model_ref(i,1),model_ref(i+1,1),model_ref(i,5),model_ref(i+1,5),r(k))
                qshear(k)=interp(model_ref(i,1),model_ref(i+1,1),model_ref(i,6),model_ref(i+1,6),r(k))
                eta(k)=interp(model_ref(i,1),model_ref(i+1,1),model_ref(i,9),model_ref(i+1,9),r(k))
                
                vsv(k)=interp(model_ref(i,1),model_ref(i+1,1),model_ref(i,4),model_ref(i+1,4),r(k))*(1+voro2(j-1,2))
                vph(k)=interp(model_ref(i,1),model_ref(i+1,1),model_ref(i,7),model_ref(i+1,7),r(k))*(1+voro2(j-1,4))
                vsh(k)=vsv(k)*sqrt(vshvsv(j-1))
                vpv(k)=vph(k)*sqrt(phi)

                xi(k)=vshvsv(j-1)
                vp_data(k)=voro2(j-1,4)
                ! write(*,*)r(k), "Here 2"
                k=k-1
                !
                
                r(k)=radius(j)
                rho(k)=interp(model_ref(i,1),model_ref(i+1,1),model_ref(i,2),model_ref(i+1,2),r(k))
                qkappa(k)=interp(model_ref(i,1),model_ref(i+1,1),model_ref(i,5),model_ref(i+1,5),r(k))
                qshear(k)=interp(model_ref(i,1),model_ref(i+1,1),model_ref(i,6),model_ref(i+1,6),r(k))
                eta(k)=interp(model_ref(i,1),model_ref(i+1,1),model_ref(i,9),model_ref(i+1,9),r(k))
                
                vsh(k)=interp(model_ref(i,1),model_ref(i+1,1),model_ref(i,8),model_ref(i+1,8),r(k))
                vpv(k)=interp(model_ref(i,1),model_ref(i+1,1),model_ref(i,3),model_ref(i+1,3),r(k))
                vsv(k)=interp(model_ref(i,1),model_ref(i+1,1),model_ref(i,4),model_ref(i+1,4),r(k))
                vph(k)=interp(model_ref(i,1),model_ref(i+1,1),model_ref(i,7),model_ref(i+1,7),r(k))
                
                xi(k)=(vsh(k)/vsv(k))**2
                vp_data(k)=0.0
                vs_data(k)=0.0
                ! write(*,*)r(k), "Here 3"
                if (r(k)<model_ref(nic_ref,1)) nic=nic+2
                if (r(k)<model_ref(noc_ref,1)) noc=noc+2  
                
                k=k-1
                
                j=j+1
                
                
            elseif ((j==npt+1).and.(model_ref(i,1)>radius(j))) then
                r(k)=model_ref(i,1)
                rho(k)=model_ref(i,2)
                qkappa(k)=model_ref(i,5)
                qshear(k)=model_ref(i,6)
                eta(k)=model_ref(i,9)
                
                vsv(k)=model_ref(i,4)*(1+voro2(j-1,2))
                vph(k)=model_ref(i,7)*(1+voro2(j-1,4))
                vsh(k)=vsv(k)*sqrt(vshvsv(j-1))
                vpv(k)=vph(k)*sqrt(phi)

                xi(k)=vshvsv(j-1)
                vp_data(k)=voro2(j-1,4)
                vs_data(k)=voro2(j-1,2)
                ! write(*,*)r(k), "Here 4"
                k=k-1
                i=i-1 
            endif
        else 
            if (model_ref(i,1)>radius(j)) then
                r(k)=model_ref(i,1)
                rho(k)=model_ref(i,2)
                qkappa(k)=model_ref(i,5)
                qshear(k)=model_ref(i,6)
                eta(k)=model_ref(i,9)
                
                vsv(k)=model_ref(i,4)*(1+voro2(j-1,2))
                vph(k)=model_ref(i,7)*(1+voro2(j-1,4))
                vsh(k)=vsv(k)*sqrt(vshvsv(j-1))
                vpv(k)=vph(k)*sqrt(phi)

                xi(k)=vshvsv(j-1)
                vp_data(k)=voro2(j-1,4)
                vs_data(k)=voro2(j-1,2)
                ! write(*,*)r(k), "Here 5"
                k=k-1
                i=i-1
            elseif (model_ref(i,1)<radius(j)) then
                
                r(k)=radius(j)
                rho(k)=interp(model_ref(i,1),model_ref(i+1,1),model_ref(i,2),model_ref(i+1,2),r(k))
                qkappa(k)=interp(model_ref(i,1),model_ref(i+1,1),model_ref(i,5),model_ref(i+1,5),r(k))
                qshear(k)=interp(model_ref(i,1),model_ref(i+1,1),model_ref(i,6),model_ref(i+1,6),r(k))
                eta(k)=interp(model_ref(i,1),model_ref(i+1,1),model_ref(i,9),model_ref(i+1,9),r(k))
                
                vsv(k)=interp(model_ref(i,1),model_ref(i+1,1),model_ref(i,4),model_ref(i+1,4),r(k))*(1+voro2(j-1,2))
                vph(k)=interp(model_ref(i,1),model_ref(i+1,1),model_ref(i,7),model_ref(i+1,7),r(k))*(1+voro2(j-1,4))
                vsh(k)=vsv(k)*sqrt(vshvsv(j-1))
                vpv(k)=vph(k)*sqrt(phi)
                xi(k)=vshvsv(j-1)
                vp_data(k)=voro2(j-1,4)
                vs_data(k)=voro2(j-1,2)
                ! write(*,*)r(k), "Here 6"
                k=k-1
                !j=j+1
                
                
                r(k)=radius(j)
                rho(k)=interp(model_ref(i,1),model_ref(i+1,1),model_ref(i,2),model_ref(i+1,2),r(k))
                qkappa(k)=interp(model_ref(i,1),model_ref(i+1,1),model_ref(i,5),model_ref(i+1,5),r(k))
                qshear(k)=interp(model_ref(i,1),model_ref(i+1,1),model_ref(i,6),model_ref(i+1,6),r(k))
                eta(k)=interp(model_ref(i,1),model_ref(i+1,1),model_ref(i,9),model_ref(i+1,9),r(k))
                vsv(k)=interp(model_ref(i,1),model_ref(i+1,1),model_ref(i,4),model_ref(i+1,4),r(k))*(1+voro2(j,2))
                vph(k)=interp(model_ref(i,1),model_ref(i+1,1),model_ref(i,7),model_ref(i+1,7),r(k))*(1+voro2(j,4))
                vsh(k)=vsv(k)*sqrt(vshvsv(j))
                vpv(k)=vph(k)*sqrt(phi)
                xi(k)=vshvsv(j)
                vp_data(k)=voro2(j,4)
                vs_data(k)=voro2(j,2)
                
                if (r(k)<model_ref(nic_ref,1)) nic=nic+2
                if (r(k)<model_ref(noc_ref,1)) noc=noc+2        
                ! write(*,*)r(k), "Here 7"
                k=k-1
                j=j+1
            
            elseif (model_ref(i,1)==radius(j)) then
                ! current Node in Voro at same depth as discont in reference model
                r(k)=radius(j)
                rho(k)=model_ref(i+1,2)
                qkappa(k)=model_ref(i+1,5)
                qshear(k)=model_ref(i+1,6)
                eta(k)=model_ref(i+1,9)
                
                vsv(k)=model_ref(i+1,4)*(1+voro2(j-1,2))
                vph(k)=model_ref(i+1,7)*(1+voro2(j-1,4))
                vsh(k)=vsv(k)*sqrt(vshvsv(j-1))
                vpv(k)=vph(k)*sqrt(phi)
                xi(k)=vshvsv(j-1)
                vp_data(k)=voro2(j-1,4)
                vs_data(k)=voro2(j-1,2)
                ! write(*,*)r(k), "Here 8"
                k=k-1
                !j=j+1
                
                
                r(k)=radius(j)
                rho(k)=model_ref(i-1,2)
                qkappa(k)=model_ref(i-1,5)
                qshear(k)=model_ref(i-1,6)
                eta(k)=model_ref(i-1,9)
                vsv(k)=model_ref(i-1,4)*(1+voro2(j,2))
                vph(k)=model_ref(i-1,7)*(1+voro2(j,4))
                vsh(k)=vsv(k)*sqrt(vshvsv(j))
                vpv(k)=vph(k)*sqrt(phi)
                xi(k)=vshvsv(j)
                vp_data(k)=voro2(j,4)
                vs_data(k)=voro2(j,2)
                
                if (r(k)<model_ref(nic_ref,1)) nic=nic+2
                if (r(k)<model_ref(noc_ref,1)) noc=noc+2        
                ! write(*,*)r(k), "Here 9"
                k=k-1
                j=j+1
                i=i-2
                nptfinal=nptfinal-2
                c=c+1


              
            endif
        endif
        
        
    enddo

    ! Account for decreasing nptfinal
    if (c.gt.0) then            
        r(1:nptfinal)=r(2*c+1:nptfinal+2*c)
        rho(1:nptfinal)=rho(2*c+1:nptfinal+2*c)
        vsv(1:nptfinal)=vsv(2*c+1:nptfinal+2*c)
        vpv(1:nptfinal)=vpv(2*c+1:nptfinal+2*c)
        vsh(1:nptfinal)=vsh(2*c+1:nptfinal+2*c)
        vph(1:nptfinal)=vph(2*c+1:nptfinal+2*c)
        qkappa(1:nptfinal)=qkappa(2*c+1:nptfinal+2*c)
        qshear(1:nptfinal)=qshear(2*c+1:nptfinal+2*c)
        eta(1:nptfinal)=eta(2*c+1:nptfinal+2*c)
        vp_data(1:nptfinal)=vp_data(2*c+1:nptfinal+2*c)
        vs_data(1:nptfinal)=vs_data(2*c+1:nptfinal+2*c)
        xi(1:nptfinal)=xi(2*c+1:nptfinal+2*c)

        r(nptfinal+2*c+1:)=0.0
        rho(nptfinal+2*c+1:)=0.0
        vsv(nptfinal+2*c+1:)=0.0
        vpv(nptfinal+2*c+1:)=0.0
        vsh(nptfinal+2*c+1:)=0.0
        vph(nptfinal+2*c+1:)=0.0
        qkappa(nptfinal+2*c+1:)=0.0
        qshear(nptfinal+2*c+1:)=0.0
        eta(nptfinal+2*c+1:)=0.0
        vp_data(nptfinal+2*c+1:)=0.0
        vs_data(nptfinal+2*c+1:)=0.0
        xi(nptfinal+2*c+1:)=0.0
    end if

end subroutine combine_linear_vp

