subroutine combine_linear_vp(model_ref,nptref,nic_ref,noc_ref,voro,npt,d_max,&
    r,rho,vpv,vph,vsv,vsh,qkappa,qshear,eta,nptfinal,nic,noc,xi,vp_data)
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
    real,dimension(malay+1) :: depth
    integer :: i,j,k,c
    
    real,dimension(mk),intent(out) :: r,rho,vpv,vph,vsv,vsh,qkappa,qshear,eta
    real,dimension(mk),intent(out) :: xi,vp_data
    real,dimension(mk) :: vs_data
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
        ! depth(i)=rearth-float(int(voro2(i,1)*1000)) ! Round to nearest meter
        depth(i)=rearth-float(int(voro2(i,1))*1000) ! Round to nearest km
    end do
    depth(npt+1)=rearth-d_max*1000
    
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
        ! write(*,*)'depth',depth(j)
        ! write(*,*)j,npt,i,depth(j)
        
        if (j>npt) then
        
            if ((j>npt+1).and.(model_ref(i,1)<depth(npt+1))) then
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
!                 vp_data(k)=(vpv(k)/vsv(k)/vpvs-1)
                vp_data(k)=0.0
                vs_data(k)=0.0
                ! write(*,*)r(k), "Here 1"

                k=k-1
                i=i-1
            
            elseif ((j==npt+1).and.(model_ref(i,1)<=depth(j))) then
                r(k)=depth(j)
                rho(k)=interp(model_ref(i,1),model_ref(i+1,1),model_ref(i,2),model_ref(i+1,2),r(k))
                qkappa(k)=interp(model_ref(i,1),model_ref(i+1,1),model_ref(i,5),model_ref(i+1,5),r(k))
                qshear(k)=interp(model_ref(i,1),model_ref(i+1,1),model_ref(i,6),model_ref(i+1,6),r(k))
                eta(k)=interp(model_ref(i,1),model_ref(i+1,1),model_ref(i,9),model_ref(i+1,9),r(k))
                
                vsv(k)=interp(model_ref(i,1),model_ref(i+1,1),model_ref(i,4),model_ref(i+1,4),r(k))*(1+voro2(j-1,2))
                vph(k)=interp(model_ref(i,1),model_ref(i+1,1),model_ref(i,7),model_ref(i+1,7),r(k))*(1+voro2(j-1,4))
                vsh(k)=vsv(k)*sqrt(vshvsv(j-1))
                vpv(k)=vph(k)*sqrt(phi)

                ! vpv(k)=vsv(k)*vpvsv(j-1)
                ! vsh(k)=vsv(k)*sqrt(vshvsv(j-1))
                ! vph(k)=vpv(k)
                xi(k)=vshvsv(j-1)
                vp_data(k)=voro2(j-1,4)
                ! write(*,*)r(k), "Here 2"
                k=k-1
                !
                
                r(k)=depth(j)
                rho(k)=interp(model_ref(i,1),model_ref(i+1,1),model_ref(i,2),model_ref(i+1,2),r(k))
                qkappa(k)=interp(model_ref(i,1),model_ref(i+1,1),model_ref(i,5),model_ref(i+1,5),r(k))
                qshear(k)=interp(model_ref(i,1),model_ref(i+1,1),model_ref(i,6),model_ref(i+1,6),r(k))
                eta(k)=interp(model_ref(i,1),model_ref(i+1,1),model_ref(i,9),model_ref(i+1,9),r(k))
                
                vsh(k)=interp(model_ref(i,1),model_ref(i+1,1),model_ref(i,8),model_ref(i+1,8),r(k))
                vpv(k)=interp(model_ref(i,1),model_ref(i+1,1),model_ref(i,3),model_ref(i+1,3),r(k))
                vsv(k)=interp(model_ref(i,1),model_ref(i+1,1),model_ref(i,4),model_ref(i+1,4),r(k))
                vph(k)=interp(model_ref(i,1),model_ref(i+1,1),model_ref(i,7),model_ref(i+1,7),r(k))
                
                xi(k)=(vsh(k)/vsv(k))**2
                ! vp_data(k)=(vpv(k)/vsv(k)/vpvs-1)
                vp_data(k)=0.0
                vs_data(k)=0.0
                ! write(*,*)r(k), "Here 3"
                if (r(k)<model_ref(nic_ref,1)) nic=nic+2
                if (r(k)<model_ref(noc_ref,1)) noc=noc+2  
                
                k=k-1
                
                j=j+1
                
                
            elseif ((j==npt+1).and.(model_ref(i,1)>depth(j))) then
                r(k)=model_ref(i,1)
                rho(k)=model_ref(i,2)
                qkappa(k)=model_ref(i,5)
                qshear(k)=model_ref(i,6)
                eta(k)=model_ref(i,9)
                
                vsv(k)=model_ref(i,4)*(1+voro2(j-1,2))
                vph(k)=model_ref(i,7)*(1+voro2(j-1,4))
                vsh(k)=vsv(k)*sqrt(vshvsv(j-1))
                vpv(k)=vph(k)*sqrt(phi)

                ! vpv(k)=vsv(k)*vpvsv(j-1)
                ! vsh(k)=vsv(k)*sqrt(vshvsv(j-1))
                ! vph(k)=vpv(k)
                xi(k)=vshvsv(j-1)
                ! vpvsv_data(k)=voro2(j-1,4)
                vp_data(k)=voro2(j-1,4)
                vs_data(k)=voro2(j-1,2)
                ! write(*,*)r(k), "Here 4"
                k=k-1
                i=i-1 
            endif
        else 
            if (model_ref(i,1)>depth(j)) then
                r(k)=model_ref(i,1)
                rho(k)=model_ref(i,2)
                qkappa(k)=model_ref(i,5)
                qshear(k)=model_ref(i,6)
                eta(k)=model_ref(i,9)
                
                vsv(k)=model_ref(i,4)*(1+voro2(j-1,2))
                vph(k)=model_ref(i,7)*(1+voro2(j-1,4))
                vsh(k)=vsv(k)*sqrt(vshvsv(j-1))
                vpv(k)=vph(k)*sqrt(phi)

                ! vpv(k)=vsv(k)*vpvsv(j-1)
                ! vsh(k)=vsv(k)*sqrt(vshvsv(j-1))
                ! vph(k)=vpv(k)
                xi(k)=vshvsv(j-1)
                ! vpvsv_data(k)=voro2(j-1,4)
                vp_data(k)=voro2(j-1,4)
                vs_data(k)=voro2(j-1,2)
                ! write(*,*)r(k), "Here 5"
                k=k-1
                i=i-1
            elseif (model_ref(i,1)<=depth(j)) then
                
                r(k)=depth(j)
                rho(k)=interp(model_ref(i,1),model_ref(i+1,1),model_ref(i,2),model_ref(i+1,2),r(k))
                qkappa(k)=interp(model_ref(i,1),model_ref(i+1,1),model_ref(i,5),model_ref(i+1,5),r(k))
                qshear(k)=interp(model_ref(i,1),model_ref(i+1,1),model_ref(i,6),model_ref(i+1,6),r(k))
                eta(k)=interp(model_ref(i,1),model_ref(i+1,1),model_ref(i,9),model_ref(i+1,9),r(k))
                
                vsv(k)=interp(model_ref(i,1),model_ref(i+1,1),model_ref(i,4),model_ref(i+1,4),r(k))*(1+voro2(j-1,2))
                vph(k)=interp(model_ref(i,1),model_ref(i+1,1),model_ref(i,7),model_ref(i+1,7),r(k))*(1+voro2(j-1,4))
                vsh(k)=vsv(k)*sqrt(vshvsv(j-1))
                vpv(k)=vph(k)*sqrt(phi)

                ! vpv(k)=vsv(k)*vpvsv(j-1)
                ! vsh(k)=vsv(k)*sqrt(vshvsv(j-1))
                ! vph(k)=vpv(k)
                xi(k)=vshvsv(j-1)
                ! vpvsv_data(k)=voro2(j-1,4)
                vp_data(k)=voro2(j-1,4)
                vs_data(k)=voro2(j-1,2)
                ! write(*,*)r(k), "Here 6"
                k=k-1
                !j=j+1
                
                
                r(k)=depth(j)
                rho(k)=interp(model_ref(i,1),model_ref(i+1,1),model_ref(i,2),model_ref(i+1,2),r(k))
                qkappa(k)=interp(model_ref(i,1),model_ref(i+1,1),model_ref(i,5),model_ref(i+1,5),r(k))
                qshear(k)=interp(model_ref(i,1),model_ref(i+1,1),model_ref(i,6),model_ref(i+1,6),r(k))
                eta(k)=interp(model_ref(i,1),model_ref(i+1,1),model_ref(i,9),model_ref(i+1,9),r(k))
                
                vsv(k)=interp(model_ref(i,1),model_ref(i+1,1),model_ref(i,4),model_ref(i+1,4),r(k))*(1+voro2(j,2))
                vph(k)=interp(model_ref(i,1),model_ref(i+1,1),model_ref(i,7),model_ref(i+1,7),r(k))*(1+voro2(j,4))
                vsh(k)=vsv(k)*sqrt(vshvsv(j))
                vpv(k)=vph(k)*sqrt(phi)
                
                ! vpv(k)=vsv(k)*vpvsv(j)
                ! vsh(k)=vsv(k)*sqrt(vshvsv(j))
                ! vph(k)=vpv(k)
                xi(k)=vshvsv(j)
                ! vpvsv_data(k)=voro2(j,4)
                vp_data(k)=voro2(j,4)
                vs_data(k)=voro2(j,2)
                
                if (r(k)<model_ref(nic_ref,1)) nic=nic+2
                if (r(k)<model_ref(noc_ref,1)) noc=noc+2        
                ! write(*,*)r(k), "Here 7"
                k=k-1
                j=j+1
                
                
            endif
        endif
        
        
    enddo
    nptfinal=nptref+2*npt
    c=1
    i=1
    do while (c.le.nptfinal)
        ! write(*,*) nptfinal, c, r(c),rho(c),vsv(c),vpv(c),vsh(c),&
        ! vph(c),qkappa(c),qshear(c),eta(c),vp_data(c), vs_data(c), xi(c)

        if (r(c).eq.r(c+5)) then
            ! Six points at same depth
            ! write(*,*)"#1",r(c+4:nptfinal)

            r(c+1:nptfinal-4)=r(c+5:nptfinal)
            rho(c+1:nptfinal-4)=rho(c+5:nptfinal)
            vsv(c+1:nptfinal-4)=vsv(c+5:nptfinal)
            vpv(c+1:nptfinal-4)=vpv(c+5:nptfinal)
            vsh(c+1:nptfinal-4)=vsh(c+5:nptfinal)
            vph(c+1:nptfinal-4)=vph(c+5:nptfinal)
            qkappa(c+1:nptfinal-4)=qkappa(c+5:nptfinal)
            qshear(c+1:nptfinal-4)=qshear(c+5:nptfinal)
            eta(c+1:nptfinal-4)=eta(c+5:nptfinal)
            vp_data(c+1:nptfinal-4)=vp_data(c+5:nptfinal)
            vs_data(c+1:nptfinal-4)=vs_data(c+5:nptfinal)
            xi(c+1:nptfinal-4)=xi(c+5:nptfinal)

            r(nptfinal-3:)=0.d0
            rho(nptfinal-3:)=0.d0
            vsv(nptfinal-3:)=0.d0
            vpv(nptfinal-3:)=0.d0
            vsh(nptfinal-3:)=0.d0
            vph(nptfinal-3:)=0.d0
            qkappa(nptfinal-3:)=0.d0
            qshear(nptfinal-3:)=0.d0
            eta(nptfinal-3:)=0.d0
            vp_data(nptfinal-3:)=0.d0
            vs_data(nptfinal-3:)=0.d0
            xi(nptfinal-3:)=0.d0
            
            ! write(*,*)"Six points"
            ! c=c+1
            nptfinal=nptfinal-4
        elseif (r(c).eq.r(c+4)) then
            ! Five points at same depth
            ! write(*,*)"#1",r(c+4:nptfinal)

            r(c+1:nptfinal-3)=r(c+4:nptfinal)
            rho(c+1:nptfinal-3)=rho(c+4:nptfinal)
            vsv(c+1:nptfinal-3)=vsv(c+4:nptfinal)
            vpv(c+1:nptfinal-3)=vpv(c+4:nptfinal)
            vsh(c+1:nptfinal-3)=vsh(c+4:nptfinal)
            vph(c+1:nptfinal-3)=vph(c+4:nptfinal)
            qkappa(c+1:nptfinal-3)=qkappa(c+4:nptfinal)
            qshear(c+1:nptfinal-3)=qshear(c+4:nptfinal)
            eta(c+1:nptfinal-3)=eta(c+4:nptfinal)
            vp_data(c+1:nptfinal-3)=vp_data(c+4:nptfinal)
            vs_data(c+1:nptfinal-3)=vs_data(c+4:nptfinal)
            xi(c+1:nptfinal-3)=xi(c+4:nptfinal)

            r(nptfinal-2:)=0.d0
            rho(nptfinal-2:)=0.d0
            vsv(nptfinal-2:)=0.d0
            vpv(nptfinal-2:)=0.d0
            vsh(nptfinal-2:)=0.d0
            vph(nptfinal-2:)=0.d0
            qkappa(nptfinal-2:)=0.d0
            qshear(nptfinal-2:)=0.d0
            eta(nptfinal-2:)=0.d0
            vp_data(nptfinal-2:)=0.d0
            vs_data(nptfinal-2:)=0.d0
            xi(nptfinal-2:)=0.d0

            ! write(*,*)"Five points"
            ! c=c+1
            nptfinal=nptfinal-3

        elseif (r(c).eq.r(c+3)) then
            ! Four points at same depth
            ! write(*,*)"#2",r(c+3:nptfinal)

            r(c+1:nptfinal-2)=r(c+3:nptfinal)
            rho(c+1:nptfinal-2)=rho(c+3:nptfinal)
            vsv(c+1:nptfinal-2)=vsv(c+3:nptfinal)
            vpv(c+1:nptfinal-2)=vpv(c+3:nptfinal)
            vsh(c+1:nptfinal-2)=vsh(c+3:nptfinal)
            vph(c+1:nptfinal-2)=vph(c+3:nptfinal)
            qkappa(c+1:nptfinal-2)=qkappa(c+3:nptfinal)
            qshear(c+1:nptfinal-2)=qshear(c+3:nptfinal)
            eta(c+1:nptfinal-2)=eta(c+3:nptfinal)
            vp_data(c+1:nptfinal-2)=vp_data(c+3:nptfinal)
            vs_data(c+1:nptfinal-2)=vs_data(c+3:nptfinal)
            xi(c+1:nptfinal-2)=xi(c+3:nptfinal)

            r(nptfinal-1:)=0.d0
            rho(nptfinal-1:)=0.d0
            vsv(nptfinal-1:)=0.d0
            vpv(nptfinal-1:)=0.d0
            vsh(nptfinal-1:)=0.d0
            vph(nptfinal-1:)=0.d0
            qkappa(nptfinal-1:)=0.d0
            qshear(nptfinal-1:)=0.d0
            eta(nptfinal-1:)=0.d0
            vp_data(nptfinal-1:)=0.d0
            vs_data(nptfinal-1:)=0.d0
            xi(nptfinal-1:)=0.d0

            ! write(*,*)"Four points"
            ! c=c+1
            nptfinal=nptfinal-2
        elseif (r(c).eq.r(c+2)) then
            ! Three points at same depth
            ! write(*,*)"#3",r(c+2:nptfinal)
            r(c+1:nptfinal-1)=r(c+2:nptfinal)
            rho(c+1:nptfinal-1)=rho(c+2:nptfinal)
            vsv(c+1:nptfinal-1)=vsv(c+2:nptfinal)
            vpv(c+1:nptfinal-1)=vpv(c+2:nptfinal)
            vsh(c+1:nptfinal-1)=vsh(c+2:nptfinal)
            vph(c+1:nptfinal-1)=vph(c+2:nptfinal)
            qkappa(c+1:nptfinal-1)=qkappa(c+2:nptfinal)
            qshear(c+1:nptfinal-1)=qshear(c+2:nptfinal)
            eta(c+1:nptfinal-1)=eta(c+2:nptfinal)
            vp_data(c+1:nptfinal-1)=vp_data(c+2:nptfinal)
            vs_data(c+1:nptfinal-1)=vs_data(c+2:nptfinal)
            xi(c+1:nptfinal-1)=xi(c+2:nptfinal)

            r(nptfinal:)=0.d0
            rho(nptfinal:)=0.d0
            vsv(nptfinal:)=0.d0
            vpv(nptfinal:)=0.d0
            vsh(nptfinal:)=0.d0
            vph(nptfinal:)=0.d0
            qkappa(nptfinal:)=0.d0
            qshear(nptfinal:)=0.d0
            eta(nptfinal:)=0.d0
            vp_data(nptfinal:)=0.d0
            vs_data(nptfinal:)=0.d0
            xi(nptfinal:)=0.d0




            ! write(*,*)"Three points"
            ! c=c+1
            nptfinal=nptfinal-1
        endif

        ! elseif (r(c).eq.r(c+1)) then
        !     ! Two points at same depth - Normal discontinuity.
        !     ! write(*,*)"Discont."
        !     c=c+1
        ! else
        !     ! Normal Layer
        !     ! write(*,*)"Layer"
        !     c=c+1
        ! endif
        c=c+1
        ! i=i+1
        ! write(*,*)c,"<=",nptfinal
    end do

end subroutine combine_linear_vp

