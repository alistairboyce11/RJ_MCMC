subroutine sum_Ps_Pp_diff(r, vpv, vph, vsv, vsh, nptfinal, nrays, rayp, tdiff, PsPp_error)

    implicit none
    include '../params_p.h'
    include '../params.h'
    real, dimension(nrays_max),intent(in) :: rayp
    integer,intent(in) :: nrays, nptfinal
    integer :: i,j, ndiscont
    real,dimension(mk),intent(in) :: r,vpv,vph,vsv,vsh
    real,dimension(mk):: depth, thick

    real :: dts, dtp, dtheta_p, dtheta_s, dcorr
    real, dimension (mk) :: rad2, thick2, pvel2, svel2
    real, dimension(nrays_max),intent(out) :: tdiff
    logical, intent(out) :: PsPp_error
    rad2(:) = 0
    thick2(:) = 0
    pvel2(:) = 0
    svel2(:) = 0
    PsPp_error = .false.

    ! Start with top and bottom of model already counted.
    ndiscont = 2
    ! Loop through model its inverted.
    do j = nptfinal,2,-1
        thick(j) = r(j) - r(j-1)
        depth(j) = rearth - r(j)
        depth(j-1) = rearth - r(j-1)

        if (depth(j).le.player*1000) then
            if (( vpv(j) .ne. vpv(j-1) ) .or. ( vsv(j) .ne. vsv(j-1) )) then
                if (thick(j) > min_thick) then
                    ! print*, r(j)/1000, depth(j)/1000, thick(j)/1000, vpv(j)/1000, vpv(j-1)/1000, vsv(j)/1000, vsv(j-1)/1000
                    ! print*, "Velocities at top and bottom of layer are not equal"
                    ! print*, "Averaging not acceptable:   (v1+v2)/2"
                    ! print*, "exiting..."
                    PsPp_error = .true.
                    return
                    ! call exit()
                else
                    ! Discontinuity found - this is okay.
                    ndiscont = ndiscont + 1
                end if
            end if
        end if


        if ((depth(j).lt.player*1000).and.(depth(j-1).lt.player*1000)) then
            ! Get correct phase:
            ! print*, r(j)/1000, depth(j)/1000, depth(j-1)/1000, thick(j)/1000, vpv(j)/1000, vsv(j)/1000
            ! Start with vertical wavespeeds vpv, vsv.
            rad2(j)=r(j)/1000
            thick2(j)=thick(j)/1000
            pvel2(j)=vpv(j)/1000
            svel2(j)=vsv(j)/1000

        else if ((depth(j).le.player*1000).and.(depth(j-1).eq.player*1000)) then
            ! print*, "Found the end exactly"
            ! print*, r(j)/1000, depth(j)/1000, depth(j-1)/1000, thick(j)/1000, vpv(j)/1000, vsv(j)/1000
            rad2(j)=r(j)/1000
            thick2(j)=thick(j)/1000
            pvel2(j)=vpv(j)/1000
            svel2(j)=vsv(j)/1000
            ! return
        
        else if ((depth(j).le.player*1000).and.(depth(j-1).gt.player*1000)) then
            ! print*, "Found end..."
            ! print*, "Next model depth: ",depth(j-1)/1000
            thick(j) = r(j) - (rearth - (player*1000.))
            ! print*, "Modified to:..."
            ! print*, r(j)/1000, depth(j)/1000, player,          thick(j)/1000, vpv(j)/1000, vsv(j)/1000

            rad2(j)=r(j)/1000
            thick2(j)=thick(j)/1000
            pvel2(j)=vpv(j)/1000
            svel2(j)=vsv(j)/1000
            ! return
        end if
    end do

    tdiff(:) = 0
    do j = 1,nrays
        do i = nptfinal,2,-1
            ! Check its all being read correctly...
            if (thick2(i) < min_thick) then
                continue
            else


!                 print*, (rayp(j)**2), 1/(svel2(i)**2), 1/(pvel2(i)**2), sqrt(1/(svel2(i)**2)), sqrt(1/(pvel2(i)**2))
                ! Flat Earth has no dependence on radius here. 
!                 tdiff(j) = tdiff(j) + thick2(i)*(sqrt((1/(svel2(i)**2)) - (rayp(j)**2)) - sqrt((1/(pvel2(i)**2)) - (rayp(j)**2)))
                ! Spherical Earth take travel time equations from TT script. Add correction for dtheta. 
                
                dts = (thick2(i)/rad2(i))*((rad2(i)/svel2(i))**2)*(1/sqrt((rad2(i)/svel2(i))**2 - rayp(j)**2))
                dtp = (thick2(i)/rad2(i))*((rad2(i)/pvel2(i))**2)*(1/sqrt((rad2(i)/pvel2(i))**2 - rayp(j)**2))
                
                dtheta_s = ((rayp(j)*thick2(i))/rad2(i))*(1/sqrt((rad2(i)/svel2(i))**2-rayp(j)**2))
                dtheta_p = ((rayp(j)*thick2(i))/rad2(i))*(1/sqrt((rad2(i)/pvel2(i))**2-rayp(j)**2))
                
                dcorr=(dtheta_p-dtheta_s)*rayp(j)
                
                tdiff(j) = tdiff(j) + (dts - dtp + dcorr)
                
                ! if (j==1) then
                !     print*, i, thick2(i), rad2(i), pvel2(i), svel2(i), rayp(j)*pi/180.0, dts, dtp, dtheta_s, dtheta_p
                ! endif

                if (isnan(tdiff(j))) then
                    ! Probably P-wavespeed issue and rayp
                    ! print*,"sum_Ps_Pp_diff FAIL found NaN"
                    ! print*, i, thick2(i), rad2(i), pvel2(i), svel2(i), rayp(j)*pi/180.0, dts, dtp, dtheta_s, dtheta_p
                    PsPp_error = .true.
                    return
                end if
                ! print*, tdiff(j)
            end if
        end do

    end do


end subroutine sum_Ps_Pp_diff