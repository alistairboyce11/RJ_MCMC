subroutine get_rayp(rp_file, nrays, rayp)

    implicit none
    include '../params_p.h'
    include '../params.h'
    real, dimension(nrays_max),intent(out) :: rayp
    integer,intent(out) :: nrays
    integer :: i
    character(len=64),intent(in) :: rp_file
    ! allocate memory      
    ! allocate ( rayp(nrays) )
    
    ! opening the file for reading
    open (4, file = rp_file, status = 'old')
    read(4,*) nrays
    do i = 1,nrays
       read(4,*) rayp(i)
    end do
    close(4)

    ! write(*,*)"IMPORTED RAY-P file, num rays: ",nrays
    ! do i = 1,nrays
    !     write(*,*) rayp(i)
    ! end do


    do i = 1,nrays
        ! Flat earth divide by radius(in km)
        ! rayp(i) = rayp(i)*180.0/(radius*pi)
        ! For spherical, remove radius factor (Eq21, Hu et al., 2015, JAES)
        rayp(i) = rayp(i)*180.0/pi
    end do
    
end subroutine get_rayp