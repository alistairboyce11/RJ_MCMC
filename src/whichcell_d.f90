subroutine whichcell_d(point,voro,npt,ind)
include '../params.h'

    real, intent(in) :: point,voro(malay,4)
    integer, intent(in) :: npt
    integer, intent(out) :: ind

    integer i
    real mmax

    mmax=-10
    ind = 1

    do i =1,npt

        if (voro(i,1)<=point) then

            if (voro(i,1)>mmax) then
                ind=i
                mmax=voro(i,1)
            endif

        endif
    enddo


    ! ind=1
    ! do i=1,npt
    ! !write(*,*)point,abs(point-voro(i,1)),abs(point-voro(ind,1))
    !     if (abs(point-voro(i,1))<abs(point-voro(ind,1))) then
    !         ind=i
    !     endif
    !  !   write(*,*)ind
    ! enddo



    return
end