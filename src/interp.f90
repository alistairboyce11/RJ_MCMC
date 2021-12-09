real FUNCTION interp(ray1,ray2,v1,v2,ray)

    !     ..Arguments..
    real,intent(in) :: ray1,ray2,v1,v2,ray
    if (ray1==ray2) then
        interp=v1
    else
        interp=v1+(ray-ray1)*(v2-v1)/(ray2-ray1)
    end if
    RETURN
end function interp