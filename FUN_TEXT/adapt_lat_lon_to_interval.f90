! adapt_lat_lon_to_interval.f90

!****************************************************************************
!
! PURPOSE:  Subroutine to determine the latitude and longitude of a new point.

!           "adapt_lat_lon_to_interval" is a subroutine to adapt input latitude and longitude values, so that
!           they lie in their specific interval for possible values.
!           In case of values that would lie out of range ([90°,-90] for latitude and [0°,360°[ for longitude)
!           corrections are applied.
!
!           Attention:
!           The corrections are only valid for lat and lon that are a maximum of one interval out of the interval!
!
!           The output of this function are the (corrected if necessary) latitude and longitude values.
!           
!           Note: creating this subroutine as elemental subroutine is not possible as a write-statement inside with output to external-file-unit or * is present!
!
! INPUT and OUTPUT:
!         lat................ latitude in [°]; output: in the interval [90°,-90]
!         lon................ longitude in [°]; output: in the interval [0°,360°[
!
!
!----------------------------------------------------------------------------
! 
! Fortran-file created by Armin Hofmeister
!   based on Matlab scripts created by Armin Hofmeister
! 
!----------------------------------------------------------------------------
! History:
! 
! 17.02.2015: create the Fortran-file based on the Matlab-file "adapt_lat_lon_to_interval.m"
! 19.02.2015: comments
! 20.04.2015: correct comments and messages
!
!****************************************************************************
 

subroutine adapt_lat_lon_to_interval( lat, &
                                      lon )
    
    ! Define modules to be used
    
    !----------------------------------------------------------------------------
    
    ! Variable definitions
    implicit none
    
    
    ! dummy arguments
    !----------------
    
    ! INPUT and OUTPUT
        
    ! input/output variable for latitude in [°]
    ! note: output value is input value, but value may be adapted to fit the interval [-90°, 90°]
    double precision, intent(in out) :: lat
        
    ! input/output variable for longitude in [°]
    ! note: output value is input value, but value may be adapted to fit the interval [0°, 360°[
    double precision, intent(in out) :: lon
    
    
    !----------------------------------------------------------------------------
    
    !============================================================================
    ! Body of the subroutine
    !============================================================================


    ! correct longitude value outside the possible interval of [0°, 360°]
    ! note: treatment of 360° = 0° follows as last step after latitude correction, this doesn't influence the longitude change of +-180° in case of crossing the pole)
    ! correction for longitude > 360°
    if (lon > 360) then
        lon= lon - 360
    ! correction for longitude <0°    
    else if (lon < 0) then
        lon= lon + 360
    end if
        
    ! correct latitude value outside the possible interval of [-90°, 90°] --> latitude crosses the pole
    ! plus additional longitude correction if latitude crosses the poles (> 90° or < -90°)
    if (lat > 90) then
        lat= 180 - lat
        ! correct longitude: -180°, if longitude is equal or above 180°
        if (lon >= 180) then
            lon= lon - 180
        ! correct longitude: +180°, if longitude is below 180°
        else if (lon < 180) then
            lon= lon + 180
        end if
    else if (lat < -90) then
        lat= -180 - lat
        ! correct longitude: -180°, if longitude is equal or above 180°
        if (lon >= 180) then
            lon= lon - 180
        ! correct longitude: +180°, if longitude is below 180°
        else if (lon < 180) then
            lon= lon + 180
        end if
    end if
        
    ! set longitude to 0° if it is 360° to get the interval [0°,360°[
    if (lon == 360) then
        lon= 0
    end if
    
    
    ! test if lat and lon are still outside the intervals
    ! Note: This case can only happen if the initial not adapted values were outside the intervals by more than one total span of the interval!
    if ( (lat < -90) .OR. (lat > 90) .OR. (lon < 0) .OR. (lon > 360) ) then
        ! write warning message
        ! note: output of values through ss --> no plus sign, f7.2 --> leading zero in case of shorter value
        write(unit= *, fmt= '(a)') 'Warning in "adapt_lat_lon_to_interval": Adapting latitude and/or longitude value to interval not successful! Adapted latitude and/or longitude are still out of possible range!'
        write(unit= *, fmt= '(a)') 'Continuing of ray-tracing would lead to erroneous or impossible grid point indices in "determine grid points"! Program stopped!'
        ! stop the program
        stop
    end if
        
end subroutine adapt_lat_lon_to_interval