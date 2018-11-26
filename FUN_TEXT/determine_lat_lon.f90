! determine_lat_lon.f90

!****************************************************************************
!
! PURPOSE:  Subroutine to determine the latitude and longitude of a new point.

!           "determine_lat_lon" is a subroutine to determine the latitude and longitude of a new point (point 2)
!           using the difference in geocentric angles (referenced to the same origin, e.g. the station point)
!           and the azimuth.
!           Attention: This means for ray-tracing that the values psi1, lat1 and lon1 must always be the
!           station position values, otherwise the observational azimuth value would have to be changed in order to get correct results
!           as the geometry alters if the starting point (first point) changes.
!   
!           Note: creating this subroutine as elemental subroutine is not possible as the inside called subroutine "adapt_lat_lon_to_interval"
!                 would then also have to be elemental, which is not possible when using a write-statement inside with output to external-file-unit or *!
!
!           This subroutine uses transformation formulas from the scriptum "Höhere Geodäsie" chapter 8.5.
!
!
! INPUT:
!         psi1... geocentric angle to point 1 in [rad]
!         psi2... geocentric angle to point 2 in [rad]
!         az..... azimuth from point 1 to point 2 in [rad] (reference 0° is north direction!)
!         lat1... latitude of point 1 in [°]
!         lon1... longitude of point 1 in [°]
!
!
! OUTPUT:
!         lat2... latitude of point 2 in [°]
!         lon2... longitude of point 2 in [°]
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
! 17.02.2015: create the Fortran-file based on the Matlab-file "determine_lat_lon.m"
! 18.06.2015: add catch of erroneous dpsi as a result of numerical issues due to elevations near to
!             90°
! 31.08.2015: exchange cotan() through cos() / sin() as gfortran does not know cotan()
! 09.09.2015: correct comments
! 14.09.2015: correct erroneous use of initilized dpsi instead of calculated dpsi
! 01.10.2015: add correction of lat2 calculation in case aux for asin(aux) is out of the possible
!             interval due to numeric reasons
! 02.10.2015: comments
! 19.11.2015: correct comments
!
!****************************************************************************
 

subroutine determine_lat_lon( psi1, &
                              psi2, &
                              az, &
                              lat1_in, &
                              lon1_in, &
                              lat2, &
                              lon2 )
    
    ! Define modules to be used
    use module_constants, only: deg2rad, rad2deg
    
    !----------------------------------------------------------------------------
    
    ! Variable definitions
    implicit none
    
    
    ! dummy arguments
    !----------------
    
    ! INPUT
    
    ! input variable for geocentric angle to point 1 in [rad]
    double precision, intent(in) :: psi1
        
    ! input variable for geocentric angle to point 2 in [rad]
    double precision, intent(in) :: psi2
        
    ! input variable for azimuth from point 1 to point 2 in [rad] (reference 0° is north direction!)
    double precision, intent(in) :: az
        
    ! input variable for latitude of point 1 in [°]
    double precision, intent(in) :: lat1_in
        
    ! input variable for longitude of point 1 in [°]
    double precision, intent(in) :: lon1_in
        
        
    ! OUTPUT
    
    ! output variable for latitude of point 2 in [°]
    double precision, intent(out) :: lat2
        
    ! output variable for longitude of point 2 in [°]
    double precision, intent(out) ::  lon2
    
    
    ! local variables
    !----------------
    
    ! local variable for latitude of point 1 in [rad]
    ! note: local version of lat1 is necessary as input "lat1_in" needs to be transformed to [rad] and inten(in)-variables can not be redefined
    double precision :: lat1
        
    ! local variable for longitude of point 1 in [rad]
    ! note: local version of lon1 is necessary as input "lon1_in" needs to be transformed to [rad] and inten(in)-variables can not be redefined
    double precision :: lon1
        
    ! variable for storing the difference in geocentric angle between the two points = arc length
    double precision :: dpsi
    
    ! variable for storing the auxiliary value for the lat2 calculation
    double precision :: aux
        
    ! auxiliary variables for calculation of lon2
    double precision :: X, Y
        
        
    ! CONSTANTS
        
    ! see "module_constants" for "deg2rad", "rad2deg"
    
    
    !----------------------------------------------------------------------------
    
    !============================================================================
    ! Body of the subroutine
    !============================================================================
    
    
    ! calculate arc length between the points 1 and 2
    ! Attention: error is made due to setting the earth radius to a constant value, which means using a
    ! sphere as earth representation!
    
    ! calculate difference in geocentric angle between the two points = arc length
    dpsi= psi2 - psi1 ! in [rad]
    
    
    ! Determine lat2 and lon2

    ! check if abs(dpsi) is very near to 0
    ! Attention: As long as dpsi is positive and not exactly 0 the calculation of lat2 and lon2 would
    !            in principle work fine.
    !            First consideration: If dpsi is very near to 0 it means that lat2 and lon2 should be
    !            almost the same as lat1 and lon1.
    !            Second consideration: As numerical issues occur if elevation angles of
    !            approximately 90° in terms of accuracy of representation of pi/2 are used for
    !            ray-tracing. Then dpsi may become a) slightly negative or b) exactly 0 leading to
    !            a) wrong lon2 by +180° or a complete failure as cot(0)= cos(0)/sin(0) is not possible!
    ! Solution: If dpsi is near to 0 the calculation of lat2 and lon2 is not necessary as the angles can
    !           be set to the values of lat1 and lon1.
    !           In case of dpsi being near 0 and negativ due to numerical issues the calculation of lon2
    !           would deliver a result that is wrong by +180°, so instead of getting lat2= 200° if
    !           lat1= 200°, lat2 would result in 20° after adapting the originally resulting longitude
    !           of 380° to the posiible longitude interval using "adapt_lat_lon_to_interval".
    ! Note: A negative value of dpsi should in principle never occur. If it is negative then only due to
    !       numerical reasons and only very near below 0 and just .
    !       A dpsi being smaller than -1e-15 will not be caught here, but it should not occur as this
    !       would mean that the ray changed azimuth direction by 180° compared to the starting point!
    ! 
    ! Explanation: An epsilon of 1d-15 is set for dpsi as an elevation angle that is differing from
    !              pi/2 by up to 6.8d-9 rad = 3.9d-7 deg, which is less than the iteration accuracy of
    !              RADIATE, results in a first dpsi of larger than 1d-14. So the ray-tracing should not
    !              become affected if the elevation is not meant to be exactly 90°.
    
    if ( abs(dpsi) < 1d-15 ) then
        ! set lat2 to the value of lat1
        lat2= lat1_in ! in [°]
    
        ! set lat2 to the value of lat1
        lon2= lon1_in ! in [°]
    
        ! note: an adaption to the interval supported by latitude and longitude is not necessary as the
        ! original values of lat1 and lon1 should already fit to the interval
    
    else
        
        ! Define transformation constants
        
        ! transformation constant from [rad] to [°]
        ! note: see "rad2deg" in "module_constants"
        
        ! transformation constant from [°] to [rad]
        ! note: see "deg2rad" in "module_constants"
        
        
        ! Transform latitude and longitude of point 1 to [rad]
        lat1= lat1_in * deg2rad ! in [rad]
        lon1= lon1_in * deg2rad ! in [rad]
    
    
        ! Determine latitude and longitude of point 2
        
        ! see scriptum Höhere Geodäsie chapter 8.5, equations (8.48)
        
        
        ! calculate latitude of point 2
        
        ! calculate auxiliary variable
        ! attention: Due to numeric reasons it can be the case that "aux" is just out of the possible interval [-1,1].
        !            The following asin() function determines lat2 and will then be NaN. Therefore a treatment of aux is necessary
        !            in this case!
        aux= sin(lat1) * cos(dpsi) + cos(lat1) * sin(dpsi) * cos(az)
        
        ! check if the auxiliary variable has a value that is out of range of [-1, 1] for the asin() function
        ! note: This can happen due to numeric reasons in case a lat2 of exactly -90° or 90° should be calculated.
        ! note: Using aint(aux) solves only numeric problems in case aux is in the intervals ]-2, -1[ and ]1, 2[.
        !       In case of aux being smaller or equal to -2 or larger or equal to 2 asin(aux) will still produce NaN since aux is then
        !       still out of the possible interval for asin(). But this case should not happen as this would mean that something went
        !       completely wrong in the calculations.
        if ( abs(aux) > 1) then
            
            ! truncate the value of aux, which solves the numeric problem for an aux that is just outside the possible interval for asin()
            ! and calculate lat2
            ! note: aint() truncates a value to a whole number, this means in other words rounding towards zero (like fix in MATLAB).
            !       Without specifying a kind value aint() preserves the input kind, but input must be a kind of real!
            !       aint() is different from anint(), which rounds towards the nearest integer!
            lat2= asin( aint(aux) ) ! in [rad]
			
        ! in case of aux being in the possible interval for asin()
        else
            ! calculate latitude of point 2 without corrections 
            lat2= asin( aux ) ! in [rad]
        end if
        
        
        ! calculate longitude of point 2
        Y= sin(az)
        X= cos(dpsi) / sin(dpsi) * cos(lat1) - sin(lat1) * cos(az) ! note: use cos()/sin() instead of cotan() as gfortran does not support cotan()
        
        ! catch erroneous longitude calculation if sin(az)=sin(180*pi/180) is not exactly zero (due to numeric reasons)
        ! set Y to 0
        if ( abs(Y) < 1d-13 ) then
            Y= 0
        end if
        
        lon2= atan2(Y, X) + lon1 ! in [rad]
        
        
        ! Transform latitude and longitude of point 2 to [°]
        
        lat2= lat2 * rad2deg ! in [rad]
        lon2= lon2 * rad2deg ! in [rad]
        
        
        ! correct values for lat2 and lon2 if they are out of the specified intervals
        ! Note: If lat and lon are still outside the interval after the adaption the program is forced to stop!
        call adapt_lat_lon_to_interval( lat2, lon2 )
        
    end if
        
end subroutine determine_lat_lon