! module_ell2xyz.f90
    

!****************************************************************************
!
! PURPOSE:  Module containing elemental subroutine to transform ellipsoidal coordinates to geocentric
!           cartesian coordinates.
!           Elemental procedure enables array input, but needs explicit interface, which can be given
!           if subroutine is included in a module.
!            
!           This subroutine transforms the ellipsoidal latitude, longitude and height coordinates of a point to
!           geocentric cartesian coordinates. Therefore WGS84 is used as the reference ellipsoid.
!           The reference ellipsoid's center defines the origin of the cartesian system.
!
!           Note: This subroutine supports the input of vectors of the same length defining geodetic coordinates
!                 for multiple points as it is an elemental subroutine.
! 
!
! INPUT:
!        lat... ellipsoidal latitude in [°]
!        lon... ellipsoidal longitude in [°]
!        h..... ellipsoidal height (height above the reference ellipsiod) in [m]
! 
! OUTPUT:
!        x... geocentric x-coordinate in [m] (x-axis has the longitude 0°)
!        y... geocentric y-coordinate in [m] (y-axis has the longitude 90°)
!        z... geocentric z-coordinate in [m] (z-axis is the rotation axis of the ellipsoid)
!
!----------------------------------------------------------------------------
! 
! Fortran-file created by Armin Hofmeister
!   based on Matlab scripts created by Armin Hofmeister
! 
!----------------------------------------------------------------------------
! History:
! 
! 02.06.2015: create the Fortran-file based on the Matlab-file "ell2xyz.m"
!
!****************************************************************************

module module_ell2xyz

contains
    
    elemental subroutine ell2xyz( lat, &
                                  lon, &
                                  h, &
                                  x, &
                                  y, &
                                  z )

        ! Define modules to be used
        use module_constants, only: a => a_WGS84, b => b_WGS84, e2 => e2_WGS84, deg2rad
        ! note: axis names and first eccentricity of WGS84 are renamed to "a" and "b" and "e2"
    
    
        !----------------------------------------------------------------------------
    
        ! Variable definitions
        implicit none
    
    
        ! dummy arguments
        !----------------
    
        ! INPUT
    
        ! variable for input of ellipsoidal latitude in [°]
        double precision, intent(in) :: lat
    
        ! variable for input of ellipsoidal longitude in [°]
        double precision, intent(in) :: lon
        
        ! variable for input of ellipsoidal height in [m]
        double precision, intent(in) :: h
    
    
        ! OUTPUT
    
        ! variable for output of geocentric x-coordinate in [m]
        double precision, intent(out) :: x
        
        ! variable for output of geocentric x-coordinate in [m]
        double precision, intent(out) :: y
        
        ! variable for output of geocentric x-coordinate in [m]
        double precision, intent(out) :: z
    
    
        ! local variables
        !----------------
    
        ! variable for storing the latitude value in [rad]
        double precision :: lat_rad
        
        ! variable for storing the longitude value in [rad]
        double precision :: lon_rad
        
        ! variable for storing the normal radius of curvature N
        double precision :: N
        
        
        ! CONSTANTS
        
        ! see "module_constants" for "a_WGS84", "b_WGS84", "e2_WGS84", "deg2rad"
    
    
        !----------------------------------------------------------------------------
    
        !============================================================================
        ! Body of the subroutine
        !============================================================================


        ! define coefficients
        
        ! conversion constant from [deg] to [rad]
        ! note: see "module_constants" for constant "deg2rad"
        
        ! define axis of reference ellipsoid WGS84 in [m]
        ! note: the constants are defined in the "module_constants"
        ! attention: their names are redefined within the above "use-statement" to "a" and "b" and "e2"
        
        
        !---------------------------------------------------------------------------------------------------
        
        ! transform ellisoidal (geodetic) latitude and longitude from [°] to [rad]
        ! note: For the later transformation to geocentric cartesian coordinates using sin and cos it is necessary
        !       to convert the ellipsoidal (geodetic) latitude and longitude from [°] to [rad].
        
        
        ! transform ellipsoidal (geodetic) latitude value from [°] to [rad]
        
        ! note: see "module_constants" for constant "deg2rad"
        ! note: extra variable for storing latitude in [rad] is needed as intent(in) variable can not be redefined
        lat_rad= lat * deg2rad
        
        ! transform ellipsoidal (geodetic) longitude value from [°] to [rad]
        ! note: extra variable for storing longitude in [rad] is needed as intent(in) variable can not be redefined
        lon_rad= lon * deg2rad
        
        
        !---------------------------------------------------------------------------------------------------
        
        ! Transformation from ellipsoidal to geocentric cartesian coordinates
        ! This section performs the transformation of the ellipsoidal to the geocentric cartesian
        ! coordinates.
        ! Equations taken from G. C. Jones (2002) "New solutions for the geodetic coordinate transformation"
        ! page 1 eq. (1-3).
        ! Or see scriptum "Terrestrische Bezugsrahmen" page 3 eq. (1.3) for W^2 and page 5 eq. (1.11) for N and page 6 eq. (1.15) for x, y, z.
        
        ! get the eccentricity e2
        ! note: see "module_constants" for constant "e2_WGS84"
                
        ! calculate the normal radius of curvature N
        N= a / sqrt( 1 - e2 * sin(lat_rad)**2 ) ! in [m]
        
        ! calculate the geocentric x-coordinate
        x= ( N + h ) * cos(lat_rad) * cos(lon_rad) ! in [m]
        
        ! calculate the geocentric y-coordinate
        y= ( N + h ) * cos(lat_rad) * sin(lon_rad) ! in [m]
        
        ! calculate the geocentric z-coordinate
        z= (N * (1 - e2) + h) * sin(lat_rad) ! in [m]
        
    end subroutine ell2xyz
end module module_ell2xyz
    