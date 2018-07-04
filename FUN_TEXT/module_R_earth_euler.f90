! module_R_earth_euler.f90

!****************************************************************************
!
! PURPOSE:  Module containing elemental subroutine to determine the radius of the Earth at a specific
!           point defined through the latitude and the azimuth (elemental procedure enabling array input).
!
!           Note: Elemental subroutine needs explicit interface. Therefore it is contained inside a module to create this interface automatically.
!           
!           "R_earth_euler" is a subroutine to determine the radius of the Earth at a specific
!           point (or various points) as the Euler radius of curvature defined through the latitude and the
!           azimuth.
! 
!           This subroutine is capable of calculating the radius for different points specified through an input
!           vector of latitudes in [°] and azimuths in [°].
!
!           Attention: Latitude and azimuth are required in [rad]!
!
!
! INPUT:
!         lat... vector or scalar containing latitude(s) in [rad] of point(s) where the Earth radius should be calculated
!         az.... vector or scalar containing azimuth(s) in [rad] of point(s) where the Earth radius should be calculated
!                The azimuth is the angle between the ellpisoidal meridian plane of the point and
!                the vertical plane containing the normal vector containing the point and the source
!                point (ray-tracing target).
!
!
! OUTPUT:
!         R_e... radius of Earth at the specific point defined through Euler radius of curvature in [m]
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
! 28.09.2015: create the Fortran-file based on the Matlab-file "R_euler.m"
! 29.09.2015: programming
!
!****************************************************************************

module module_R_earth_euler
    
contains
    
    elemental subroutine R_earth_euler( lat, &
										az, &
                                        R_e )
    
        ! Define modules to be used
        use module_constants, only: a => a_WGS84, b => b_WGS84, e2s => e2s_WGS84, cel => cel_WGS84
        ! note: axis names of WGS84 are renamed to "a" and "b",
        !       second eccentricity name of WGS84 is renamed to "e2s" and
        !       polar curvature radius of WGS84 is renamed to "cel"
    
        !----------------------------------------------------------------------------
    
        ! Variable definitions
        implicit none
    
    
        ! dummy arguments
        !----------------
    
        ! INPUT
    
        ! input variable for latitude in [rad] of point where the Earth radius should be calculated
        ! attention: value in [rad]!
        double precision, intent(in) :: lat
		
		! input variable for azimuth in [rad] to the source point for ray-tracing
        ! attention: value in [rad]!
        double precision, intent(in) :: az
    
    
        ! OUTPUT
    
        ! output variable for radius of Earth at the specific point defined through Euler radius of curvature in [m]
        double precision, intent(out) :: R_e
        
        
        ! local variables
        !----------------
        
        ! define variable for the auxiliary parameter V
        double precision :: V
        
        ! define variable for the meridional curvature radius M
        double precision :: M
        
        ! define variable for the radius of curvature N in the prime vertical
        double precision :: N
        
        
        ! CONSTANTS
        
        ! see "module_constants" for "e2s_WGS84", "cel_WGS84"
    
        !----------------------------------------------------------------------------
    
        !============================================================================
        ! Body of the subroutine
        !============================================================================

        ! Define constants
        
        ! define second eccentricity of reference ellipsoid WGS84
        ! note: the constant is defined in the "module_constants"
        ! attention: the name is redefined within the above "use-statement" to "e2s"
        
        ! define polar curvature radius of reference ellipsoid WGS84
        ! note: the constant is defined in the "module_constants"
        ! attention: the name is redefined within the above "use-statement" to "cel"


        ! Calculate radius of the Earth using formula for the Euler radius of curvature

        ! calculate auxiliary parameter V
		! see scriptum Terrestrische Bezugsrahmen (2011), equation (1.3) on page 3
		V= sqrt(1 + e2s * cos(lat)**2)
        
        ! calculate meridional curvature radius M in [m]
		! see scriptum Terrestrische Bezugsrahmen (2011), equation (1.10a) on page 5
		M= cel / V**3
        
        ! calculate radius of curvature N in the prime vertical in [m]
		! see scriptum Terrestrische Bezugsrahmen (2011), equation (1.11) on page 5
		N= cel / V
        
        ! calculate Euler radius of curvature in [m]
		! see scriptum Terrestrische Bezugsrahmen (2011), equation (1.27) on page 11
		R_e = M * N / (M * sin(az)**2 + N * cos(az)**2)

    end subroutine R_earth_euler
    
end module module_R_earth_euler