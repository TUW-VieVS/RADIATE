! bilinear_interpolation.f90

!****************************************************************************
!
! PURPOSE:  Subroutine for bilinear interpolation.
!           
!           "bilinear_interpolation" is a subroutine to determine the value of a parameter at a
!           specific point in a two dimensional grid of (ellipsoidal) coordinates.
!           The four nearest grid points around the point of interest are used to calculate the bilinear
!           interpolation.
!           The formulas are taken from Hobiger et al. 2008, Fast and accurate ray-tracing algorithms for
!           real-time space geodetic applications using numerical weather models, equation (18-19) on page 7
!           
!           The output of this function is the value calculated at the specific point of interest (POI).
!
!
! INPUT:
!         lat_POI....... latitude of the point of interest in [°]
!         lon_POI....... longitude of the point of interest in [°]
!         dint_lat...... grid interval for latitude in [°]
!         dint_lon...... grid interval for longitude in [°]
!         N_lat1_lon1... value of parameter at grid point with lat1 and lon1
!         N_lat1_lon2... value of parameter at grid point with lat1 and lon2
!         N_lat2_lon2... value of parameter at grid point with lat2 and lon2
!         N_lat2_lon1... value of parameter at grid point with lat2 and lon1
! 
! 
! OUTPUT:
!         N_POI....... value of parameter at the point of interest (POI)
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
! 04.02.2015: create the Fortran-file based on the Matlab-file "bilinear_interpolation.m"
! 09.02.2015: comments
! 11.05.2015: correct the calculation of xi and eta by using modulo() which is correct as mod() is wrong
!
!****************************************************************************
    
    
subroutine bilinear_interpolation( lat_POI, &
                                   lon_POI, &
                                   dint_lat, &
                                   dint_lon, &
                                   N_lat1_lon1, &
                                   N_lat1_lon2, &
                                   N_lat2_lon2, &
                                   N_lat2_lon1, &
                                   N_POI )
    
    ! Define modules to be used
    
    !----------------------------------------------------------------------------
    
    ! Variable definitions
    implicit none
    
    
    ! dummy arguments
    !----------------
    
    ! INPUT
    
    ! input variable for latitude of the point of interest in [°]
    double precision, intent(in) :: lat_POI
    
    ! input variable for longitude of the point of interest in [°]
    double precision, intent(in) :: lon_POI
    
    ! input variable for the grid resolution in latitude in [°]
    double precision, intent(in) :: dint_lat
    
    ! input variable for the grid resolution in longitude in [°]
    double precision, intent(in) :: dint_lon
    
    ! input variable for storing value of parameter at grid point with lat1 and lon1
    double precision, intent(in) :: N_lat1_lon1
    
    ! input variable for storing value of parameter at grid point with lat1 and lon2
    double precision, intent(in) :: N_lat1_lon2
    
    ! input variable for storing value of parameter at grid point with lat2 and lon2
    double precision, intent(in) :: N_lat2_lon2
    
    ! input variable for storing value of parameter at grid point with lat2 and lon1
    double precision, intent(in) :: N_lat2_lon1
    
    
    ! OUTPUT
    
    ! output variable for 
    double precision, intent(out) :: N_POI
    
    
    ! local variables
    !----------------
    
    ! variable for auxiliary parameter "eta"
    double precision :: eta
    
    ! variable for auxiliary parameter "xi"
    double precision :: xi
    
    
    !----------------------------------------------------------------------------
    
    !============================================================================
    ! Body of the subroutine
    !============================================================================
    
    ! calculate auxiliary parameters eta and xi
    ! see Hobiger et al. 2008, between equation (18) and (19) on page 7
    ! using the grid resolution allows to calculate the differences in latitude or longitude from the
    ! grid point to the POI using the modulus
    ! the difference between two grid points is represented directly through the resolution
    ! note: modulo() in Fortran works like mod() in MATLAB (if non-zero result then the sign of the second argument is retained)
    !       mod() in Fortran is rem() in MATLAB and would deliver result with wrong sign in case of negative latitude values
    eta= ( modulo( lat_POI, dint_lat ) ) / dint_lat
    xi= ( modulo( lon_POI, dint_lon ) ) / dint_lon
    
    ! calculate bilinear interpolation and receive parameter at POI
    ! see Hobiger et al. 2008, equation (18) on page 7
    N_POI= (1 - xi) * (1 - eta) * N_lat1_lon1 + xi * (1 - eta) * N_lat1_lon2 + xi * eta * N_lat2_lon2 + (1 - xi) * eta * N_lat2_lon1

end subroutine bilinear_interpolation