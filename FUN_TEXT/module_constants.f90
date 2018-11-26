! module_constants.f90

!****************************************************************************
!
! PURPOSE:  Module for the definition of constants.
!
!----------------------------------------------------------------------------
! 
! Fortran-file created by Armin Hofmeister
! 
!----------------------------------------------------------------------------
! History:
! 
! 29.01.2015: create the Fortran-file
! 12.02.2015: add constant for station name length, source name length, session name length, grib epoch name
!             add constants for pi, deg2rad and rad2deg
!             add constants for ray-tracing: "accuracy_elev", "limit_iterations_elev"
! 16.02.2015: add axis of reference ellipsoid WGS84
! 18.02.2015: add "epsilon_refr_ind": epsilon for checking if refractive index is almost equal to 1
! 26.02.2015: add comments
!             add "accuracy_layer", "limit_iterations_layer"
! 02.03.2015: correct variable "accuracy_layer", correct comments
! 13.05.2015: add constant "c" for speed of light
! 02.06.2015: add constant "e2_WGS84" and "e2s_WGS84" for first and second eccentricity e2 and e2' and
!             "cel_WGS84" for the polar curvature radius of reference ellipsoid WGS84
! 11.06.2015: comments
! 10.09.2015: add constants "len_epologname", "len_epolog_filename", "len_session_index_filename"
! 17.09.2015: add constants "len_year", "len_month", "len_day", "len_hour"
! 17.12.2015: add temperature scale offset for transformation from [°C] to [K]
!
!****************************************************************************

module module_constants
    
    ! Variable definitions
    implicit none
    
    save ! Ensures that the content of the memory spaces between different inclusions to different program units remains unchanged, i.e. the value of the
         ! variable in the next call will be the same as when the subroutine was left at the last call and the value does not get lost after the subroutine.
         ! Attention: The value of a variable can be changed during the subroutine!
    
    
    ! define constant for length of station names
    integer, parameter :: len_statname= 8
    
    ! define constant for length of source names
    integer, parameter :: len_souname= 8
    
    ! define constant for length of session names
    integer, parameter :: len_sessname= 14
    
    ! define constants for year string lengths
    integer, parameter :: len_year= 4
    
    ! define constants for month string lengths
    integer, parameter :: len_month= 2
    
    ! define constants for day string lengths
    integer, parameter :: len_day= 2
    
    ! define constants for hour string lengths
    integer, parameter :: len_hour= 2
    
    ! define constant for length of grib epoch name
    integer, parameter :: len_gribepochname= len_year + len_month + len_day + len_hour
    
    ! define constant for length of epolog names
    ! note: length is the sum of the length of the word "epolog1_" or "eplog2_" and the length of the session name and one "_" and the grib epoch name
    integer, parameter :: len_epologname= 8 + len_sessname + 1 + len_gribepochname
    
    ! define constant for length of epolog filenames
    ! note: length is the sum of the length of the parameter "len_epologname" and the file extension ".txt"
    integer, parameter :: len_epolog_filename= len_epologname + 4
    
    ! define constant for length of session index filenames
    ! note: length is the sum of the length of the word "session_index1_" or session_index2_", the length of the session name and the file extension ".txt"
    integer, parameter :: len_session_index_filename= 15 + len_sessname + 4
    
    !----------------------------------------------------------------------------------------------
    
    ! define pi
    ! note: use pi= 4 * atan(1.0d0) as done in Matlab
    !       setting 1.0d0 instead of just 1 is necessary as atan() requires input of real kind and result is of the same kind, here double precision
    double precision, parameter :: pi = 4 * atan(1.0d0)
    
    ! define transformation constant from [°] to [rad]
    double precision, parameter :: deg2rad = pi / 180
    
    ! define transformation constant from [rad] to [°]
    double precision, parameter :: rad2deg = 180 / pi
    
    !----------------------------------------------------------------------------------------------
    
    ! define constant for normal gravity
    double precision, parameter :: gn= 9.80665d0 ! in [m/s^2]
    
    !----------------------------------------------------------------------------------------------
    
    ! define parameters of reference ellipsoid WGS84
    
    ! equatorial axis of WGS84 reference ellipsoid in [m]
    double precision, parameter :: a_WGS84= 6378137.0d0
    
    ! polar axis of WGS84 reference ellipsoid in [m]
    double precision, parameter :: b_WGS84= 6356752.3142d0
    
    ! first eccentricity e2 of WGS84 reference ellipsoid
    ! see also scriptum Terrestrische Bezugsrahmen (2011), equation (1.2a) on page 3
    double precision, parameter :: e2_WGS84= 1 - b_WGS84**2 / (a_WGS84**2)
    
    ! second eccentricity e2s of WGS84 reference ellipsoid
    ! see also scriptum Terrestrische Bezugsrahmen (2011), equation (1.2b) on page 3
    double precision, parameter :: e2s_WGS84= a_WGS84**2 / (b_WGS84**2) - 1
    
    ! polar curvature radius of WGS84 reference ellipsoid
    ! see scriptum Terrestrische Bezugsrahmen (2011), equation (1.10a) on page 5
    double precision, parameter :: cel_WGS84= a_WGS84**2 / b_WGS84
    
    !----------------------------------------------------------------------------------------------
    
    ! define constant for iteration accuracy of outgoing elevation angle
    ! note: setting accuracy for iteration of elevation to 1e-7 rad
    !       leads to a difference in hydrostatic slant delay at 3° elevation smaller than 0.1 mm
    double precision, parameter :: accuracy_elev= 1d-7 ! in [rad]
    
    ! define constant for breaking limit via number of loops in the while loop for iteration of outgoing elevation angle
    integer, parameter :: limit_iterations_elev= 10
    
    ! define constant for iteration accuracy (mean coordinate accuracy) of next intersection point (e.g. in case of Thayer approach)
    double precision, parameter :: accuracy_layer= 1d-7 ! in [°]
    
    ! define constant for breaking limit via number of loops in the while loop for iteration of the next intersection
    ! point position
    integer, parameter :: limit_iterations_layer= 10
    
    !----------------------------------------------------------------------------------------------
    
    ! define constant for epsilon for checking if refractive index is almost equal to 1
    double precision, parameter :: epsilon_refr_ind= 1d-10
    
    !----------------------------------------------------------------------------------------------
    
    ! define constant for speed of light
    double precision, parameter :: c= 299792458d0 ! in [m/s]
    
    !----------------------------------------------------------------------------------------------
    
    ! define temperature scale offset for transformation from [°C] to [K]
    ! note: [K] = [°C] + 273.15
    !       [°C] = [K] - 273.15
    double precision, parameter :: C2K= 273.15d0
    
end module module_constants