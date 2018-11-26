! module_get_RayTrace2D_ref_pwl_global.f90

!****************************************************************************
!
! PURPOSE:  Module containing subroutine to do refined piece-wise linear ray-tracing.
!    
!           Note: Created as module to avoid explicit inferface for the subroutine as needed because of assumed-shape arrays in subroutine.
!           
!           "get_RayTrace2D_ref_pwl_global" is a subroutine to determine the path of the ray trough the atmosphere for
!           a specific initial elevation angle e0. The refined piece-wise linear propagation approach is used
!           (see Hobiger et al. 2008, Fast and accurate ray-tracing algorithms for real-time space geodetic
!           applications using numerical weather models). In this approach the azimuth of the ray is kept
!           constant. Therefore this ray-tracing method is a 2D approach.
!           
!           The outputs of this function are calculated total, hydrostatic and wet zenith and total,
!           hydrostatic and wet slant delays.
!           Furthermore the estimated elevation at the station, the ray-traced outgoing elevation angle as well
!           as his difference to the outgoing elevation from the azel-file (calculated, theoretical) angle
!           are reported.
!           All variables are stored in the structure "delay".
!           
!           --- Version Info ---
!           This version of the function is capable of using a global input grid, which contains the
!           refractive indices. Additionally input of values for refractive indices at station position
!           is necessary to do ray-tracing.
!           --------------------
!
!
! INPUT:
!         az............ azimuth of the observation in [rad]
!         e_outgoing.... outgoing (vacuum) elevation angle in [rad]
!         stat_lat...... latitude of station in [°]
!         stat_lon...... longitude of station in [°]
!         stat_height... (ellipsoidal) height of station in [m]
!         ind_lat1lon1_stat... index of point lat1lon1 in grid used for bilinear interpolation of values at station position
!         ind_lat1lon2_stat... index of point lat1lon2 in grid used for bilinear interpolation of values at station position
!         ind_lat2lon2_stat... index of point lat2lon2 in grid used for bilinear interpolation of values at station position
!         ind_lat2lon1_stat... index of point lat2lon1 in grid used for bilinear interpolation of values at station position
!         nr_h_lev...... total number of available height levels (from vertical interpolation/extrapolation)
!         h_lev_all..... (ellipsoidal) height levels in which the intersection points with the ray path should be estimated; in [m]
!                        attention: station height level is missing in this vector and station height lies somewhere in between
!                                   the input levels
!         start_lev..... index of first height level in "h_lev_all" above station height (needed for
!                        ray-tracing start above station)
!         grid_lat...... grid containing the latitude nodes around the station in [°]
!         grid_lon...... grid containing the longitude nodes around the station in [°]
!         grid_size..... size of the grid
!         dint_lat...... grid intervall for latitude in [°]
!         dint_lon...... grid intervall for longitude in [°]
!         n_stat........ total refractive index at station height and position
!         n_h_stat...... hydrostatic refractive index at station height and position
!         n_w_stat...... wet refractive index at station height and position
!         n............. gridded values of the total refractive index
!         n_h........... gridded values of the hydrostatic refractive index
!         n_w........... gridded values of the wet refractive index
!         start_and_global_check... logical value telling if input grid has the desired starting
!                                   values in latitude and longitude and if it is a global grid in
!                                   latitude and longitude.
!                                   This information is needed when bilinear interpolation is done.
!                                   .TRUE. ... all checks are true, .FALSE. ... at least one check is false
! 
! 
! OUTPUT:
!         delay......... structure containing the following variables:
!           % dz_total........ zenith total delay in [m]
!           % dz_h............ zenith hydrostatic delay in [m]
!           % dz_w............ zenith wet delay in [m]
!           % ds_total_geom... slant total delay including geometric bending effect in [m]
!           % ds_total........ slant total delay in [m]
!           % ds_h_geom....... slant hydrostatic delay including geometric bending effect in [m]
!           % ds_h............ slant hydrostatic delay in [m]
!           % ds_w............ slant wet delay in [m]
!           % e_stat.......... elevation angle at the station in [rad]
!           % e_outgoing_rt... (iteratively) ray-traced outgoing elevation angle in [rad]
!           % diff_e.......... difference: outgoing (theoretical) - outgoing (ray-traced) elevation angle in [rad]
!           % dgeo............ geometric bending effect in [m]
!           % mf_total_geom... value for total mapping factor (includes treatment of geometric bending effect)
!           % mf_h_geom....... value for hydrostatic mapping factor (includes treatment of geometric bending effect)
!           % mf_w............ value for wet mapping factor
!           % break_elev...... logical signaling if a break in the while loop for calculating the
!                              outgoing elevation angle has occured (.FLASE.= no break, .TRUE.= break)
!           not specified in this subroutine:
!               % break_layer..... logical signaling if a break in the while loop for calculating the
!                                  next intersection point has occured (at least for one intersection point) (.FLASE.= no break, .TRUE.= break)
!
!----------------------------------------------------------------------------
! 
! Fortran-file created by Armin Hofmeister
!   based on Matlab scripts created by Armin Hofmeister
! 
!----------------------------------------------------------------------------
! History:
! 
! 18.02.2015: create the Fortran-file based on the Matlab subroutine "get_RayTrace2D_ref_pwl_global.m"
! 19.02.2015: programming
! 25.02.2015: programming
! 26.02.2015: programming
! 02.03.2015: comments
! 07.05.2015: correct comments
! 13.05.2015: correct comments
! 02.06.2015: correct comments
! 09.06.2015: correct setting "break_elev" in case accuracy has been reached at last permitted
!             iteration loop
! 07.09.2015: enhance program code by avoiding unnecessary duplicate calculations of zenith refractive index and zenith delay
! 09.09.2015: correct comments
! 29.09.2015: replace the subroutine for the gaussian curvature radius by the euler radius of curvature
! 10.11.2015: move calculation of second ray point into loop and move persistent definitions at
!             station point out of the loop
! 19.11.2015: avoid unnecessary assignments of e_stat and e_outgoing_rt during iteration loop
!             correct comments
! 17.12.2015: add the total mapping factor calculation

! Changes by Daniel Landskron:
! 05.02.2018: epolog % raytrace is not globally stored anymore because it requires huge disk space and is not
!             needed outside of the get_RayTrace2D_* subroutines
!
!****************************************************************************

module module_get_RayTrace2D_ref_pwl_global

contains

    subroutine get_RayTrace2D_ref_pwl_global( az, &
                                              e_outgoing, &
                                              stat_lat, &
                                              stat_lon, &
                                              stat_height, &
                                              ind_lat1lon1_stat, &
                                              ind_lat1lon2_stat, &
                                              ind_lat2lon2_stat, &
                                              ind_lat2lon1_stat, &
                                              nr_h_lev_all, &
                                              h_lev_all, &
                                              start_lev, &
                                              grid_lat, &
                                              grid_lon, &
                                              grid_size, &
                                              dint_lat, &
                                              dint_lon, &
                                              n_stat, &
                                              n_h_stat, &
                                              n_w_stat, &
                                              n, &
                                              n_h, &
                                              n_w, &
                                              start_and_global_check, &
                                              delay )
        
        ! Define modules to be used
        use module_type_definitions, only: delay_type, raytrace_type
        
        use module_constants, only: deg2rad, rad2deg, accuracy_elev, limit_iterations_elev
        
        use module_R_earth_euler
        
        use module_get_bilint_value
        
        use module_get_ref_pwl_delay
        
        use module_get_ref_pwl_delay_zenith

        !----------------------------------------------------------------------------
    
        ! Variable definitions
        implicit none
    
    
        ! dummy arguments
        !----------------
    
        ! INPUT
    
        ! variable for storing azimuth of the observation in [rad]
        double precision, intent(in) :: az
        
        ! variable for storing outgoing (vacuum) elevation angle in [rad]
        double precision, intent(in) :: e_outgoing
        
        ! variable for storing latitude of station in [°]
        double precision, intent(in) :: stat_lat
    
        ! variable for storing longitude of station in [°]
        double precision, intent(in) :: stat_lon
    
        ! variable for storing (ellipsoidal) height of station in [m]
        double precision, intent(in) :: stat_height
        
        ! variable for storing index of point lat1lon1 in grid used for bilinear interpolation of values at station position
        integer, dimension(2), intent(in) :: ind_lat1lon1_stat
        
        ! variable for storing index of point lat1lon2 in grid used for bilinear interpolation of values at station position
        integer, dimension(2), intent(in) :: ind_lat1lon2_stat
        
        ! variable for storing index of point lat2lon2 in grid used for bilinear interpolation of values at station position
        integer, dimension(2), intent(in) :: ind_lat2lon2_stat
        
        ! variable for storing index of point lat2lon1 in grid used for bilinear interpolation of values at station position
        integer, dimension(2), intent(in) :: ind_lat2lon1_stat
        
        ! variable for storing total number of available height levels (from vertical interpolation/extrapolation)
        integer, intent(in) :: nr_h_lev_all
        
        ! variable for storing (ellipsoidal) height levels in which the intersection points with the ray path should be estimated; in [m]
        ! Attention: Station height level is missing in this vector and station height lies somewhere in between the input levels.
        double precision, dimension(nr_h_lev_all), intent(in) :: h_lev_all
        
        ! variable for storing index of first height level in "h_lev_all" above station height (needed for ray-tracing start above station)
        integer, intent(in) :: start_lev
        
        ! variable for storing grid containing the latitude nodes around the station in [°]
        double precision, dimension(:, :), intent(in) :: grid_lat
        
        ! variable for storing grid containing the longitude nodes around the station in [°]
        double precision, dimension(:, :), intent(in) :: grid_lon
        
        ! variable for storing size of the grid
        integer, dimension(:), intent(in) :: grid_size
        
        ! variable for storing grid interval in latitude in [°]
        double precision, intent(in) :: dint_lat
        
        ! variable for storing grid interval in longitude in [°]
        double precision, intent(in) :: dint_lon
        
        ! variable for storing mean value of the total refractive index between station height and
        ! first interpolated level above at the station's horizontal position
        double precision, intent(in) :: n_stat
                    
        ! variable for storing mean value of the hydrostatic refractive index between station height and
        ! first interpolated level above at the station's horizontal position
        double precision, intent(in) :: n_h_stat
                    
        ! variable for storing mean value of the wet refractive index between station height and
        ! first interpolated level above at the station's horizontal position
        double precision, intent(in) :: n_w_stat
                    
        ! variable for storing gridded values of the total refractive index
        ! (mean between the two original height levels)
        double precision, dimension(:, :, :), intent(in) :: n
                    
        ! variable for storing gridded values of the hydrostatic refractive index
        ! (mean between the two original height levels)
        double precision, dimension(:, :, :), intent(in) :: n_h
                    
        ! variable for storing gridded values of the wet refractive index
        ! (mean between the two original height levels)
        double precision, dimension(:, :, :), intent(in) :: n_w
                    
        ! variable for storing logical value telling if input grid has the desired starting
        ! values in latitude and longitude and if it is a global grid in latitude and longitude.
        ! Note: This information is needed when bilinear interpolation is done.
        !       .TRUE. ... all checks are true, .FALSE. ... at least one check is false
        logical, intent(in) :: start_and_global_check
        
        
        ! OUTPUT
        
        ! type for storing results from ray-tracing concerning the delays, the elevation angles, mapping factors and possible break in ray-tracing iteration
        type(delay_type), intent(out) :: delay
        
    
        ! local variables
        !----------------
        
        ! type for storing results from ray-tracing concerning the ray path
        type(raytrace_type) :: raytrace
        
        ! variable for storing the geocentric heights of height levels (h_lev with added radius of the Earth at station position)
        double precision, dimension(:), allocatable :: R
        
        ! variable for storing the geocentric heights of height levels at half heigth between to consecutive levels (including geocentric station height as first level)
        double precision, dimension(:), allocatable :: R_m
        
        ! variable for storing height differences between two consecutive levels
        double precision, dimension(:), allocatable :: dh
        
        ! variable for storing the linear distances (ray-part) between two intersection points in consecutive height levels
        double precision, dimension(:), allocatable :: s
        
        ! variable for storing the the geocentric y-coordinates of the ray points (= intersection points)
        double precision, dimension(:), allocatable :: y
        
        ! variable for storing the the geocentric z-coordinates of the ray points (= intersection points)
        double precision, dimension(:), allocatable :: z
        
        ! variable for the interpolated (at intersection point) total refractive index at mid height level
        double precision, dimension(:), allocatable :: n_m
        
        ! variable for the interpolated (at intersection point) hydrostatic refractive index at mid height level
        double precision, dimension(:), allocatable :: n_m_h
        
        ! variable for the interpolated (at intersection point) wet refractive index at mid height level
        double precision, dimension(:), allocatable :: n_m_w
        
        ! variable for the interpolated total refractivity along the zenith direction at mid height level
        double precision, dimension(:), allocatable :: n_m_z
        
        ! variable for the interpolated hydrostatic refractivity along the zenith direction at mid height level
        double precision, dimension(:), allocatable :: n_m_h_z
        
        ! variable for the interpolated wet refractivity along the zenith direction at mid height level
        double precision, dimension(:), allocatable :: n_m_w_z
        
        ! variable for the slant total delay between two consecutive levels
        double precision, dimension(:), allocatable :: ds_total_l
        
        ! variable for the slant hydrostatic delay between two consecutive levels
        double precision, dimension(:), allocatable :: ds_h_l
        
        ! variable for the slant wet delay between two consecutive levels
        double precision, dimension(:), allocatable :: ds_w_l
        
        ! variable for the zenith total delay between two consecutive levels
        double precision, dimension(:), allocatable :: dz_total_l
        
        ! variable for the zenith hydrostatic delay between two consecutive levels
        double precision, dimension(:), allocatable :: dz_h_l
        
        ! variable for the zenith wet delay between two consecutive levels
        double precision, dimension(:), allocatable :: dz_w_l
        
        ! variable for storing the a priori bending effect
        double precision :: ap_bend
        
        ! variable for storing the first theta (location dependend elevation angle) = "theoretical" outgoing elevation angle + a priori bending effect
        double precision :: theta_start
        
        ! variable for storing the index of the ray-trace level starting with 1 at the station level
        integer :: ind_trace
        
        ! variable for counting the number of iteration loops when calculating the outgoing elevation angle
        integer :: loop_elev
        
        ! define variable for storing the index of the lower level for refractive index calculation at mid-height level
        integer :: level_lower
        
        ! define variable for storing the index of the upper level for refractive index calculation at mid-height level
        integer :: level_upper
        
        ! define variables for total, hydrostatic and wet slant refractive indices values at the horizontal position of the intersection point
        ! in the original height level below the current mid-height level
        double precision :: n_lower, n_h_lower, n_w_lower
        
        ! define variables for total, hydrostatic and wet zenith refractive indices values at the horizontal position of the station
        ! in the original height level below the current mid-height level
        double precision :: n_z_lower, n_h_z_lower, n_w_z_lower
        
        ! define variables for total, hydrostatic and wet slant refractive indices values at the horizontal position of the intersection point
        ! in the original height level above the current mid-height level
        double precision :: n_upper, n_h_upper, n_w_upper
        
        ! define variables for total, hydrostatic and wet zenith refractive indices values at the horizontal position of the station
        ! in the original height level above the current mid-height level
        double precision :: n_z_upper, n_h_z_upper, n_w_z_upper
        
        ! define variables for mean slant total refractive index (RI) for the lower/upper level determined by vertical
        ! exponential decay
        double precision :: n_m_lower, n_m_upper
        
        ! loop variable
        integer :: i
        
        
        ! CONSTANTS
        
        ! see "module_constants" for "deg2rad", "rad2deg", "accuracy_elev", "limit_iterations_elev"
        
        !----------------------------------------------------------------------------
        
        !============================================================================
        ! Body of the subroutine
        !============================================================================
        
        
        ! Define transformation constants
        
        ! transformation constant from [°] to [rad]
        ! note: see "deg2rad" in "module_constants"
        ! transformation constant from [rad] to [°]
        ! note: see "rad2deg" in "module_constants"
        
        
        ! Define iteration conditions and message status
        
        ! define iteration accuracy of outgoing elevation angle in [rad]
        ! note: see "accuracy_elev" in "module_constants"
        
        ! define breaking limit via number of loops in the while loop for iteration of outgoing elevation angle
        ! note: see "limit_iterations_elev" in "module_constants"
        
        
        ! set status of message for breaking while loop when calculating outgoing elevation angle
        delay % break_elev= .FALSE.
        
        
        ! Calculate intersection points of the ray path with the height levels and determine values necessary for calculating the delays
        
        ! calculate radius of the Earth at station position using formula for the Euler radius of curvature
        ! attention: latitude and azimuth are required in [rad]
        call R_earth_euler( stat_lat * deg2rad, &
                            az, & ! note: already in [rad]
                            raytrace % R_e )
        
        ! get all height levels starting from station height up to last interpolated height level
        
        ! determine output of number of height levels starting at station level up to highest supported height level
        raytrace % nr_h_lev = nr_h_lev_all - start_lev + 2 ! note: +2 as station height level and "start_lev" are also needed
        
        ! allocate the output variable for storing the height levels from station level up to highest level
        allocate( raytrace % h_lev(raytrace % nr_h_lev) )
        
        ! determine the needed height level values
        ! station height level has to be added extra! Omit interpolation levels below station!
        raytrace % h_lev= [ stat_height, h_lev_all(start_lev:nr_h_lev_all) ] ! in [m]
        
        ! calulate geocentric heights of height levels
        
        ! allocate variable for storing the geocentric heights
        allocate( R(raytrace % nr_h_lev) )
        
        ! see scriptum Atmospheric Effects in Geodesy 2012, equation (2.60) on page 36
        ! Note: Adding the Earth radius calculated at the station's latitude to all height levels is sufficient,
        !       although it is an approximation as the radius changes as the ray path alters its latitude!
        !       Therefore a strict solution would need an iterative step to do ray-tracing, re-determine
        !       the radius and do ray-tracing and so on until a certain accuracy would be reached.
        R= raytrace % R_e + raytrace % h_lev ! in [m]
        
        
        ! Refined piece-wise linear approach
        ! see Hobiger et al. 2008, Fast and accurate ray-tracing algorithms for real-time space geodetic
        ! applications using numerical weather models, equations (26-31) on page 8-9
        
        ! define new height levels used for determining the intersection points:
        
        ! allocate variable for storing the geocentric heights at half height (intermediate level)
        allocate( R_m(raytrace % nr_h_lev) )
        
        ! add geocentric station height level as first height level
        R_m(1)= R(1)
        
        ! calculate new geocentric height levels at half height
        R_m(2:raytrace % nr_h_lev)= ( R(1:raytrace % nr_h_lev - 1) + R(2:raytrace % nr_h_lev) ) / 2
        
        ! determine height differences between two consecutive levels
        ! note: size = nr_h_lev - 1
        allocate( dh(raytrace % nr_h_lev - 1) )
        dh= ( R_m(2:raytrace % nr_h_lev) - R_m(1:raytrace % nr_h_lev - 1) )
        
        
        ! convert half-height levels back to ellipsoidal heights and assign height levels of intersection points (half height levels) to "raytrace"-structure
        ! note: this also includes the (real) station height level
        raytrace % h_lev= R_m - raytrace % R_e ! in [m]
        
        
        ! allocate variables for the calculation of the intersection points
        ! note: some variables are also part of the output by the "raytrace" structure
        
        ! allocate all necessary output variables in "raytrace" using "raytrace % nr_h_lev" = number of levels from station level to highest height level
        ! note: setting the variables to NaN is not necessary as sizes are chosen correctly!
        allocate( raytrace % theta(raytrace % nr_h_lev), & ! define variable for elevation (location dependend)
                  raytrace % e(raytrace % nr_h_lev), & ! define variable for elevation (fixed reference to station position)
                  raytrace % anggeo(raytrace % nr_h_lev), & ! define variable for the geocentric angle
                  raytrace % ray_lat(raytrace % nr_h_lev), & ! define variable for latitude of ray points (intersection points with the height levels)
                  raytrace % ray_lon(raytrace % nr_h_lev), & ! define variable for longitude of ray points (intersection points with the height levels)
                  raytrace % ind_lat1lon1_trace(raytrace % nr_h_lev, 2), & ! define variables for the indices of the grid points used for interpolation, note: second dimensions allocated with size=2 as for storing the two indices of lat and lon
                  raytrace % ind_lat1lon2_trace(raytrace % nr_h_lev, 2), & ! define variables for the indices of the grid points used for interpolation, note: second dimensions allocated with size=2 as for storing the two indices of lat and lon
                  raytrace % ind_lat2lon2_trace(raytrace % nr_h_lev, 2), & ! define variables for the indices of the grid points used for interpolation, note: second dimensions allocated with size=2 as for storing the two indices of lat and lon
                  raytrace % ind_lat2lon1_trace(raytrace % nr_h_lev, 2), & ! define variables for the indices of the grid points used for interpolation, note: second dimensions allocated with size=2 as for storing the two indices of lat and lon
                  raytrace % n_total(raytrace % nr_h_lev), & ! define variable for the interpolated (at intersection point) total refractive indices
                  raytrace % n_h(raytrace % nr_h_lev), & ! define variable for the interpolated (at intersection point) hydrostatic refractive indices
                  raytrace % n_w(raytrace % nr_h_lev), & ! define variable for the interpolated (at intersection point) wet refractive indices
                  raytrace % n_total_z(raytrace % nr_h_lev), & ! define variable for the interpolated total refractive indices along the zenith direction
                  raytrace % n_h_z(raytrace % nr_h_lev), & ! define variable for the interpolated hydrostatic refractive indices along the zenith direction
                  raytrace % n_w_z(raytrace % nr_h_lev) ) ! define variable for the interpolated wet refractive indices along the zenith direction
        
        
        ! allocate remaining variables for the calculation of the intersection points
        ! note: setting the variables to NaN is not necessary as sizes are chosen correctly!
        
        ! allocate variable for linear distances
        ! note: size= raytrace % nr_h_lev - 1 as the distances connect two levels and therefore their number is -1 of number of height levels
        allocate( s(raytrace % nr_h_lev - 1) )
        
        ! allocate variables for the geocentric coordinates of the ray points (= intersection points)
        allocate( y(raytrace % nr_h_lev), &
                  z(raytrace % nr_h_lev) )
        
        ! allocate variables for the refractive indices at mid-height level (at mean height between two consecutive original levels)
        ! note: size= raytrace % nr_h_lev, as the station level is included as the first mid height level
        allocate( n_m(raytrace % nr_h_lev), & ! note: this variable could be replaced by raytrace % n_total to save one variable declaration, but program code readability would be worse
                  n_m_h(raytrace % nr_h_lev), & ! note: this variable could be replaced by raytrace % n_h to save one variable declaration, but program code readability would be worse
                  n_m_w(raytrace % nr_h_lev), & ! note: this variable could be replaced by raytrace % n_w to save one variable declaration, but program code readability would be worse
                  n_m_z(raytrace % nr_h_lev), & ! note: this variable could be replaced by raytrace % n_total_z to save one variable declaration, but program code readability would be worse
                  n_m_h_z(raytrace % nr_h_lev), & ! note: this variable could be replaced by raytrace % n_h_z to save one variable declaration, but program code readability would be worse
                  n_m_w_z(raytrace % nr_h_lev) ) ! note: this variable could be replaced by raytrace % n_w_z to save one variable declaration, but program code readability would be worse
                  
        ! allocate variables for the delays between two consecutive levels
        ! note: size= raytrace % nr_h_lev - 1 as the delays connect two levels and therefore their number is -1 of number of height levels
        allocate( ds_total_l(raytrace % nr_h_lev - 1), &
                  ds_h_l(raytrace % nr_h_lev - 1), &
                  ds_w_l(raytrace % nr_h_lev - 1), &
                  dz_total_l(raytrace % nr_h_lev - 1), &
                  dz_h_l(raytrace % nr_h_lev - 1), &
                  dz_w_l(raytrace % nr_h_lev - 1) )
        
        
        ! calculate a priori bending effect
        ! see Hobiger et al. 2008, Fast and accurate ray-tracing algorithms for real-time space geodetic
        ! applications using numerical weather models, equation (32) on page 9
        ap_bend= 0.02d0 * exp(-stat_height / 6000) / tan(e_outgoing) * deg2rad ! in [rad], conversion is necessary!
        
        ! set initial theta (= location dependend) elevation angle
        ! first theta = "theoretical" outgoing elevation angle + a priori bending effect
        ! see scriptum Atmospheric Effects in Geodesy 2012, equation (2.61) on page 36
        theta_start= e_outgoing + ap_bend ! in [rad]
        
        ! set latitude of first intersection point (= station latitude)
        raytrace % ray_lat(1)= stat_lat ! in [°]
        
        ! set longitude of first intersection point (= station longitude)
        raytrace % ray_lon(1)= stat_lon ! in [°]
        
        
        ! get value of refractive index for the first height level (=station position and level)
        
        ! get total, hydrostatic and wet refractive indices value at the horizontal position of the first
        ! intersection point (= station position)
        
        ! set index for ray-trace level of current observation starting at the station
        ind_trace= 1
        
        ! get total, hydrostatic and wet refractive indices value at the horizontal position of the first
        ! intersection point (= station position)
        n_m(ind_trace)= n_stat
        n_m_h(ind_trace)= n_h_stat
        n_m_w(ind_trace)= n_w_stat
        
        ! assign indices of the grid points used for bilinear interpolation at station position and level
        raytrace % ind_lat1lon1_trace(ind_trace,:)= ind_lat1lon1_stat
        raytrace % ind_lat1lon2_trace(ind_trace,:)= ind_lat1lon2_stat
        raytrace % ind_lat2lon2_trace(ind_trace,:)= ind_lat2lon2_stat
        raytrace % ind_lat2lon1_trace(ind_trace,:)= ind_lat2lon1_stat
        
        
        ! define geocentric coordinates of the first ray point = first intersection point
        ! see scriptum Atmospheric Effects in Geodesy 2012, equation (2.63) on page 36
        z(1)= R_m(1) ! in [m]
        y(1)= 0 ! in [m]
        
        ! define the first geocentric angle
        ! see scriptum Atmospheric Effects in Geodesy 2012, equation (2.64) on page 36
        raytrace % anggeo(1)= 0 ! in [rad]
        
        
        ! define first value for "diff_e" (difference between "theoretical" outgoing elevation angle and
        ! ray-traced outgoing elevation angle) used for iteration decision
        ! note: make sure the initial value is higher than the value of "accuracy_elev"
        ! note: see "accuracy_elev" in "module_constants"
        delay % diff_e= 100 * accuracy_elev ! in [rad]
        
        ! initialize variable for counting the loop number for calculating the outgoing elevation angle
        loop_elev= 0
        
        ! loop in order to iterate outgoing elevation angle at the station
        ! note: see "accuracy_elev" in "module_constants"
        iteration_loop: do while ( abs(delay % diff_e) > accuracy_elev ) ! loop until absolute value is lower
            
            ! set first theta (= location dependend) elevation angle) to value retrieved for start (from
            ! initial setting or previous iteration loop)
            raytrace % theta(1)= theta_start
            
            ! set first elevation angle (fixed reference to station position)
            ! see scriptum Atmospheric Effects in Geodesy 2012, equation (2.61) on page 36
            raytrace % e(1)= raytrace % theta(1)
            
            
            ! loop to calculate all remaining ray points using the mid-height levels
            do i= 1, raytrace % nr_h_lev - 1 ! -1 as i+1 is needed inside the loop
                
                ! calculate linear distance to the next intersection point
                ! see scriptum Atmospheric Effects in Geodesy 2012, equation (2.67) on page 37
                s(i)= -R_m(i) * sin(raytrace % theta(i)) + sqrt( R_m(i+1)**2 - R_m(i)**2 * cos(raytrace % theta(i))**2 ) ! in [m]
                
                ! calculate geocentric coordinates of the next intersection point
                ! see scriptum Atmospheric Effects in Geodesy 2012, equation (2.68) on page 37
                z(i+1)= z(i) + s(i) * sin(raytrace % e(i)) ! in [m]
                y(i+1)= y(i) + s(i) * cos(raytrace % e(i)) ! in [m]
                
                ! define the geocentric angle to the next intersection point
                ! see scriptum Atmospheric Effects in Geodesy 2012, equation (2.69) on page 37
                raytrace % anggeo(i+1)= atan(y(i+1) / z(i+1)) ! in [rad]
                
                
                ! calculate latitude and longitude of the next intersection point
                
                ! call subroutine to determine the latitude and longitude of the next intersection point from
                ! geocentric angles and azimuth (incl. check of values lying in the possible intervals)
                ! note: as we want to have a fixed azimuth to "az", it is necessary always to use the station's coordinates and geocentric angle as the first point
                call determine_lat_lon( raytrace % anggeo(1), &
                                        raytrace % anggeo(i+1), &
                                        az, &
                                        stat_lat, &
                                        stat_lon, &
                                        raytrace % ray_lat(i+1), &
                                        raytrace % ray_lon(i+1) )
                
                
                ! determine refractive indices for the next intersection point at the original height levels
                ! below and above the present mid-level
                
                ! index for ray-trace level of current observation
                ind_trace= i + 1
                
                !!! get values for the upper level
                
                ! define level of interpolated data used for bilinear interpolation
                ! attention: In order to get values of correct upper level: actually it is start_lev + (i+1) - 2 = start_lev + i - 1,
                !            where (i+1) is used to get to "next" level and -2 is needed to reduce the data level to the upper level as station level
                !            and start_level introduce each an additional level number in the ray-tracing algorithm (see calculation of raytrace % nr_h_lev)
                level_upper= start_lev + i - 1 ! start_lev+i-1 finds data level above current mid height level
        
                ! slant
                ! get total, hydrostatic and wet slant refractive indices value at the horizontal position of the next intersection point
                ! in the original height level above the current mid-height level
                ! note: see subroutine for optional arguments, use keywords for all arguments as specified in the subroutine to get correct call if some optional arguments are missing in between!
                call get_bilint_value( v1= n, &
                                       v2= n_h, &
                                       v3= n_w, &
                                       level= level_upper, &
                                       POI_lat= raytrace % ray_lat(ind_trace), &
                                       POI_lon= raytrace % ray_lon(ind_trace), &
                                       dint_lat= dint_lat, &
                                       dint_lon= dint_lon, &
                                       grid_lat= grid_lat, &
                                       grid_lon= grid_lon, &
                                       grid_size= grid_size, &
                                       start_and_global_check= start_and_global_check, &
                                       v1_bilint= n_upper, &
                                       v2_bilint= n_h_upper, &
                                       v3_bilint= n_w_upper, &
                                       ind_lat1lon1_out= raytrace % ind_lat1lon1_trace(ind_trace,:), &
                                       ind_lat1lon2_out= raytrace % ind_lat1lon2_trace(ind_trace,:), &
                                       ind_lat2lon2_out= raytrace % ind_lat2lon2_trace(ind_trace,:), &
                                       ind_lat2lon1_out= raytrace % ind_lat2lon1_trace(ind_trace,:) )
                
                
                
                !!! get values for the lower level
                
                ! check if the next intersection point is the second point
                ! note: This is necessary as an exception for the determination of the refractive indices is
                !       needed at the second intersection point.
                if (ind_trace == 2) then
                    
                    ! slant
                    ! get total, hydrostatic and wet slant refractive indices value at the horizontal position of the second intersection point
                    ! in the original height level below the current mid level (=first level)
            
                    ! Attention: Here the refractive indices of the station height level would be needed as grid and
                    ! not only the value at the exact horizontal station position, because the second intersection
                    ! point is usually dislocated compared to the station position.
                    ! Due to RAM concerns when working with many different stations, it is not possible to provide
                    ! (global) gridded data at all station height levels.
                    ! For this reason an approximation is used and the refractive indices from the station level at
                    ! the horizontal position of the station is used to describe the value at the lower level at the
                    ! horizontal position of the second intersection point!
                    n_lower= n_stat
                    n_h_lower= n_h_stat
                    n_w_lower= n_w_stat
                    
                ! for all other intersection points the refractive indices can be determined in the strict
                ! way
                else
                    
                    ! define level of interpolated data used for bilinear interpolation
                    ! attention: In order to get values of correct lower level: actually it is start_lev + (i+1) - 3 = start_lev + i - 2,
                    !            where (i+1) is used to get to "next" level and -3 is needed to reduce the data level to the lower level as station level
                    !            and start_level introduce each an additional level number in the ray-tracing algorithm (see calculation of raytrace % nr_h_lev)
                    level_lower= start_lev + i - 2 ! start_lev+i-2 in order to find data level below current mid height level
                
                    ! slant
                    ! get total, hydrostatic and wet slant refractive indices value at the horizontal position of the next intersection point
                    ! in the original height level below the current mid-height level
                    ! note: see subroutine for optional arguments, use keywords for all arguments as specified in the subroutine to get correct call if some optional arguments are missing in between!
                    call get_bilint_value( v1= n, &
                                           v2= n_h, &
                                           v3= n_w, &
                                           level= level_lower, &
                                           POI_lat= raytrace % ray_lat(ind_trace), &
                                           POI_lon= raytrace % ray_lon(ind_trace), &
                                           dint_lat= dint_lat, &
                                           dint_lon= dint_lon, &
                                           ind_lat1lon1_in= raytrace % ind_lat1lon1_trace(ind_trace,:), & ! input of indices to avoid the call of subroutine to determine them as indices of the current intersection point have already
                                           ind_lat1lon2_in= raytrace % ind_lat1lon2_trace(ind_trace,:), & ! been determined before for the upper level refractive indices determination,
                                           ind_lat2lon2_in= raytrace % ind_lat2lon2_trace(ind_trace,:), & ! grid_lat, grid_lon, grid_size and start_and_global_check are therefore not needed
                                           ind_lat2lon1_in= raytrace % ind_lat2lon1_trace(ind_trace,:), &
                                           v1_bilint= n_lower, &
                                           v2_bilint= n_h_lower, &
                                           v3_bilint= n_w_lower )
                    
                end if
                
                
                ! call subroutine to calculate "n_m_lower" and "n_m_upper" and total, hydrostatic and wet slant delays
                ! set index for ray-trace level
                !ind_trace= i + 1 ! recall not necessary --> commented
                
                call get_ref_pwl_delay( n_m, n_m_h, n_m_w, & ! whole vectors
                                        ds_total_l, ds_h_l, ds_w_l, & ! whole vectors
                                        n_lower, n_h_lower, n_w_lower, & ! scalars
                                        n_upper, n_h_upper, n_w_upper, & ! scalars
                                        R, & ! whole vector
                                        ind_trace, & ! scalar
                                        s, & ! whole vector
                                        n_m_lower, n_m_upper ) ! scalars
                
                
                ! calculate elevation angle of ray path at next intersection point (position specific elevation angle due to earth curvature)
                ! see Hobiger et al. 2008, equation (26) on page 8
                raytrace % theta(i+1)= acos( n_m_lower / n_m_upper * cos(raytrace % theta(i) + raytrace % anggeo(i+1) - raytrace % anggeo(i) ) )
                
                ! calculate elevation angle (fixed reference at ray point)
                ! see scriptum Atmospheric Effects in Geodesy 2012, equation (2.66) on page 36
                raytrace % e(i+1)= raytrace % theta(i+1) - raytrace % anggeo(i+1)
                
            end do
            
            ! difference between "theoretical" outgoing elevation angle and ray-traced outgoing elevation angle
            ! note: the ray-traced outgoing elevation angle is the elevation angle value of the uppermost ray-traced level
            delay % diff_e= e_outgoing - raytrace % e(raytrace % nr_h_lev) ! in [rad]
            
            ! raise counter for loop_elev
            loop_elev= loop_elev + 1
            
            ! break while loop if number of iterations has reached defined limit
            ! note: see "limit_iterations_elev" in "module_constants"
            if (loop_elev >= limit_iterations_elev) then
                ! check if accuracy is still not reached with last iteration
                ! note: break of loop is done afterwards anyway, but setting of "break_elev" to .TRUE. may be skipped
                if ( abs(delay % diff_e) > accuracy_elev ) then
                    ! set message variable to .TRUE. to indicate that break has occurred
                    delay % break_elev= .TRUE.
                end if
                
                ! exit the do while loop ("iteration_loop")
                exit iteration_loop
            end if
            
            ! determine new starting elevation angle at the station based on the difference between
            ! "theoretical" and ray-traced outgoing elevation angle
            theta_start= theta_start + delay % diff_e ! in [rad]
            
        end do iteration_loop
        
        ! get starting elevation angle at station
        delay % e_stat= raytrace % theta(1) ! in [rad]
            
        ! define outgoing elevation angle (location independent)
        ! note: this is the elevation angle value of the uppermost ray-traced level
        delay % e_outgoing_rt= raytrace % e(raytrace % nr_h_lev) ! in [rad]
        
        
        ! Determine the zenith refractive indices and zenith delays between the levels
        ! note: This step is done outside the ray-tracing loop since no iteration is needed in the zenith direction.
        
        ! set index for ray-trace level of current observation starting at the station
        ind_trace= 1
        
        ! zenith refractive indices at first level are equal to slant refractive indices at first level
        ! (station level) and therefore equal to values calculated for the station level
        n_m_z(ind_trace)= n_stat
        n_m_h_z(ind_trace)= n_h_stat
        n_m_w_z(ind_trace)= n_w_stat
        
        
        ! determine refractive indices for the second intersection point at the original height levels
        ! below and above the present mid-level
            
        ! index for ray-trace level of current observation
        ind_trace= 2
            
        ! get values for the lower level
            
        ! zenith
        ! get total, hydrostatic and wet zenith refractive indices value at the horizontal position of
        ! the station in the original height level below the current mid level (=first level)
        ! --> assign value of refractive index at station level and position
        n_z_lower= n_stat
        n_h_z_lower= n_h_stat
        n_w_z_lower= n_w_stat
            
            
        ! get values for the upper level
            
        ! define level of interpolated data used for bilinear interpolation
        level_upper= start_lev ! first original level above station level is "start_lev"
        
        ! zenith 
        ! get total, hydrostatic and wet zenith refractive indices value at the horizontal position of
        ! the station in the original height level above the current mid-height level (= second level)
        ! note: see subroutine for optional arguments, use keywords for all arguments as specified in the subroutine to get correct call if some optional arguments are missing in between!
        call get_bilint_value( v1= n, &
                                v2= n_h, &
                                v3= n_w, &
                                level= level_upper, &
                                POI_lat= stat_lat, &
                                POI_lon= stat_lon, &
                                dint_lat= dint_lat, &
                                dint_lon= dint_lon, &
                                ind_lat1lon1_in= ind_lat1lon1_stat, & ! input of indices to avoid the call of subroutine to determine them as indices have already been determined in the subroutine "calc_refr_ind_at_stations"
                                ind_lat1lon2_in= ind_lat1lon2_stat, & ! grid_lat, grid_lon, grid_size and start_and_global_check are therefore not needed
                                ind_lat2lon2_in= ind_lat2lon2_stat, &
                                ind_lat2lon1_in= ind_lat2lon1_stat, &
                                v1_bilint= n_z_upper, &
                                v2_bilint= n_h_z_upper, &
                                v3_bilint= n_w_z_upper )
        
        ! call subroutine to calculate total, hydrostatic and wet zenith refractive indices and delays
        call get_ref_pwl_delay_zenith( n_m_z, n_m_h_z, n_m_w_z, & ! whole vectors
                                       dz_total_l, dz_h_l, dz_w_l, & ! whole vectors
                                       n_z_lower, n_h_z_lower, n_w_z_lower, & ! scalars
                                       n_z_upper, n_h_z_upper, n_w_z_upper, & ! scalars
                                       ind_trace, & ! scalar
                                       dh ) ! whole vector
        
        
        ! loop over all remaining mid-height levels
        loop_zenith: do i= 2, raytrace % nr_h_lev - 1 ! -1 as i+1 is needed inside the loop
                
            ! determine refractive indices for the next intersection point at the original height levels
            ! below and above the present mid-level
                
            ! index for ray-trace level of current observation
            ind_trace= i + 1
                
                
            !!! get values for the lower level
                
            ! define level of interpolated data used for bilinear interpolation
            ! attention: In order to get values of correct lower level: actually it is start_lev + (i+1) - 3 = start_lev + i - 2,
            !            where (i+1) is used to get to "next" level and -3 is needed to reduce the data level to the lower level as station level
            !            and start_level introduce each an additional level number in the ray-tracing algorithm (see calculation of raytrace % nr_h_lev)
            level_lower= start_lev + i - 2 ! start_lev+i-2 in order to find data level below current mid height level
                
            ! zenith
            ! get total, hydrostatic and wet zenith refractive indices value at the horizontal position of the next intersection point
            ! in the original height level below the current mid-height level
            ! note: see subroutine for optional arguments, use keywords for all arguments as specified in the subroutine to get correct call if some optional arguments are missing in between!
            call get_bilint_value( v1= n, &
                                   v2= n_h, &
                                   v3= n_w, &
                                   level= level_lower, &
                                   POI_lat= stat_lat, &
                                   POI_lon= stat_lon, &
                                   dint_lat= dint_lat, &
                                   dint_lon= dint_lon, &
                                   ind_lat1lon1_in= ind_lat1lon1_stat, & ! input of indices to avoid the call of subroutine to determine them as indices have already been determined in the subroutine "calc_refr_ind_at_stations"
                                   ind_lat1lon2_in= ind_lat1lon2_stat, & ! grid_lat, grid_lon, grid_size and start_and_global_check are therefore not needed
                                   ind_lat2lon2_in= ind_lat2lon2_stat, &
                                   ind_lat2lon1_in= ind_lat2lon1_stat, &
                                   v1_bilint= n_z_lower, &
                                   v2_bilint= n_h_z_lower, &
                                   v3_bilint= n_w_z_lower )
                
                
            !!! get values for the upper level
                
            ! define level of interpolated data used for bilinear interpolation
            ! attention: In order to get values of correct upper level: actually it is start_lev + (i+1) - 2 = start_lev + i - 1,
            !            where (i+1) is used to get to "next" level and -2 is needed to reduce the data level to the upper level as station level
            !            and start_level introduce each an additional level number in the ray-tracing algorithm (see calculation of raytrace % nr_h_lev)
            level_upper= start_lev + i - 1 ! start_lev+i-1 finds data level above current mid height level
            
            ! zenith
            ! get total, hydrostatic and wet zenith refractive indices value at the horizontal position of the next intersection point
            ! in the original height level above the current mid-height level
            ! note: see subroutine for optional arguments, use keywords for all arguments as specified in the subroutine to get correct call if some optional arguments are missing in between!
            call get_bilint_value( v1= n, &
                                    v2= n_h, &
                                    v3= n_w, &
                                    level= level_upper, &
                                    POI_lat= stat_lat, &
                                    POI_lon= stat_lon, &
                                    dint_lat= dint_lat, &
                                    dint_lon= dint_lon, &
                                    ind_lat1lon1_in= ind_lat1lon1_stat, & ! input of indices to avoid the call of subroutine to determine them as indices have already been determined in the subroutine "calc_refr_ind_at_stations"
                                    ind_lat1lon2_in= ind_lat1lon2_stat, & ! grid_lat, grid_lon, grid_size and start_and_global_check are therefore not needed
                                    ind_lat2lon2_in= ind_lat2lon2_stat, &
                                    ind_lat2lon1_in= ind_lat2lon1_stat, &
                                    v1_bilint= n_z_upper, &
                                    v2_bilint= n_h_z_upper, &
                                    v3_bilint= n_w_z_upper )
                
            ! call subroutine to calculate total, hydrostatic and wet zenith refractive indices and delays
            call get_ref_pwl_delay_zenith( n_m_z, n_m_h_z, n_m_w_z, & ! whole vectors
                                           dz_total_l, dz_h_l, dz_w_l, & ! whole vectors
                                           n_z_lower, n_h_z_lower, n_w_z_lower, & ! scalars
                                           n_z_upper, n_h_z_upper, n_w_z_upper, & ! scalars
                                           ind_trace, & ! scalar
                                           dh ) ! whole vector
            
        end do loop_zenith
        
        
        ! Calculate zenith and slant delays ( sum of individual delays between two levels )
        
        ! calculate zenith total delay
        ! note: sum(array [, dim]) without specifying the dimension calculates the sum of the whole array
        delay % dz_total= sum(dz_total_l) ! in [m]
        ! calculate zenith hydrostatic delay
        delay % dz_h= sum(dz_h_l) ! in [m]
        ! calculate zenith wet delay
        delay % dz_w= sum(dz_w_l) ! in [m]
        
        ! calculate slant total delay
        delay % ds_total=sum(ds_total_l) ! in [m]
        ! calculate slant hydrostatic delay
        delay % ds_h=sum(ds_h_l) ! in [m]
        ! calculate slant wet delay
        delay % ds_w=sum(ds_w_l) ! in [m]
        
        
        ! Calculate geometric bending effect
        
        ! see scriptum Atmospheric Effects in Geodesy 2012, equation (2.75) on page 38
        ! note: as there are only vectors and a scalar, it is not necessary to specify a dimension for the sum() function
        delay % dgeo= sum( s - cos(raytrace % e(1:raytrace % nr_h_lev-1) - delay % e_outgoing_rt) * s) ! in [m]
        
        
        ! Calculate delays including geometric bending effect
        
        ! slant hydrostatic delay + geometric bending effect
        delay % ds_h_geom = delay % ds_h + delay % dgeo ! in [m]
        
        ! calculate slant total delay + geometric bending effect
        delay % ds_total_geom = delay % ds_total + delay % dgeo ! in [m]
        
        
        ! Calculate mapping factor values
        
        ! calculate value for total mapping factor (includes treatment of geometric bending effect)
        delay % mf_total_geom = delay % ds_total_geom / delay % dz_total
        
        ! calculate value for hydrostatic mapping factor (includes treatment of geometric bending effect)
        ! see Böhm and Schuh (2003) "Vienna Mapping Functions" in Working Meeting on European VLBI for Geodesy and Astrometry
        ! equation (A26)
        delay % mf_h_geom = delay % ds_h_geom / delay % dz_h
        
        ! calculate value for wet mapping factor
        ! see Böhm and Schuh (2003) "Vienna Mapping Functions" in Working Meeting on European VLBI for Geodesy and Astrometry
        ! equation (A27)
        delay % mf_w = delay % ds_w / delay % dz_w
        
        
        ! Assign output variables to output structure variable "raytrace" and transform the values
        ! note: this is done only for variables that are still missing or for variables that need to be transformed e.g. from [rad] to [deg]
        raytrace % az= az * rad2deg ! in [°], "rad2deg" specified in "module_constants"
        raytrace % theta= raytrace % theta * rad2deg ! in [°]
        raytrace % e= raytrace % e * rad2deg ! in [°]
        raytrace % anggeo= raytrace % anggeo * rad2deg ! in [°]
        
        raytrace % n_total= n_m
        raytrace % n_h= n_m_h
        raytrace % n_w= n_m_w
        raytrace % n_total_z= n_m_z
        raytrace % n_h_z= n_m_h_z
        raytrace % n_w_z= n_m_w_z
        
    end subroutine get_RayTrace2D_ref_pwl_global
    
end module module_get_RayTrace2D_ref_pwl_global