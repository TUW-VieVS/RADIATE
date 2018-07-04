! module_calc_refr_ind_at_stations.f90

!****************************************************************************
!
! PURPOSE:  Module containing subroutine to calculate the refractive index at specific station positions.
!           
!           Note: Created as module to avoid explicit inferface for the subroutine as needed because of assumed-shape arrays in subroutine.
!           
!           "calc_refr_ind_at_stations" is a subroutine to  calculate the refractive index
!           at the specific station positions at their actual heights.
!
!           Prior to this step interpolation of the meterological values is done, which will be used for
!           calculating the refractive index at the newly interpolated height = station levels and positions.
!           
!           Two step approach of calculating n_h and n_w is the fastest method (at least in MATLAB) under minimum memory usage!
! 
! 
! INPUT:
!         POIs_lat................ ellips. latitude of POIs in [°], interval [-90°,90°]
!         POIs_lon................ ellips. longitude of POIs in [°], interval [0°,360°[
!         POIs_h_ell.............. ellips. heigth of POIs in [m]
!         meteo_int.............. structure containing (vertically interpolated, highly resolved)
!                                  meteorological values at specific ellipsoidal height levels used
!                                  for interpolating the values at the station heights
!          contents:
!           % h........................ interpolation/extrapolation height levels in [m] to increase the
!                                       vertical resolution of the meteorological parameters
!           % nr_h_lev................. number of interpolation/extrapolation height levels
!           % upper_limit_ECMWF........ uppermost height level with data from the ECMWF grib-file, in [m]
!           % upper_limit_atm.......... uppermost height level assuming the end of the atmosphere, in [m]
!           % p........................ interpolated/extrapolated pressure values in [hPa]
!           % T........................ interpolated/extrapolated temperature values in [K]
!           % wvpr..................... interpolated/extrapolated water vapour pressure values in [hPa]
!           % n_h...................... hydrostatic refractive index for the interpolated/extrapolated
!                                       height levels
!           % n_w...................... wet refractive index for the interpolated/extrapolated height
!                                       levels
!           % n_total.................. array containing the calculated total refractive index for the
!                                       interpolated/extrapolated height levels
!
!         dint_lat............ grid intervall for latitude in [°]
!         dint_lon............ grid intervall for longitude in [°]
!         grid_lat............ grid containing the latitude nodes around the station in [°]
!         grid_lon............ grid containing the longitude nodes around the station in [°]
!         grid_size........... size of the grid
!         start_and_global_check... logical value telling if input grid has the desired starting
!                                   values in latitude and longitude and if it is a global grid in
!                                   latitude and longitude.
!                                   This information is needed when bilinear interpolation is done.
!                                   .TRUE. ... all checks are true, .FALSE. ... at least one check is false
!         suspend_raytr........... variable which reports, if ray-tracing is not possible, because POI is not
!                                  supported by the grid or has not been found in the catalogue
!                                  . FALSE. ... ray-tracing will be done for this specific station
!                                  .TRUE. ... ray-tracing will be suspended for this specific station
! 
! 
! OUTPUT:
!         meteo_int... structure (per station) containing the interpolated meteorological values at each exact
!                       station position and height
!          contents:
!            .p........... pressure in [hPa]
!            .T........... temperature in [K]
!            .wvpr........ water vapour pressure in [hPa]
!            .n_h......... hydrostatic refractive index
!            .n_w......... wet refractive index
!            .n_total..... total refractive index
!            .start_lev... index of first height level in meteo_int above station height (needed for
!                          ray-tracing start above station)
!            .ind_lat1lon1..... index of point lat1lon1 in grid used for bilinear interpolation of values at station position
!            .ind_lat1lon2..... index of point lat1lon2 in grid used for bilinear interpolation of values at station position
!            .ind_lat2lon2..... index of point lat2lon2 in grid used for bilinear interpolation of values at station position
!            .ind_lat2lon1..... index of point lat2lon1 in grid used for bilinear interpolation of values at station position
!
!----------------------------------------------------------------------------
! 
! Fortran-file created by Armin Hofmeister
!   based on Matlab scripts created by Armin Hofmeister
! 
!----------------------------------------------------------------------------
! History:
! 
! 28.01.2015: create the Fortran-file based on the Matlab-file "calc_refr_ind_at_stations.m"
! 29.01.2015: programming
! 03.02.2015: programming
! 04.02.2015: programming
! 09.02.2015: programming
! 20.04.2015: correct an error message and comment
! 05.05.2015: declare dimension of input variables for POIs depending on "POIs_lat"
! 13.05.2015: comments
! 02.06.2015: correct comments
! 06.10.2015: add check if station height is lower than the first height level from the vertical
!             interpolation
!
!****************************************************************************

module module_calc_refr_ind_at_stations

contains
    
    subroutine calc_refr_ind_at_stations( POIs_lat, &
                                          POIs_lon, &
                                          POIs_h_ell, &
                                          meteo_int, &
                                          dint_lat, &
                                          dint_lon, &
                                          grid_lat, &
                                          grid_lon, &
                                          grid_size, &
                                          start_and_global_check, &
                                          suspend_raytr, &
                                          meteo_stat )

        ! Define modules to be used
        use module_type_definitions, only: meteo_int_type, meteo_stat_type
        use module_constants, only: gn
        use module_meteo_constants
        use module_get_bilint_value
        use, intrinsic :: ieee_arithmetic ! intrinsic module to e.g. specify NaN and Inf
        
        !----------------------------------------------------------------------------
    
        ! Variable definitions
        implicit none
    
    
        ! dummy arguments
        !----------------
    
        ! INPUT
    
        ! vector for input vector of ellips. latitudes of POIs
        double precision, dimension(:), intent(in) :: POIs_lat
    
        ! vector for input vector of ellips. longitudes of POIs
        ! note: size is the number of input "POIs_lat"
        double precision, dimension(size(POIs_lat)), intent(in) :: POIs_lon
    
        ! vector for input vector of ellips. heights of POIs
        ! note: size is the number of input "POIs_lat"
        double precision, dimension(size(POIs_lat)), intent(in) :: POIs_h_ell
    
        ! input of structure type for storing vertically interpolated meteorological data and some additional data
        type(meteo_int_type), intent(in) :: meteo_int
    
        ! input variable for storing the grid resolution in latitude
        double precision, intent(in) :: dint_lat
        
        ! input variable for storing the grid resolution in longitude
        double precision, intent(in) :: dint_lon
    
        ! input variable for storing grid of latitude values
        ! note: dimension = nlat, nlon
        double precision, dimension(:, :), intent(in) :: grid_lat
        
        ! input variable for storing grid of longitude values
        ! note: dimension = nlat, nlon
        double precision, dimension(:, :), intent(in) :: grid_lon
        
        ! input variable for storing the grid size in latitude and longitude
        integer, dimension(:), intent(in) :: grid_size
        
        ! input variable for check if starting values of first grid point (latitude and longitude) apply to requirements
        ! as well as the grid coverage overall (latitude and longitude) is global
        logical, intent(in) ::  start_and_global_check
    
        ! input variable that reports if ray-tracing for a specific observing station needs to be suspended (default: .FALSE. --> no suspension == ray-tracing will be done)
        logical, dimension(:), intent(in) :: suspend_raytr
    
    
        ! OUTPUT
    
        ! output structure for storing interpolated meteorological data at the station position
        ! note: structure needs to be allocated outside in "get_meteo_undulation_refr_glob" with size = number of POIs
        type(meteo_stat_type), dimension(:), intent(out) :: meteo_stat
    
    
        ! local variables
        !----------------
    
        ! CONSTANTS
    
        ! Define constants and variables concerning meteorological fundamentals
        
        ! see "module_meteo_constants"
    
    
        ! Define gravitational acceleration in [m/s^2]
        ! this gravity value is needed for interpolating the pressure
        ! define variable for normal gravity, which is sufficient (no latitude or height corrections
        ! necessary for sub-mm ray-tracing precision)
        ! see "module_constants"
    
    
        ! OTHER LOCAL VARIABLES
    
        ! define variable for storing the "NaN" value
        double precision :: NaN
    
        ! define variable for storing number of POIs
        integer :: nr_POIs
    
        ! define variable for storing the (loop) index of the current station (POI)
        integer :: ind_POI
        
        ! define variable for checking which height levels are below or equal to the current station height level
        logical, dimension(:), allocatable :: check_lower
        
        ! index of height level below or equal to the station height of current POI
        integer :: ind_lower
        
        ! index of first height level just above to the station height of current POI
        integer :: ind_upper
        
        ! define variable for pressure at station position in lower height level
        double precision :: p_lower_lev
        
        ! define variable for water vapour pressure at station position in lower height level
        double precision :: wvpr_lower_lev
        
        ! define variable for temperature at station position in lower height level
        double precision :: T_lower_lev
        
        ! define variable for pressure at station position in upper height level
        double precision :: p_upper_lev
        
        ! define variable for water vapour pressure at station position in upper height level
        double precision :: wvpr_upper_lev
        
        ! define variable for temperature at station position in upper height level
        double precision :: T_upper_lev
        
        ! define variable for lower height level value
        double precision :: h_lower_lev
        
        ! define variable for upper height level value
        double precision :: h_upper_lev
        
        ! define variable for mean ellipsoidal height between upper and lower level
        double precision :: h_mean_uplow
        
        ! define variable for storing virtual temperature at station position in lower height level
        double precision :: Tv_lower_lev
    
        ! define variable for storing virtual temperature at station position in upper height level
        double precision :: Tv_upper_lev
        
        ! define variable for storing the auxiliary parameter C for logarithmic interpolation of the water vapour pressure for the station height level
        double precision :: C
        
        ! define variable for storing the density of dry air for the station height level
        double precision :: rho_d_stat
    
        ! define variable for storing the density of water vapour for the station height level
        double precision :: rho_w_stat
    
        ! define variable for storing the inverse compressibility factor of water vapour for the station height level
        double precision :: Zw_inv_stat
    
        !----------------------------------------------------------------------------
    
        !============================================================================
        ! Body of the subroutine
        !============================================================================
        
        ! determine the NaN-value as signaling NaN (producing exceptions when used in calculation)
        NaN= ieee_value(0.0d0, ieee_signaling_nan) ! 0.0d0 specifies double precision kind
        
        ! Define number of stations (POIs) for which interpolation of the meteorological values will be done
        nr_POIs= size(POIs_lat)
        
        ! allocate variable for checking if height levels are lower or equal to station height
        ! note: the dimension is equal to the number of height levels in meteo_int % h
        allocate(check_lower(meteo_int % nr_h_lev))
    
        ! Interpolate meteorological values at the station position at exact height level
        ! This section uses gridded (global) meteorological values at specific height levels to interpolate
        ! the values for input stations at their position (horizontal and height).
    
        ! loop over all stations
        do ind_POI= 1, nr_POIs
    
            ! check if ray-tracing for this station can be done
            ! check if "suspend_raytr" is .FALSE.
            if (.NOT. suspend_raytr(ind_POI)) then
                
                ! check if the station height is lower than the first height level of meteo_int % h
                ! --> not possible to calculate the meteorological values at the station height
                if ( POIs_h_ell(ind_POI) < meteo_int % h(1) ) then
                    ! report error message
                    write(unit= *, fmt= '(a, f0.3, a, f0.3, a)') 'Error: Height of a station (', POIs_h_ell(ind_POI), ' m) is too small for calculating meteorological data at the station as it is lower than the lowest vertically interpolated height level (', meteo_int % h(1), ' m)! Program stopped!'
                    ! stop the program
                    stop
                end if
                
                
                ! get index of height level from provided values that is just below the current station
                ! attention: meteo_int % h must be sorted ascending otherwise wrong upper index will be received!
                
                ! determine, which height levels are below or equal to the station height
                check_lower= ( meteo_int % h <= POIs_h_ell(ind_POI) )
                
                ! determine index of uppermost height level that is still below or equal to station height
                ! note/ attention: the index can be found by counting the number of height levels that are below or equal to the station height
                !                  as indexing starts with 1 and the input height levels are required to be sorted ascending
                ind_lower= count(check_lower)
                
                ! determine index of height level that is just above the current station
                ind_upper= ind_lower + 1
                
                !assign value "start_lev" to structure of current POI
                meteo_stat(ind_POI) % start_lev= ind_upper
                
                ! check if "ind_upper" is bigger or equal to number of height levels supported by meteo_int % h
                ! --> in this case there are not sufficient meteorological data (height levels) to do ray-tracing
                if ( ind_upper >= meteo_int % nr_h_lev ) then
                    ! report error message
                    write(unit= *, fmt= '(a)') 'Error: Height of a station is too large in order to do ray-tracing with remaining interpolated height levels! Program stopped!'
                    ! stop the program
                    stop
                end if
                
                
                ! get the meteorological values received from the vertically interpolated values for the specific lower and
                ! upper level at the exact horizontal station position
                
                ! determine values at horizontal station position in the lower level
                ! get lower level pressure  in [hPa]
                ! get lower level water vapour pressure in [hPa]
                ! get lower level temperature in [K]
                ! note: see subroutine for optional arguments, use keywords for all arguments as specified in the subroutine to get correct call if some optional arguments are missing in between!
                call get_bilint_value( v1= meteo_int % p, &
                                       v2= meteo_int % wvpr, &
                                       v3= meteo_int % T, &
                                       level= ind_lower, &
                                       POI_lat= POIs_lat(ind_POI), &
                                       POI_lon= POIs_lon(ind_POI), &
                                       dint_lat= dint_lat, &
                                       dint_lon= dint_lon, &
                                       grid_lat= grid_lat, &
                                       grid_lon= grid_lon, &
                                       grid_size= grid_size, &
                                       start_and_global_check= start_and_global_check, &
                                       v1_bilint= p_lower_lev, &
                                       v2_bilint= wvpr_lower_lev, &
                                       v3_bilint= T_lower_lev, &
                                       ind_lat1lon1_out= meteo_stat(ind_POI) % ind_lat1lon1, &
                                       ind_lat1lon2_out= meteo_stat(ind_POI) % ind_lat1lon2, &
                                       ind_lat2lon2_out= meteo_stat(ind_POI) % ind_lat2lon2, &
                                       ind_lat2lon1_out= meteo_stat(ind_POI) % ind_lat2lon1 )
                
                ! determine values at horizontal station position in the upper level
                ! get upper level pressure  in [hPa]
                ! get upper level water vapour pressure in [hPa]
                ! get upper level temperature in [K]
                ! note: see subroutine for optional arguments, use keywords for all arguments as specified in the subroutine to get correct call if some optional arguments are missing in between! 
                call get_bilint_value( v1= meteo_int % p, &
                                       v2= meteo_int % wvpr, &
                                       v3= meteo_int % T, &
                                       level= ind_upper, &
                                       POI_lat= POIs_lat(ind_POI), &
                                       POI_lon= POIs_lon(ind_POI), &
                                       dint_lat= dint_lat, &
                                       dint_lon= dint_lon, &
                                       ind_lat1lon1_in= meteo_stat(ind_POI) % ind_lat1lon1, & ! input of indices to avoid the call of subroutine to determine them again as indices should be the same as in the lower level
                                       ind_lat1lon2_in= meteo_stat(ind_POI) % ind_lat1lon2, & ! grid_lat, grid_lon, grid_size and start_and_global_check are therefore not needed
                                       ind_lat2lon2_in= meteo_stat(ind_POI) % ind_lat2lon2, &
                                       ind_lat2lon1_in= meteo_stat(ind_POI) % ind_lat2lon1, &
                                       v1_bilint= p_upper_lev, &
                                       v2_bilint= wvpr_upper_lev, &
                                       v3_bilint= T_upper_lev )
                
                                       ! note: usage of optional input of indices lat1lon1 etc. as values should be the same as for subroutine call above in lower level
                                       !       repeated output of these indices is therefore not necessary
                
                
                ! get values for lower and upper height levels
                
                ! set lower height level
                h_lower_lev= meteo_int % h(ind_lower)
                ! set upper height level
                h_upper_lev=meteo_int % h(ind_upper)
                
                ! Calculate the mean ellipsoidal height between upper and lower level
                h_mean_uplow= (h_lower_lev + h_upper_lev) / 2
                
                
                ! interpolation of meteorological values (p, T, wvpr) at the actual station height
                
                ! --> if the actual interpolation height level (= station height) is lower or equal than the mean height between the upper and lower level
                ! --> interpolation for the pressure is done using the pressure and virtual temperature values from the lower level 
                if ( POIs_h_ell(ind_POI) <= h_mean_uplow ) then
                    
                    ! calculate the virtual temperature using meteorological values from the lower level
                    ! see scriptum Atmospheric Effects in Geodesy 2012, equation (1.67) on page 15 or
                    ! see Nafisi et al. 2012, Ray-traced tropospheric delays in VLBI analysis , equation (21) on page 7
                    ! + Kleijer 2004 equation (4.13) + (4.3))
                    Tv_lower_lev= T_lower_lev * p_lower_lev / (p_lower_lev - (1 - Mw / Md) * wvpr_lower_lev) ! in [K]
                    
                    ! interpolate the pressure --> linear interpolation of the logarith-transformed pressure values (=exponential interpolation for the pressure domain)
                    ! see Nafisi et al. 2012, Ray-traced tropospheric delays in VLBI analysis , equation (20) on page 7
                    meteo_stat(ind_POI) % p = p_lower_lev * exp( -(POIs_h_ell(ind_POI) - h_lower_lev) * gn / (Rd * Tv_lower_lev) ) ! in [hPa]
                    
                ! --> if the actual interpolation height level is higher than the mean height level between the upper and lower level
                ! --> interpolation for the pressure is done using the pressure and virtual temperature values from the upper level
                else
                    
                    ! calculate the virtual temperature using meteorological values from the upper level
                    ! see scriptum Atmospheric Effects in Geodesy 2012, equation (1.67) on page 15 or
                    ! see Nafisi et al. 2012, Ray-traced tropospheric delays in VLBI analysis , equation (21) on page 7
                    ! + Kleijer 2004 equation (4.13) + (4.3))
                    Tv_upper_lev= T_upper_lev * p_upper_lev / (p_upper_lev - (1 - Mw / Md) * wvpr_upper_lev) ! in [K]
                    
                    ! interpolate the pressure --> linear interpolation of the logarith-transformed pressure values (=exponential interpolation for the pressure domain)
                    ! see Nafisi et al. 2012, Ray-traced tropospheric delays in VLBI analysis , equation (20) on page 7
                    meteo_stat(ind_POI) % p= p_upper_lev * exp( -(POIs_h_ell(ind_POI) - h_upper_lev) * gn / (Rd * Tv_upper_lev) ) ! in [hPa]
                    
                end if
                
                ! interpolate the temperature --> linear interpolation
                meteo_stat(ind_POI) % T= T_upper_lev + (T_lower_lev - T_upper_lev) / (h_lower_lev - h_upper_lev) * (POIs_h_ell(ind_POI) - h_upper_lev)
                
                ! interpolate the water vapour pressure
                ! --> logarithmic interpolation, if possible otherwise use linear interpolation
                ! for logarithmic interpolation eqautions see
                ! Nafisi et al. 2012, Ray-traced tropospheric delays in VLBI analysis , equation (22) on page 7
                ! or Böhm and Schuh 2003, Vienna mapping functions
                
                ! check if the auxiliary parameter C for logarithmic interpolation can be calculated correctly
                ! in case that wvpr_lower_lev or wvpr_upper_lev are 0 or if wvpr_lower_lev is equal to wvpr_upper_lev
                ! C would deliver 0 or Inf --> no correct logarithmic interpolation possible
                ! --> switch to linear interpolation in these cases
                
                ! check if linear interpolation is necessary
                ! note: in case of wvpr_lower_lev == wvpr_upper_lev the result of the interpolation is wvpr_lower_lev, but an extra if case would not be more efficient
                if ( (wvpr_lower_lev == 0) .OR. (wvpr_upper_lev == 0) .OR. (wvpr_lower_lev == wvpr_upper_lev) ) then

                    ! calculate the water vapour pressure for the station height using linear
                    ! interpolation method
                    meteo_stat(ind_POI) % wvpr= wvpr_lower_lev + (wvpr_upper_lev - wvpr_lower_lev) * (POIs_h_ell(ind_POI) - h_lower_lev) / (h_upper_lev - h_lower_lev) ! in [hPa]

                ! for all other cases use logarithmic interpolation
                else
                    
                    ! calculate the auxiliary parameter C for logarithmic interpolation
                    C= (h_upper_lev - h_lower_lev) / log(wvpr_upper_lev / wvpr_lower_lev)
                    
                    ! calculate the water vapour pressure for the station height using logarithmic
                    ! interpolation method
                    meteo_stat(ind_POI) % wvpr= wvpr_lower_lev * exp((POIs_h_ell(ind_POI) - h_lower_lev) / C) ! in [hPa]

                end if
                
                ! calculate the density of dry air for the station height level, [Md/R=1/Rd=Rd^-1]
                ! see scriptum Atmospheric Effects in Geodesy 2012, equation (2.53) on page 34
                rho_d_stat= (meteo_stat(ind_POI) % p - meteo_stat(ind_POI) % wvpr) / (Rd * meteo_stat(ind_POI) % T) ! in [hPa*s^2/m^2]
                ! calculate the density of water vapour for the station height level, [Mw/R=1/Rw=Rw^-1]
                ! see scriptum Atmospheric Effects in Geodesy 2012, equation (2.54) on page 34
                rho_w_stat= meteo_stat(ind_POI) % wvpr / (Rw * meteo_stat(ind_POI) % T) ! in [hPa*s^2/m^2]
                
                ! calculate inverse compressibility factor of water vapour for the station height level
                ! see Kleijer 2004, Troposphere Modeling and Filtering for Precise GPS Leveling, equation  (4.31) on page 28
                Zw_inv_stat= 1 + meteo_stat(ind_POI) % wvpr * (1 + 3.7d-4 * meteo_stat(ind_POI) % wvpr) * (-2.37321d-3 + 2.23366d0 / meteo_stat(ind_POI) % T - 710.792d0 / (meteo_stat(ind_POI) % T **2) + 7.75141d4 / (meteo_stat(ind_POI) % T **3))
                
                ! calculate the refractive index for the interpolation height level
                ! see Kleijer 2004, Troposphere Modeling and Filtering for Precise GPS Leveling, equation  (4.36) on page 29
                ! calculate the hydrostatic refractive index, [R/Md=Rd] (first step: only total rho is assigned)
                meteo_stat(ind_POI) % n_h= rho_d_stat + rho_w_stat
                ! calculate the wet refractive index (first step: only assign Zw_inv)
                meteo_stat(ind_POI) % n_w= Zw_inv_stat
                
            ! in case of suspending the ray-tracing for the current station
            else
                ! set variables for station where ray-tracing will be suspeded to "NaN"
                meteo_stat(ind_POI) % ind_lat1lon1= NaN
                meteo_stat(ind_POI) % ind_lat1lon2= NaN
                meteo_stat(ind_POI) % ind_lat2lon2= NaN
                meteo_stat(ind_POI) % ind_lat2lon1= NaN
                meteo_stat(ind_POI) % p= NaN
                meteo_stat(ind_POI) % T= NaN
                meteo_stat(ind_POI) % wvpr= NaN
                meteo_stat(ind_POI) % n_h= NaN
                meteo_stat(ind_POI) % n_w= NaN
                meteo_stat(ind_POI) % start_lev= NaN
                
                ! note: n_total will be automatically set to NaN, when it is calculated using the NaN values from n_h and n_w
                !meteo_stat(ind_POI) % n_total= NaN
                                          
            end if ! end of if for check of suspend ray-tracing
    
        end do ! end of loop over all stations
        
                                      
        ! Finish calculation of refractivities
        ! note: calculations are done for each station separate, i.e. the operations are executed on the scalars in the meteo_stat(index) structure for each index separately
        ! attention: in case of variables initialized as "NaN" in case that ray-tracing is skipped for a specific station "NaN" value will remain after the calculations!
        
        ! calculate the hydrostatic refractive index, [R/Md=Rd] (second step: multiplicate with constants and convert from refractivity to refractive index)
        meteo_stat % n_h= 1 + (meteo_stat % n_h * k1 * Rd) * 1d-6
        
        ! calculate the wet refractive index (second step: multiplicate assigned ZW_inv with rest of formula and convert from refractivity to refractive index)
        meteo_stat % n_w= 1 + ( meteo_stat % n_w * (k2s * meteo_stat % wvpr / meteo_stat % T + k3 * meteo_stat % wvpr / (meteo_stat % T **2))) * 1d-6
        
        ! calculate total refractive index
        meteo_stat % n_total= meteo_stat % n_h +  meteo_stat % n_w -1 ! -1 is necessary, because n_h_stat and n_w_stat are already converted to refractive indices and one +1 is two much if the two values are summed up!
        
    end subroutine calc_refr_ind_at_stations

end module module_calc_refr_ind_at_stations