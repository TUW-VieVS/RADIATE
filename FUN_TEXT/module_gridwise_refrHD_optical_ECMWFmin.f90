! module_gridwise_refrHD_optical_ECMWFmin.f90

!****************************************************************************
!
! PURPOSE:  Module containing subroutine to calculate the optical refractive index profiles for a grid
!            
!            "gridwise_refrHD_optical_ECMWFmin" is a subroutine to calculate the optical
!            refractive index profiles using the meteorologic values derived for example from a grib-file from the ECMWF.
!            The vertical resolution (height level resolution) from the ECMWF data will be increased and
!            above the highest ECMWF level extrapolation using a standard atmosphere will be carried out.
!            The lowest interpolation height level is determined according to the minimum height of all observing stations
!            in a specific epoch in order to avoid interpolation for levels that will never be used later during the ray-tracing.
!            The highest interpolation height level is defined (e.g. 84 000 m).
!            
!            The outputs of this function are the calculated refractive indices provided in an array.
!            
!            --- Version Info ---
!            This version of the function uses the mean height between the upper
!            and the lower level to determine the end where the lower height level ECMWF values are used for
!            interpolation of the pressure values. The other meteorological values are interpolated using
!            approaches that utilize the values from the lower and the upper level.
!            
!            The interpolation is done together for the whole horizontal grid.
!            Up to the minimum height of the ECMWF uppermost pressure level the ECMWF NWM data is used.
!            Above a standard atmosphere model is used for extrapolation of the meteorological parameters
!            and the extrapolation is again done for all grid points at one level at once.
!            
!            Two step approach of calculating n_h and n_w is the fastest method in MATLAB under minimum memory usage!
!            --------------------
!
!
! INPUT:
!         min_stat_h.... minimum ellipsoidal height in [m] of all stations for which ray-tracing will be done in a specific epoch
!         lat_grid...... (ellipsoidal) latitude of grid points in [°]
!         grid_size..... variable that determines the size in latitude and longitude of the
!                        grid [length_latitude, length_longitude]
!         nr_pres_lev... number of pressure levels retrieved from the ECMWF grib-file
!         gph........... geopotential heights of gridded points in [m]
!         p............. pressure in [hPa]
!         wvpr.......... water vapour pressure in [hPa]
!         T............. temperature in [K]
!         undulation.... geoid undulation for the grid points in [m]
! 
! 
! OUTPUT:
!         h_lev_lowest2upper_limit... interpolation/extrapolation height levels in [m] to increase the
!                                     vertical resolution of the meteorological parameters
!         nr_h_lev_lowest2upper_limit... total number of interpolated/extrapolated height levels
!         upper_limit_ECMWF........ uppermost height level with data from the ECMWF grib-file, in [m]
!         upper_limit_atm.......... uppermost height level assuming the end of the atmosphere, in [m]
!         p_int.................... interpolated/extrapolated pressure values in [hPa]
!         T_int.................... interpolated/extrapolated temperature values in [K]
!         wvpr_int................. interpolated/extrapolated water vapour pressure values in [hPa]
!         n_h_int.................. interpolated/extrapolated optical hydrostatic refractive index for the
!                                   height levels
!         n_w_int.................. interpolated/extrapolated optical wet refractive index for the height
!                                   levels
!         n_total_int.............. array containing the calculated total optical refractive index for all
!                                   interpolated/extrapolated height levels (with increased
!                                   resolution) covering the whole grid-area.
!
!----------------------------------------------------------------------------
! 
! Fortran-file created by Janina Boisits
!   based on Matlab and Fortran scripts created by Armin Hofmeister
! 
!----------------------------------------------------------------------------
! History:
! 
! 27.02.2018: create the Fortran-file based on the Fortran module "module_gridwise_refrHD_ECMWFmin.f90"
!             adapt module description
! 17.04.2018: add constants (wavelength, group refractive indices)
!             adapt equations for optical refractivity computation
!
!****************************************************************************

module module_gridwise_refrHD_optical_ECMWFmin

contains

    subroutine gridwise_refrHD_optical_ECMWFmin( min_stat_h, &
                                         lat_grid, &
                                         grid_size, &
                                         nr_pres_lev, &
                                         gph, &
                                         p, &
                                         wvpr, &
                                         T, &
                                         undulation, &
                                         h_lev_lowest2upper_limit, &
                                         nr_h_lev_lowest2upper_limit, &
                                         upper_limit_ECMWF, &
                                         upper_limit_atm, &
                                         p_int, &
                                         T_int, &
                                         wvpr_int, &
                                         n_h_int, &
                                         n_w_int, &
                                         n_total_int )

        ! Define modules to be used
        use module_constants, only: deg2rad, gn
        use module_meteo_constants
        use module_gph2horth
        use module_standard_atm
        use module_mean_total2D

        !----------------------------------------------------------------------------
    
        ! Variable definitions
        implicit none
    
    
        ! dummy arguments
        !----------------
    
        ! INPUT
        
        ! variable for input of minimum ellipsoidal height of all stations for which ray-tracing will be done in a specific epoch
        double precision, intent(in) :: min_stat_h
        
        ! variable for input of (ellipsoidal) latitude of grid points in [°]
        double precision, dimension(:, :), intent(in) :: lat_grid
    
        ! variable for input of grid dimension in latitude and longitude [size_latitude, size_longitude]
        integer, dimension(:), intent(in) :: grid_size
    
        ! variable for input of the number of pressure levels retrieved from the ECMWF grib-file
        integer, intent(in) :: nr_pres_lev
    
        ! variable for input of the geopotential heights of the gridded points in [m]
        double precision, dimension(:, :, :), intent(in) :: gph
    
        ! variable for input of the pressure in [hPa]
        double precision, dimension(:), intent(in) :: p
    
        ! variable for input of the water vapour pressure in [hPa]
        double precision, dimension(:, :, :), intent(in) :: wvpr
    
        ! variable for input of the temperature in [K]
        double precision, dimension(:, :, :), intent(in) :: T
    
        ! variable for input of the geoid undulation for the grid points in [m]
        double precision, dimension(:, :), intent(in) :: undulation
    
    
        ! OUTPUT
    
        ! variable for storing the interpolation/extrapolation height levels in [m]
        ! note: contains all height levels beginning with the lowest interpolation level up to the upper limit of the atmosphere
        double precision, dimension(:), allocatable, intent(out) :: h_lev_lowest2upper_limit
        
        ! variable for storing the number of total height levels (lowest_int_lev to upper_limit_atm)
        integer, intent(out) :: nr_h_lev_lowest2upper_limit
    
        ! variable for storing the uppermost height level with data from the ECMWF grib-file, in [m]
        double precision, intent(out) :: upper_limit_ECMWF
    
        ! variable for storing the upper limit of the atmosphere for the extrapolation of meteorologic values
        double precision, intent(out) :: upper_limit_atm
    
        ! variable for storing the interpolated/extrapolated gridded pressure values in [hPa]
        double precision, dimension(:, :, :), allocatable, intent(out) :: p_int
    
        ! variable for storing the interpolated/extrapolated gridded temperature values in [K]
        double precision, dimension(:, :, :), allocatable, intent(out) :: T_int
    
        ! variable for storing the interpolated/extrapolated gridded water vapour pressure values in [hPa]
        double precision, dimension(:, :, :), allocatable, intent(out) :: wvpr_int
    
        ! variable for storing the interpolated/extrapolated gridded hydrostatic refractive index
        double precision, dimension(:, :, :), allocatable, intent(out) :: n_h_int
    
        ! variable for storing the interpolated/extrapolated gridded wet refractive index
        double precision, dimension(:, :, :), allocatable, intent(out) :: n_w_int
    
        ! variable for storing the interpolated/extrapolated gridded total refractive index
        double precision, dimension(:, :, :), allocatable, intent(out) :: n_total_int
        
    
        ! local variables
        !----------------
    
        ! CONSTANTS
    
        ! Define constants and variables concerning meteorological fundamentals
        ! see "module_meteo_constants"
    
        ! see "module_constants" for deg2rad", "gn"
    
        ! most commonly used wavelength for optical observations and wave number
        double precision, parameter :: lambda = 0.532d0 ! in [mym]
        double precision, parameter :: sigma = 1 / lambda ! in [1/mym]
        
        ! group refractive index of dry air
        ! see Mendes & Pavlis 'High-accuracy zenith delay prediction at optical wavelengths'
        double precision, parameter :: Ngaxs = 1d-2 * (k1_opt * (k0_opt+sigma**2) / (k0_opt-sigma**2)**2 + k3_opt * (k2_opt+sigma**2) / (k2_opt-sigma**2)**2) * Cco2
        
        ! group refractive index of water vapour
        ! see Mendes & Pavlis 'High-accuracy zenith delay prediction at optical wavelengths'
        double precision, parameter :: Ngws = 1d-2 * cf * (w0 + 3 *w1*sigma**2 + 5*w2*sigma**4 + 7*w3*sigma**6)


        ! OTHER LOCAL VARIABLES
        
        ! variable for transformed input of (ellipsoidal) latitude of grid points from [°] to [rad]
        ! attention: The size is determined by the grid dimension through the input variable "grid_size"
        !            in order to avoid stack overflow in case of more than one array as following here allocatable option needs to be set!
        !            Side effect: deallocation after use is available.
        double precision, dimension(:, :), allocatable :: lat_grid_rad
    
        ! variable for storing the grid of the geopotential heights for the highest ECMWF level
        ! attention: The size is determined by the grid dimension through the input variable "grid_size"
        !            in order to avoid stack overflow in case of more than one array as following here allocatable option needs to be set!
        !            Side effect: deallocation after use is available.
        double precision, dimension(:, :), allocatable :: gph_max_ECMWF_lev
    
        ! variable for storing the grid of the orthometric heights for the highest ECMWF level
        !            in order to avoid stack overflow in case of more than one array as following here allocatable option needs to be set!
        !            Side effect: deallocation after use is available.
        double precision, dimension(:, :), allocatable :: h_orth_max_ECMWF_lev
    
        ! variable for storing the grid of the ellipsoidal heights for the highest ECMWF level
        !            in order to avoid stack overflow in case of more than one array as following here allocatable option needs to be set!
        !            Side effect: deallocation after use is available.
        double precision, dimension(:, :), allocatable :: h_ell_max_ECMWF_lev
    
        ! define varaible for the size of the variables for storing the step widths and the interval ends
        integer :: nr_intervals
    
        ! define variable for the general definition of possibilities for the different step widths for the height levels
        double precision, dimension(:), allocatable :: stepw_general
    
        ! define variable for the general definition of possibilities for the ends of the different height level intervals where the step widths change
        double precision, dimension(:), allocatable :: int_end_general
        
        ! define variable for storing logical vector of test which interval end heights that are above or equal to the minimum station height
        logical, dimension(:), allocatable :: log_above_intervals
        
        ! define variable for the number of found interval end heights that are above or equal to the minimum station height
        integer :: nr_above_intervals
        
        ! define variable for the values of the interval end heights that are above or equal to the minimum station height
        double precision, dimension(:), allocatable :: above_intervals
        
        ! define variable for the values of the step widths according to the interval end heights that are above or equal to the minimum station height
        double precision, dimension(:), allocatable :: above_stepws
        
        ! define variable for the index to the single interval end height value that is just above or equal to the minimum station height
        integer :: ind_int_end_above
        
        ! define variable for the value of the single interval end height value that is just above or equal to the minimum station height
        double precision :: int_end_above
        
        ! define variable for the value of the step width according to the single interval end height value that is just above or equal to the minimum station height
        double precision :: stepw_above
        
        ! define variable for the number of height steps which have to be subtracted to find the lowest interpolation height level from int_end_above
        integer ::nr_steps_subtr
        
        ! define variable for the lowest height level for which interpolation will be done
        double precision :: lowest_int_lev
    
        ! define variable for storing logical vector of test which step widths and end level heights are needed
        logical, dimension(:), allocatable :: needed_int_end
    
        ! define variable for the number of needed step widths and height levels where the step widths change (their number must be equal!)
        integer :: nr_needed_int_end
    
        ! define variable for the really needed different step widths for the height levels
        double precision, dimension(:), allocatable :: stepw_needed
    
        ! define variable for the really needed ends of the different height level intervals where the step widths change
        double precision, dimension(:), allocatable :: int_end_needed
    
        ! define variables for loop indices
        integer :: i, j
    
        ! define variable for storing temporal height level value
        double precision :: h_temp
    
        ! define variable for storing the number of height levels per step width
        integer, dimension(:), allocatable :: nr_lev_currstep
    
        ! define variable for index of current height level
        integer :: ind_lev_lastint
    
        ! define variable for storing logical vector of test which height levels are lower or equal to the ECMWF upper limit
        logical, dimension(:), allocatable :: check_lev_lowest2ECMWF
    
        ! define variable for storing the number of height levels from lowest level to ECMWF upper limit
        integer :: nr_h_lev_lowest2ECMWF
    
        ! define variable for interpolation loop index (pressure level index)
        integer :: level
    
        ! initialize index for interpolated height level
        integer :: ind_interp
        
        ! define variable for storing the values of part 1 of the gravity calculation for the pressure interpolation
        double precision, dimension(:, :), allocatable :: g_prep
        
        ! define variable for storing the values of the gravity for the pressure interpolation from the lower pressure level
        double precision, dimension(:, :), allocatable :: g_lower
        
        ! define variable for storing the values of the gravity for the pressure interpolation from the upper pressure level
        double precision, dimension(:, :), allocatable :: g_upper
    
        ! define variable for storing geopotential heights for lower level
        double precision, dimension(:, :), allocatable :: gph_lower_lev
    
        ! define variable for storing geopotential heights for upper level
        double precision, dimension(:, :), allocatable :: gph_upper_lev
    
        ! define variable for storing orthometric heights for lower level
        double precision, dimension(:, :), allocatable :: h_orth_lower_lev
    
        ! define variable for storing orthometric heights for upper level
        double precision, dimension(:, :), allocatable :: h_orth_upper_lev
    
        ! define variable for storing ellipsoidal heights for lower level
        double precision, dimension(:, :), allocatable :: h_ell_lower_lev
    
        ! define variable for storing ellipsoidal heights for upper level
        double precision, dimension(:, :), allocatable :: h_ell_upper_lev
    
        ! define variable for storing mean ellipsoidal height of the grid points in lower level
        double precision :: h_ell_lower_lev_mean
    
        ! define variable for storing mean ellipsoidal height of the grid points in upper level
        double precision :: h_ell_upper_lev_mean
    
        ! define variable for storing pressure of lower level
        double precision :: p_lower_lev
    
        ! define variable for storing  pressure of upper level
        double precision :: p_upper_lev
    
        ! define variable for storing water vapour pressure of lower level
        double precision, dimension(:, :), allocatable :: wvpr_lower_lev
    
        ! define variable for storing water vapour pressure of upper level
        double precision, dimension(:, :), allocatable :: wvpr_upper_lev
    
        ! define variable for storing temperature of lower level
        double precision, dimension(:, :), allocatable :: T_lower_lev
    
        ! define variable for storing temperature of upper level
        double precision, dimension(:, :), allocatable :: T_upper_lev
    
        ! define variable for storing virtual temperature of lower level
        double precision, dimension(:, :), allocatable :: Tv_lower_lev
    
        ! define variable for storing virtual temperature of upper level
        double precision, dimension(:, :), allocatable :: Tv_upper_lev
    
        ! define variable for storing mean ellipsoidal height of between mean upper and mean lower level
        double precision :: h_ell_mean_uplow
    
        ! define variable for storing the auxiliary parameters C for logarithmic interpolation of the water vapour pressure for the interpolation height level
        double precision, dimension(:, :), allocatable :: C
    
        ! define variable for storing the density of dry air for the interpolation height level
        !double precision, dimension(:, :), allocatable :: rho_d_int
    
        ! define variable for storing the density of water vapour for the interpolation height level
        double precision, dimension(:, :), allocatable :: rho_w_int
    
        ! define variable for storing the total density of air for the interpolation height level
        double precision, dimension(:, :), allocatable :: rho_total_int
    
        ! define variable for storing the temperature in degree celsius for the interpolation height level
        double precision, dimension(:, :), allocatable :: t_int_deg
    
        ! define variable for storing the ratio of water vapour pressure and total pressure for the interpolation height level
        double precision, dimension(:, :), allocatable :: xw
    
        ! define variable for storing the inverse compressibility factor of water vapour for the interpolation height level
        !double precision, dimension(:, :), allocatable :: Zw_inv_int
    
        ! define variable for storing the compressibility factor of air for the interpolation height level
        double precision, dimension(:, :), allocatable :: Z_int
    
        ! define variable for storing extrapolated temperature value
        double precision :: T_extrap
    
        ! define variable for storing extrapolated pressure value
        double precision :: p_extrap
    
        ! define variable for storing extrapolated water vapour pressure value
        double precision :: wvpr_extrap
    
        !----------------------------------------------------------------------------
    
        !============================================================================
        ! Body of the subroutine
        !============================================================================
    
        
        ! Define limits for the atmosphere model
    
        ! upper limit of the atmosphere for the extrapolation of meteorologic values
        ! note: d0 at the end of all numbers is necessary to receive double precision and not integer type!
        upper_limit_atm = 84000d0 ! in [m], must be higher than 36000 m otherwise possible error source when defining height levels (see below)
                                  ! must be max. 84852 m otherwise error in function "standard_atm"
                                  ! wvpr for heights above are set to 0 --> be careful when defining ECMWF limit
        
        ! check if the minimum station height is larger than upper_limit_atm
        ! note: In this case the program is stopped since the vertical interpolation can not be carried out.
        !       When calculating the meteorological data at the specific station heights each station height will be checked for being able to do
        !       ray-tracing in case the individual heights are too high
        if ( min_stat_h > upper_limit_atm ) then
            ! report error message
            write(unit= *, fmt= '(a)') 'Error: Minimum station height is too large for doing vertical interpolation! Program stopped!'
            ! stop the program
            stop
        end if
        
        ! determine upper limit of the ECMWF profiles
        ! --> calculate the minimum ellipsoidal height of the grid points in the uppermost ECMWF level to
        ! define the ECMWF limit used for defining the interpolation heights and for setting the limit up to which the
        ! ECMWF data is used for interpolation
    
        ! allocate variables
        allocate( lat_grid_rad(grid_size(1), grid_size(2)), &
                  gph_max_ECMWF_lev(grid_size(1), grid_size(2)), &
                  h_orth_max_ECMWF_lev(grid_size(1), grid_size(2)), &
                  h_ell_max_ECMWF_lev(grid_size(1), grid_size(2)) )
    
        
        ! transform latitude grid-values from [°] to [rad]
        lat_grid_rad= lat_grid * deg2rad ! in [rad]
        
    
        ! get geopotential heights for the highest ECMWF level
        gph_max_ECMWF_lev= gph(nr_pres_lev, :, :) ! not necessary in Fortran: a call of a function like squeeze() = variant application of reshape() in MATLAB, which removes singleton dimensions (here the pressure level dimension)
    
    
        ! transform geopotential (=dynamic) heights to orthometric heights
        ! It is necessary to transform the geopotential heights retrieved from the ECMWF
        ! geopotentials of the grib-file in a first step to orthometric heights. Then using the geoid
        ! undulation retrieved from the EGM2008 ellipsoidal (=geometric) heights can be calculated. The
        ! ellipsoidal heights define the reference system used for ray tracing.
     
        ! This section uses an approximation of the conversion between geopotential and orthometric
        ! (=geometric) heights to avoid an iterative calculation
        ! Therefore the subroutine "gph2horth" is called, which uses a formula from Kraus for the
        ! conversion in one step.
        ! calculate orthometric heights of the gridded points in the highest ECMWF level
        ! --> h_orth_max_ECMWF_lev in [m], using gph in [m] and lat_grid in [°]
        call gph2horth( gph_max_ECMWF_lev, &
                        lat_grid, &
                        h_orth_max_ECMWF_lev )
    
        ! transform orthometric heights to ellipsoidal (=geometric) heights
        h_ell_max_ECMWF_lev= h_orth_max_ECMWF_lev + undulation ! in [m]
    
        ! calculate the minimum ellipsoidal height of the grid points in the uppermost ECMWF level to define the
        ! ECMWF limit used for defining the interpolation heights and for setting the limit up to which the
        ! ECMWF data is used for interpolation
        upper_limit_ECMWF= minval(h_ell_max_ECMWF_lev) ! in [m], result: omitting optional dim parameter leads to minimum of whole array
    
        ! deallocate variables only used for getting minimum ellipsoidal height of the grid points in the uppermost ECMWF level,
        ! which are not needed any more
        deallocate( gph_max_ECMWF_lev, &
                    h_orth_max_ECMWF_lev, &
                    h_ell_max_ECMWF_lev )
    
    
        ! Define height levels for refinement of vertical resolution (interpolation and extraploation steps)
        
        ! Define step widths for the height levels used for interpolation and extrapolation
        ! see declaration part
    
        ! the following tables show the step widths for the height levels which are used to increase the
        ! vertical resolution of the meteorological parameters
    
        ! step widths, see Rocken et al. 2001, Improved mapping of tropospheric delays:
        ! A:  10 m: [    0m or below,  2000m]
        ! B:  20 m: ] 2000m,           6000m]
        ! C:  50 m: ] 6000m,          16000m]
        ! D: 100 m: ]16000m,          36000m]
        ! E: 500 m: ]36000m, upper_limit_atm]
    
        ! define size of the variables for storing the step widths and the interval ends
        nr_intervals= 5
        
        ! allocate variables for general possible step widths and height levels were step width changes
        ! attention: stepw_general and int_end_general must have the same size!
        allocate( stepw_general(nr_intervals), &
                  int_end_general(nr_intervals) )
    
        ! define variables for the different step widths
        stepw_general= [ 10d0, 20d0, 50d0, 100d0, 500d0 ]
        
        ! define variables for the end of the different height level intervals where the step widths change
        ! attention: The values in int_end_general must be strictly ascending otherwise the program will fail!
        ! note: The last entry must be upper_limit_atm.
        int_end_general= [ 2000d0, 6000d0, 16000d0, 36000d0, upper_limit_atm ]
        
        
        ! define lowest height level in [m] for which interpolation will be done according to the input of the
        ! minimum ellipsoidal height of the stations that are observing at the current epoch of the grib-data
        ! note: The lowest interpolation level is the interpolation level that is the last below or equal to
        !       the minimum station height, which fullfills the following condition:
        !           The desired interpolation level must be a level that can be constructed by starting at the
        !           one interval end boundary that is just above or equal to the minimum station height and 
        !           subsequent subtracting of the specific height steps of the interval until the height is
        !           just below or equal to the minimum station height.
        !           This setup ensures that independent of the minimum station height in all epochs of grib-data
        !           the exact same interpolation height levels are used.
        
        
        ! find the one interval end, which is just above or equal to the minimum station height
                
        ! allocate variable (logical vector of same size as stew_general and int_end_general)
        allocate(log_above_intervals(nr_intervals))
        
        ! determine a logical vector, which has the values .TRUE. in case the minimum station height is lower than or equal to the interval ends
        log_above_intervals= ( int_end_general >= min_stat_h )
    
        ! determine the number of interval ends that are larger than or equal to min_stat_h
        nr_above_intervals= count(log_above_intervals)
    
        ! allocate the variables containing the interval ends that are above the minimum station height and get the according step widths for these intervals
        ! note: above_intervals and above_stepw must have the same size!
        allocate( above_intervals(nr_above_intervals), &
                  above_stepws(nr_above_intervals) )
    
        ! find and assign the interval ends that are larger than or equal to min_stat_h and the according step widths for these intervals
        ! note: the mask of log_above_intervals leads to assigning only the .TRUE. elements to the result
        above_intervals= pack(int_end_general, log_above_intervals)
        above_stepws= pack(stepw_general, log_above_intervals)
        
        ! get the index of the lowest interval end that is above or equal to the minimum station height
        ! note: minloc(array, dim) returns the index of the first appearance of the minimum value in array order of the dimension dim.
        ind_int_end_above= minloc(above_intervals, 1)
        
        ! get the lowest interval end that is above or equal to the minimum station height and the according step width for this interval
        int_end_above= above_intervals(ind_int_end_above)
        stepw_above= above_stepws(ind_int_end_above)
        
        
        ! determine the the lowest interpolation level
        ! by subtracting the number of height steps, so that lowest_int_lev will just be lower than or equal to the minimum station height
        
        ! determine the number of steps which need to be subtracted in order to reach lowest_int_lev
        ! note: Use ceiling() to get the number of steps that are just enough to get below or equal to the minimum station height as
        !       this returns the smallest integer larger or equal to the input.
        ! note: In case the interval end is equal to the minimum station height the resulting number of subtracted steps is 0.0 and no problem will occur.
        !       In case of numeric reasons and therefore not being exactly 0.0 the lowest interpolation height level will be set one height level step too low,
        !       which will not cause any problems. This level will then just be unnecessary in the later ray-tracing.
        nr_steps_subtr= ceiling( (int_end_above - min_stat_h ) / stepw_above )
        
        
        ! determine the lowest interpolation level by subtracting the number of height steps from the interval end above or equal to the minimum station height
        lowest_int_lev= int_end_above - nr_steps_subtr * stepw_above ! in [m]
        
    
        ! define all height levels used for interpolation and extrapolation starting at the lowest interpolation height level
        ! note: In general the following determination of the needed interval ends and step widths will deliver the same results
        !       as the detemined values for the definition of lowest_int_lev. But in case of a minimum station height that is equal to an
        !       interval end, this is not the case since only the interval ends that are really above the lowest_int_lev are needed.
        !       Therefore the needed intervals and step widths are determined explicitely for lowest_int_lev.
    
        ! determine which step widths and end level heights are needed
        ! allocate variable (logical vector of same size as stew_general and int_end_general)
        allocate(needed_int_end(nr_intervals))
        needed_int_end= int_end_general > lowest_int_lev
    
        ! determine the number of needed step widths and height levels where the step widths change (their number must be equal!)
        nr_needed_int_end= count(needed_int_end)
    
        ! allocate the variables containing the really needed step widths and end level heights
        ! note: stepw_needed and int_end_needed must have the same size!
        allocate( stepw_needed(nr_needed_int_end), &
                  int_end_needed(nr_needed_int_end))
    
        ! find and assign the needed step widths and interval end levels used for setting up the interpolation height levels
        ! note: the mask of needed_int_end leads to assigning only the .TRUE. elements to the result
        stepw_needed= pack(stepw_general, needed_int_end)
        int_end_needed= pack(int_end_general, needed_int_end)
    
    
        ! determine number of total height levels = size of the variable for all height levels beginning with the lowest interpolation level
    
        ! initialize temporal variable for storing current height level
        ! set it as lowest height level
        h_temp= lowest_int_lev
    
        ! initialize the counter of the total height levels
        nr_h_lev_lowest2upper_limit= 1
    
        ! allocate the variable that stores the number of height levels per step width
        allocate( nr_lev_currstep(nr_needed_int_end) )
    
        ! loop over individual needed step widths to calculate the number of height levels per step width
        do i= 1, nr_needed_int_end
        
            ! calculate the number of height levels for current step width
            ! note: ceiling() delivers type integer, integer number greater or equal to argument
            nr_lev_currstep(i)= ceiling( (int_end_needed(i) - h_temp) / stepw_needed(i) )
        
            ! calculate last height level of current step width (real end of step interval)
            h_temp= h_temp + nr_lev_currstep(i) * stepw_needed(i)
        
            ! sum up the number of total height levels (add number of height levels for currewt step width)
            nr_h_lev_lowest2upper_limit= nr_h_lev_lowest2upper_limit + nr_lev_currstep(i)
        
        end do
    
        ! allocate variable for all height levels beginning with the lowest interpolation level
        ! using the above determined total number of height levels
        allocate( h_lev_lowest2upper_limit(nr_h_lev_lowest2upper_limit) )
    
        ! set first (lowest height) level
        h_lev_lowest2upper_limit(1)= lowest_int_lev
    
        ! initialize index of current height level
        ind_lev_lastint= 1
    
        ! determine values of the different height levels
        ! loop over individual needed step widths to calculate the height levels
        do i= 1, nr_needed_int_end
        
            ! calculate height level
            h_lev_lowest2upper_limit( ind_lev_lastint + 1 : ind_lev_lastint + nr_lev_currstep(i) )= [( h_lev_lowest2upper_limit( ind_lev_lastint ) + stepw_needed(i) * j, j= 1, nr_lev_currstep(i) )]
        
            ! determine index of last already assigned height level
            ind_lev_lastint= ind_lev_lastint + nr_lev_currstep(i)
        
        end do
    
    
        ! define number of height levels in specific intervals
    
        ! determine number of height levels from lowest interpolation level up to the ECMWF upper limit
        ! allocate variable (logical vector of size= nr_h_lev_lowest2upper_limit)
        allocate(check_lev_lowest2ECMWF(nr_h_lev_lowest2upper_limit))
        check_lev_lowest2ECMWF= ( h_lev_lowest2upper_limit <= upper_limit_ECMWF )
    
        ! determine number of height levels from lowest interpolation level up to the ECMWF upper limit
        nr_h_lev_lowest2ECMWF= count(check_lev_lowest2ECMWF)
    
    
        ! Initialize all variables for which interpolation (vertical refinement) will be done and that will be stored for all height levels
    
        ! total number of height levels for the interpolation are defined through the variable "h_lev_lowest2upper_limit"
    
        ! allocate all necessary variables with their dimensions (3-dimensional variables)
        allocate( p_int(nr_h_lev_lowest2upper_limit, grid_size(1), grid_size(2)), & ! interpolated pressure
                  T_int(nr_h_lev_lowest2upper_limit, grid_size(1), grid_size(2)), & ! interpolated temperature
                  wvpr_int(nr_h_lev_lowest2upper_limit, grid_size(1), grid_size(2)), & ! interpolated water vapour pressure
                  n_h_int(nr_h_lev_lowest2upper_limit, grid_size(1), grid_size(2)), & ! hydrostatic refractive index
                  n_w_int(nr_h_lev_lowest2upper_limit, grid_size(1), grid_size(2)), & ! wet refractive index
                  n_total_int(nr_h_lev_lowest2upper_limit, grid_size(1), grid_size(2)) ) ! total refractive index
    
        ! set all values to NaN
        ! note: this is not needed and slow
    
        ! allocate all necessary variables with their dimensions (2-dimensional variables)
        allocate( g_prep(grid_size(1), grid_size(2)), & ! variable for part 1 of gravity for pressure interpolation
                  g_lower(grid_size(1), grid_size(2)), & ! variable for gravity for the pressure interpolation from the lower pressure level
                  g_upper(grid_size(1), grid_size(2)), & ! vvariable for gravity for the pressure interpolation from the upper pressure level
                  gph_lower_lev(grid_size(1), grid_size(2)), & ! variable for geopotential heights of lower level
                  gph_upper_lev(grid_size(1), grid_size(2)), & ! variable for geopotential heights of upper level
                  h_orth_lower_lev(grid_size(1), grid_size(2)), & ! variable for orthometric heights of lower level
                  h_orth_upper_lev(grid_size(1), grid_size(2)), & ! variable for orthometric heights of upper level
                  h_ell_lower_lev(grid_size(1), grid_size(2)), & ! variable for ellipsoidal heights of lower level
                  h_ell_upper_lev(grid_size(1), grid_size(2)), & ! variable for ellipsoidal heights of upper level
                  wvpr_lower_lev(grid_size(1), grid_size(2)), & ! variable for water vapour pressure of lower level
                  wvpr_upper_lev(grid_size(1), grid_size(2)), & ! variable for water vapour pressure of upper level
                  T_lower_lev(grid_size(1), grid_size(2)), & ! variable for temperature of lower level
                  T_upper_lev(grid_size(1), grid_size(2)), & ! variable for temperature of upper level
                  Tv_lower_lev(grid_size(1), grid_size(2)), & ! variable for virtual temperature of lower level
                  Tv_upper_lev(grid_size(1), grid_size(2)), & ! variable for virtual temperature of upper level
                  C(grid_size(1), grid_size(2)), & ! variable for auxiliary parameter C for logarithmic interpolation of the water vapour pressure for the interpolation height level
                  !rho_d_int(grid_size(1), grid_size(2)), & ! variable for density of dry air for the interpolation height level
                  rho_w_int(grid_size(1), grid_size(2)), & ! variable for density of water vapour for the interpolation height level
                  rho_total_int(grid_size(1), grid_size(2)), & ! variable for total density of wet air for the interpolation height level
                  t_int_deg(grid_size(1), grid_size(2)), & ! variable for temperature in degree celsius for the interpolation height level
                  xw(grid_size(1), grid_size(2)), & ! variable for the ratio of water vapour pressure and total pressure for the interpolation height level
                  !Zw_inv_int(grid_size(1), grid_size(2)) ) ! variable for inverse compressibility factor of water vapour for the interpolation height level
                  Z_int(grid_size(1), grid_size(2)) ) ! variable for the compressibility factor of wet air for the interpolation height level
    
    
        ! Refinement of the vertical profiles from lowest interpolation level up to the ECMWF upper limit
        ! This section is necessary to improve the vertical resolution of the meteorological data received
        ! from the ECMWF grib-file (downwards to lowest interpolation level, but also upwards to the ECMWF upper limit).
        ! Therefore an interpolation of the meteorological values is carried out using the previously
        ! defined height levels for which the values should be interpolated (extrapolated).
     
        ! Total pressure, water vapor pressure and temperature have non-linear contributions in the formula 
        ! for calculation of the refractive index. Therefore it wouldn't be suitable to interpolate refractive index
        ! values directly, which are needed for ray-tracing.
        ! Instead each contributing meteorological parameter has to be interpolated on his own. Afterwards
        ! refractive index can be calculated with the interpolated values.
    
        ! calculate gravitational acceleration for pressure interpolation (first part)
        ! see scriptum Atmospheric Effects in Geodesy 2012, equation (1.64) on page 14
        
        ! prepare gravity calculation
        ! first part: without height dependent term
        ! see "module_constants" for normal gravity constant "gn"
        ! note: since Fortran 90: two arrays that are conformable (same shape) are processed element-wise in a common expression
        !                         array and scalar in a statement lead to element-wise calculation of array elements
        g_prep= gn * ( 1 - 0.0026373d0 * cos( 2 * lat_grid_rad ) + 0.0000059d0 * ( cos( 2 * lat_grid_rad ) )**2 )
        
        ! initialize index for interpolation
        ind_interp= 1
    
        ! treatment of all pressure levels in a loop in order to do interpolation
        ! the following loop uses the data received from the grib-file
        ! pressure levels retrieved from the grib-file have already been and should be sorted in such a way that the lowest height = max. pressure as first entry
        interpolation_loop: do level= 1, nr_pres_lev-1 ! only till nr_pres_lev-1, because level+1 is used in the loop for the upper level
        
            ! get geopotential heights for the current = lower level
            !note: not necessary in Fortran: a call of a function like squeeze() = variant application of reshape() in MATLAB, which removes singleton dimensions (here the pressure level dimension)
            gph_lower_lev= gph(level, :, :)
            ! get geopotential heights for the upper level
            gph_upper_lev= gph(level+1, :, :)
        
            ! Transform geopotential (=dynamic) heights to orthometric heights
            ! It is necessary to transform the geopotential heights retrieved from the ECMWF
            ! geopotentials in the grib-file in a first step to orthometric heights. Then using the geoid
            ! undulation retrieved from the EGM2008 ellipsoidal (=geometric) heights can be calculated. The
            ! ellipsoidal heights define the reference system used for ray tracing.
        
            ! This section uses an approximation of the conversion between geopotential and orthometric
            ! (=geometric) heights to avoid an iterative calculation
            ! Therefore the subroutine "gph2horth" is called which uses a formula from Kraus for the
            ! conversion in one step.
        
            ! calculate orthometric heights of the gridded points of the lower level
            ! --> h_orth_lower_lev in [m], using gph in [m] and lat_grid in [°]
            call gph2horth(gph_lower_lev, &
                           lat_grid, &
                           h_orth_lower_lev)
        
            ! calculate orthometric heights of the gridded points of the upper level
            ! --> h_orth_upper_lev in [m], using gph in [m] and lat_grid in [°]
            call gph2horth(gph_upper_lev, &
                           lat_grid, &
                           h_orth_upper_lev)
        
            ! Transform orthometric heights to ellipsoidal (=geometric) heights
            ! for the lower level
            h_ell_lower_lev=h_orth_lower_lev + undulation ! in [m]
            ! for the upper level
            h_ell_upper_lev=h_orth_upper_lev + undulation ! in [m]
        
        
            ! Calculate the mean ellipsoidal height of the grid points in the specific level
            ! for the lower level
            call mean_total2D( h_ell_lower_lev, h_ell_lower_lev_mean ) ! in [m], result: scalar (two times mean --> mean in both directions)
            ! for the upper level
            call mean_total2D( h_ell_upper_lev, h_ell_upper_lev_mean ) ! in [m], result: scalar (two times mean --> mean in both directions)
        
        
            ! get the meteorological values received from the grib-file for the specific lower and upper level
        
            ! get lower level pressure
            p_lower_lev= p(level) ! in [hPa]
            ! get upper level pressure
            p_upper_lev= p(level+1) ! in [hPa]
        
            ! get lower level water vapour pressure
            wvpr_lower_lev= wvpr(level, :, :) ! in [hPa]
            ! get upper level water vapour pressure
            wvpr_upper_lev= wvpr(level+1, :, :) ! in [hPa]
        
            ! get lower level temperature
            T_lower_lev= T(level, :, :) ! in [K]
            ! get upper level temperature
            T_upper_lev= T(level+1, :, :) ! in [K]
        
        
            ! calculate the virtual temperature using meteorological values from the lower level
            ! see scriptum Atmospheric Effects in Geodesy 2012, equation (1.67) on page 15 or
            ! see Nafisi et al. 2012, Ray-traced tropospheric delays in VLBI analysis , equation (21) on page 7
            ! + Kleijer 2004 equation (4.13) + (4.3))
            ! note: since Fortran 90: two arrays that are conformable (same shape) are processed element-wise in a common expression
            !                         array and scalar in a statement lead to element-wise calculation of array elements
            Tv_lower_lev= T_lower_lev * p_lower_lev / (p_lower_lev - (1 - Mw / Md) * wvpr_lower_lev) ! in [K]
        
            ! calculate the virtual temperature using meteorological values from the upper level
            ! see scriptum Atmospheric Effects in Geodesy 2012, equation (1.67) on page 15 or
            ! see Nafisi et al. 2012, Ray-traced tropospheric delays in VLBI analysis , equation (21) on page 7
            ! + Kleijer 2004 equation (4.13) + (4.3))
            ! note: since Fortran 90: two arrays that are conformable (same shape) are processed element-wise in a common expression
            !                         array and scalar in a statement lead to element-wise calculation of array elements
            Tv_upper_lev= T_upper_lev * p_upper_lev / (p_upper_lev - (1 - Mw / Md) * wvpr_upper_lev) ! in [K]
        
            
            ! calculate gravitational acceleration for pressure interpolation (second part)
            ! see scriptum Atmospheric Effects in Geodesy 2012, equation (1.64) on page 14
            
            ! prepare gravity calculation
            ! second part: apply height dependent term
            ! note: since Fortran 90: two arrays that are conformable (same shape) are processed element-wise in a common expression
            !                         array and scalar in a statement lead to element-wise calculation of array elements
            ! for use of pressure interpolation with the lower pressure level
            g_lower= g_prep * (1-3.14d-7 * h_ell_lower_lev) ! in [m/s^2]
            ! for use of pressure interpolation with the upper pressure level
            g_upper= g_prep * (1-3.14d-7 * h_ell_upper_lev) ! in [m/s^2]
    
        
            ! Calculate the mean ellipsoidal height of between upper and lower level
            h_ell_mean_uplow= (h_ell_lower_lev_mean + h_ell_upper_lev_mean) / 2
        
        
            ! refinement of meteorological values received from the grib-file
        
            ! start loop for interpolation
            ! --> until the mean upper level of the grib-file is reached interpolation is done for all height levels below (defined in the variable "h_lev_lowest2upper_limit")
            inner_interpolation_loop: do while (h_lev_lowest2upper_limit(ind_interp) <= h_ell_upper_lev_mean) ! last for-loop delivers mean value of last ECMWF pressure level as value for h_ell_upper_lev_mean
                                                                                                        ! --> upper_limit_ECMWF is the minimum and therefore always lower or equal
                ! --> if the actual interpolation height level is lower or equal than the mean height level between the current upper and lower level of the grib-file
                ! --> interpolation for the pressure is done using the pressure and virtual temperature values from the lower level 
                if (h_lev_lowest2upper_limit(ind_interp) <= h_ell_mean_uplow) then
                
                    ! interpolate the pressure --> linear interpolation of the logarithm-transformed pressure values (=exponential interpolation for the pressure domain)
                    ! see Nafisi et al. 2012, Ray-traced tropospheric delays in VLBI analysis , equation (20) on page 7
                    ! note: since Fortran 90: two arrays that are conformable (same shape) are processed element-wise in a common expression
                    !                         array and scalar in a statement lead to element-wise calculation of array elements
                    p_int(ind_interp, :, :)= p_lower_lev * exp(-(h_lev_lowest2upper_limit(ind_interp) - h_ell_lower_lev) * g_lower / (Rd * Tv_lower_lev)) ! in [hPa]
                
                ! --> if the actual interpolation height level is higher than the mean height level between the current upper and lower level of the grib-file
                ! --> interpolation for the pressure is done using the pressure and virtual temperature values from the upper level
                else
                    ! interpolate the pressure --> linear interpolation of the logarith-transformed pressure values (=exponential interpolation for the pressure domain)
                    ! see Nafisi et al. 2012, Ray-traced tropospheric delays in VLBI analysis , equation (20) on page 7
                    ! note: since Fortran 90: two arrays that are conformable (same shape) are processed element-wise in a common expression
                    !                         array and scalar in a statement lead to element-wise calculation of array elements
                    p_int(ind_interp, :, :)= p_upper_lev * exp(-(h_lev_lowest2upper_limit(ind_interp) - h_ell_upper_lev) * g_upper / (Rd * Tv_upper_lev)) ! in [hPa]
            
                end if
            
                ! interpolate the temperature --> linear interpolation
                T_int(ind_interp, :, :)= T_upper_lev + (T_lower_lev - T_upper_lev) / (h_ell_lower_lev - h_ell_upper_lev) * (h_lev_lowest2upper_limit(ind_interp) - h_ell_upper_lev)
            
                ! for testing switch of pressure levels
                !print *, level, ind_interp
            
            
                ! interpolate the water vapour pressure
                ! --> logarithmic interpolation, if possible otherwise use linear interpolation
                ! for logarithmic interpolation equations see
                ! Nafisi et al. 2012, Ray-traced tropospheric delays in VLBI analysis , equation (22) on page 7
                ! or Böhm and Schuh 2003, Vienna mapping functions
            
                ! check if the auxiliary parameter C for logarithmic interpolation can be calculated correctly
                ! in case that wvpr_lower_lev or wvpr_upper_lev are 0 or if wvpr_lower_lev is equal to wvpr_upper_lev
                ! C would deliver 0 or Inf --> no correct logarithmic interpolation possible
                ! --> switch to linear interpolation in these cases
            
                ! set C to NaN
                ! note: this is not necessary as C values will be calculated and taken only for those elements where needed for logarithmic interpolation of wvpr
                !       and setting to NaN is slow
            
                ! check if linear interpolation is necessary
                where ( (wvpr_lower_lev == 0) .OR. (wvpr_upper_lev == 0) .OR. (wvpr_lower_lev == wvpr_upper_lev) )
                
                    ! calculate the water vapour pressure for the interpolation height level using linear
                    ! interpolation method
                    ! note: this overwrites values from prevoius loop (only at positions where calculation takes place, remaining positions will be overwritten when calculating logarithmic interpolated values)
                    wvpr_int(ind_interp, :, :)= wvpr_lower_lev + (wvpr_upper_lev - wvpr_lower_lev) * (h_lev_lowest2upper_limit(ind_interp) - h_ell_lower_lev) / (h_ell_upper_lev - h_ell_lower_lev) ! in [hPa]
            
                ! for all other cases use logarithmic interpolation
                elsewhere
                
                    ! calculate the auxiliary parameter C for logarithmic interpolation
                    C= (h_ell_upper_lev - h_ell_lower_lev) / log(wvpr_upper_lev / wvpr_lower_lev)
                
                    ! calculate the water vapour pressure for the interpolation height level using logarithmic
                    ! interpolation method
                    ! note: this overwrites values from prevoius loop (only at positions where calculation takes place, remaining positions will be overwritten when calculating linear interpolated values)
                    wvpr_int(ind_interp, :, :)= wvpr_lower_lev * exp( (h_lev_lowest2upper_limit(ind_interp) - h_ell_lower_lev) / C) ! in [hPa]
                
                end where
            
            
                ! calculate the density of dry air for the interpolation height level, [Md/R=1/Rd=Rd^-1]
                ! see scriptum Atmospheric Effects in Geodesy 2012, equation (2.53) on page 34
                !rho_d_int= ( p_int(ind_interp, :, :) - wvpr_int(ind_interp, :, :) ) / ( Rd * T_int(ind_interp, :, :) ) ! in [hPa*s^2/m^2]

                ! calculate auxiliary paramerters
                ! see Mendes & Pavlis 'High-accuracy zenith delay prediction at optical wavelengths'
                t_int_deg = T_int(ind_interp, :, :) - 273.15d0 ! in [°C]
                xw = wvpr_int(ind_interp, :, :) / p_int(ind_interp, :, :)
                Z_int = 1 - (p_int(ind_interp, :, :) / T_int(ind_interp, :, :)) * (a0 + a1*t_int_deg + a2*t_int_deg**2 + (b0 + b1*t_int_deg) * xw + (c0 + c1*t_int_deg) * xw**2) + (p_int(ind_interp, :, :) / T_int(ind_interp, :, :))**2 * (d0 + e0*xw**2)
                        
                
                ! calculate the density of water vapour for the interpolation height level, [Mw/R=1/Rw=Rw^-1]
                ! see scriptum Atmospheric Effects in Geodesy 2012, equation (2.54) on page 34
                rho_w_int = (Mw * xw * p_int(ind_interp, :, :)) / (R * Z_int * T_int(ind_interp, :, :)) ! in [hPa*s^2/m^2]
            
                ! calculate the total density of air for the interpolation height level
                ! see Mendes & Pavlis 'High-accuracy zenith delay prediction at optical wavelengths'
                rho_total_int = Md / (R * Z_int) * (p_int(ind_interp, :, :) / T_int(ind_interp, :, :) - (1 - Mw/Md) * wvpr_int(ind_interp, :, :) / T_int(ind_interp, :, :)) ! in [hPa*s^2/m^2]
            
                ! calculate inverse compressibility factor of water vapour for the interpolation height level
                ! see Kleijer 2004, Troposphere Modeling and Filtering for Precise GPS Leveling, equation  (4.31) on page 28
                !Zw_inv_int= 1 + wvpr_int(ind_interp, :, :) * ( 1 + 3.7d-4 * wvpr_int(ind_interp, :, :) ) * ( -2.37321d-3 + 2.23366d0 / T_int(ind_interp, :, :) - 710.792d0 / ( T_int(ind_interp, :, :) **2 ) + 7.75141d4 / ( T_int(ind_interp, :, :) **3 ) )
            
                ! calculate the refractive index for the interpolation height level
                ! see Mendes & Pavlis 'High-accuracy zenith delay prediction at optical wavelengths'
                ! calculate the hydrostatic refractive index (first step: only assign non constant terms)
                n_h_int(ind_interp, :, :)= rho_total_int
                ! calculate the wet refractive index
                n_w_int(ind_interp, :, :)= Ngws * (rho_w_int/rho_ws) - Ngaxs * (Td_k/Pd) * (Mw/Md) * (Zd/Z_int) * (wvpr_int(ind_interp, :, :) / T_int(ind_interp, :, :))
            
                ! raise index (for interpolation step)
                ind_interp=ind_interp+1;
            
                ! Stop the while loop if the next interpolation height level would be above
                ! the last entry in "h_lev_lowest2ECMWF".
                ! Otherwise an index error in the next call of the while statement will occur.
                ! note: This if expression must contain the question if the variable "ind_interp" is higher than a possible index, otherwise an exception error would occur
                !       as an index higher than possible would be used in the next loop cycle!
                !       In the specific case here the if also asks if the next height level would exceed the index number of the last level that is below or equal to the ECMWF upper limit.
                !       So two things can be checked here at once, which is necessary for all possible error sources in the while-loop.
                if (ind_interp > nr_h_lev_lowest2ECMWF) then
                    exit ! continue execution after end do of current do loop (= after inner_interpolation_loop)
                end if
            
            end do inner_interpolation_loop
            
            ! Stop the interpolation loop over the ECMWF pressure height levels if the next
            ! interpolation height level would be above the last entry in "h_lev_lowest2ECMWF".
            ! Otherwise an index error in the while statement may occur in case a new for-loop cycle
            ! is started.
            ! This case can happen if the used last entry in "h_lev_lowest2ECMWF" that is just
            ! below the minimum ECMWF height level is reached at a point where the for-loop is below
            ! its last cycle which is the last but one pressure level of the ECMWF.
            if (ind_interp > nr_h_lev_lowest2ECMWF) then
                exit ! continue execution after end do of current do loop (= after interpolation_loop)
            end if
        
        end do interpolation_loop
    
    
    
        ! Extrapolation of the vertical profiles above the ECMWF upper limit to the upper limit of the atmosphere
        ! This section uses the subroutine "standard_atm" to calculate temperature and pressure values for the
        ! height levels above the upper height limit of the ECMWF data.
    
        ! continue interpolation for height levels that are above the defined ECMWF upper limit
        ! (minimum height of all grid points in the uppermost ECMWF pressure level)
    
        ! loop over all height levels in "h_lev_lowest2upper_limit" that are above the ECMWF limit (decided by index for height level)
        ! in order to call the subroutine "standard_atm" for determining temperature and pressure
        ! at these height levels 
        extrapolation_loop: do ind_interp= (nr_h_lev_lowest2ECMWF + 1), nr_h_lev_lowest2upper_limit
            ! call the subroutine "standard_atm" for getting temperature and pressure at the desired height level
            ! note: see subroutine for optional arguments
            call standard_atm(h_lev_lowest2upper_limit(ind_interp), & ! input height level in [m]
                              T_extrap, & ! output T in [K]
                              p_extrap, & !  output p in [hPa]
                              .FALSE. ) ! unit type: .FALSE. for metric units
        
            ! assign extrapolated temperature and pressure to grid at extrpolation heigth level
            T_int(ind_interp, :, :)= T_extrap
            p_int(ind_interp, :, :)= p_extrap
        
            ! determine the water vapour pressure for the extrapolation height level
            ! --> humidity components of the air disappear in these heights above the ECMWF upper limit and
            ! the water vapour pressure is therefore set to 0
            wvpr_int(ind_interp, :, :)= 0 ! in [hPa]
            
            ! calculate the density of dry air for the extrapolation height level, [Md/R=1/Rd=Rd^-1]
            ! see scriptum Atmospheric Effects in Geodesy 2012, equation (2.53) on page 34
            !rho_d_int= ( p_int(ind_interp, :, :) - wvpr_int(ind_interp, :, :) ) / ( Rd * T_int(ind_interp, :, :) ) ! in [hPa*s^2/m^2]
            
            ! calculate auxiliary paramerters
            ! see Mendes & Pavlis 'High-accuracy zenith delay prediction at optical wavelengths'
            t_int_deg = T_extrap - 273.15d0 ! in [°C]
            xw = wvpr_int(ind_interp, :, :) / p_extrap
            Z_int = 1 - (p_extrap / T_extrap) * (a0 + a1*t_int_deg + a2*t_int_deg**2 + (b0 + b1*t_int_deg) * xw + (c0 + c1*t_int_deg) * xw**2) + (p_extrap / T_extrap)**2 * (d0 + e0*xw**2)
            
            ! calculate the density of water vapour for the interpolation height level
            ! see Mendes & Pavlis 'High-accuracy zenith delay prediction at optical wavelengths'
            rho_w_int= (Mw * xw * p_extrap) / (R * Z_int * T_extrap) ! in [hPa*s^2/m^2]
            
            ! calculate the total density of air for the interpolation height level
            ! see Mendes & Pavlis 'High-accuracy zenith delay prediction at optical wavelengths'
            rho_total_int= Md / (R * Z_int) * (p_extrap / T_extrap - (1 - Mw/Md) * wvpr_int(ind_interp, :, :) / T_extrap) ! in [hPa*s^2/m^2]
            
            ! calculate inverse compressibility factor of water vapour for the extrapolation height level
            ! see Kleijer 2004, Troposphere Modeling and Filtering for Precise GPS Leveling, equation  (4.31) on page 28
            !Zw_inv_int= 1 + wvpr_int(ind_interp, :, :) * ( 1 + 3.7d-4 * wvpr_int(ind_interp, :, :) ) * ( -2.37321d-3 + 2.23366d0 / T_int(ind_interp, :, :) - 710.792d0 / ( T_int(ind_interp, :, :) **2 ) + 7.75141d4 / ( T_int(ind_interp, :, :) **3 ) )
            
            ! calculate the refractive index for the extrapolation height level
            ! see Mendes & Pavlis 'High-accuracy zenith delay prediction at optical wavelengths'
            ! calculate the hydrostatic refractive index (first step: only assign non constant terms)
            n_h_int(ind_interp, :, :)= rho_total_int
            ! calculate the wet refractive index
            n_w_int(ind_interp, :, :)= Ngws * (rho_w_int/rho_ws) - Ngaxs * (Td_k/Pd) * (Mw/Md) * (Zd/Z_int) * (wvpr_int(ind_interp, :, :) / T_extrap)
        
        end do extrapolation_loop
    
    
        ! Finish calculation of refractive indices
    
        ! calculate the hydrostatic refractive index (second step: multiply with constants and convert from refractivity to refractive index)
        n_h_int= 1 + ( n_h_int * Ngaxs * (Td_k/Pd) * Zd * Rd ) * 1d-6
    
        ! calculate the wet refractive index (convert from refractivity to refractive index)
        n_w_int= 1 + n_w_int * 1d-6
    
        ! calculate total refractive index
        n_total_int= n_h_int + n_w_int - 1 ! -1 is necessary, because n_h_int and n_w_int are already converted to refractive indices and one +1 is two much if the two values are summed up!
    
    
        ! deallocate all unnecessary variables
        deallocate( gph_lower_lev, & ! variable for geopotential heights of lower level
                    gph_upper_lev, & ! variable for geopotential heights of upper level
                    h_orth_lower_lev, & ! variable for orthometric heights of lower level
                    h_orth_upper_lev, & ! variable for orthometric heights of upper level
                    h_ell_lower_lev, & ! variable for ellipsoidal heights of lower level
                    h_ell_upper_lev, & ! variable for ellipsoidal heights of upper level
                    wvpr_lower_lev, & ! variable for water vapour pressure of lower level
                    wvpr_upper_lev, & ! variable for water vapour pressure of upper level
                    T_lower_lev, & ! variable for temperature of lower level
                    T_upper_lev, & ! variable for temperature of upper level
                    Tv_lower_lev, & ! variable for virtual temperature of lower level
                    Tv_upper_lev, & ! variable for virtual temperature of upper level
                    C, & ! variable for auxiliary parameter C for logarithmic interpolation of the water vapour pressure for the interpolation height level
                    !rho_d_int, & ! variable for density of dry air for the interpolation height level
                    rho_w_int, & ! variable for density of water vapour for the interpolation height level
                    rho_total_int, & ! variable for total density of wet air for the interpolation height level
                    t_int_deg, & ! variable for temperature in degree celsius for the interpolation height level
                    xw, & ! variable for the ratio of water vapour pressure and total pressure for the interpolation height level
                    !Zw_inv_int ) ! variable for inverse compressibility factor of water vapour for the interpolation height level
                    Z_int ) ! variable for the compressibility factor of wet air for the interpolation height level

        
    end subroutine gridwise_refrHD_optical_ECMWFmin
    
end module module_gridwise_refrHD_optical_ECMWFmin