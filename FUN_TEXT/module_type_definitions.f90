! module_type_definitions.f90 

!****************************************************************************
!
! PURPOSE:  Module for defining types (structures) used in the RADIATE program
!
!----------------------------------------------------------------------------
! 
! Fortran-file created by Armin Hofmeister
! 
!----------------------------------------------------------------------------
! History:
! 
! 03.11.2014: create the Fortran-file
! 04.11.2014: enhance the Fortran-file
! 05.11.2014: programming
! 06.11.2014: comments
! 10.11.2014: programming
! 11.11.2014: programming
! 17.11.2014: programming
! 24.11.2014: programming
! 26.11.2014: programming, change the dimensional architecture of the content of the observations_type
! 27.11.2014: programming
! 01.12.2014: programming
! 02.12.2014: add grib_filename to epolog structure
! 03.12.2014: add enhanced date_type as epoch_type and add it to epolog structure
! 12.01.2015: add gribdata_type for storing grib-data from textfile
!             add profile_type and grid_type
! 13.01.2015: add meteo_type and add this to profile_type
! 20.01.2015: add meteo_int_type and add this to profile_type
! 29.01.2015: add meteo_stat_type and add this to profile_type
! 03.02.2015: add "nr_h_lev" to "meteo_int" structure
! 10.02.2015: add delay_type and add this to epolog_type
! 12.02.2015: programming
!             use constant to define length of station names
!             use constant to define length of source names
!             add raytrace_type and add this to epolog_type
! 16.02.2015: add "nr_h_lev" to "delay" structure
! 20.04.2015: correct comments
! 21.04.2015: add rdlog_type and rdlog_epoch_type for storing ray-tracing data
! 23.04.2015: add duplicate_obs_data_type
! 04.05.2015: comments
! 05.05.2015: add gridres_type and add it also to epolog_type
! 07.05.2015: add "cleanup" to parameter_type
! 12.05.2015: comments
! 02.06.2015: correct comments
! 26.08.2015: add comment
! 19.11.2015: correct comments
! 17.12.2015: changes due to azel file extension by water vapour pressure
!             add the total mapping factor to delay_type
!             add the meteo_stat_out_type for the output to the .radiate-file
!             add it to the epolog structure
! 21.12.2015: add comments

! Changes by Daniel Landskron:
! 19.11.2017: also enabled for individual internal creation of uniform azel-files by adding a 9th input argument
! 05.02.2018: epolog % raytrace is not globally stored anymore because it requires huge disk space and is not
!             needed outside of the get_RayTrace2D_* subroutines
!
!****************************************************************************
    
module module_type_definitions
    
    ! Define modules to be used
    use module_date_type_definition ! contains the date_type
    
    use module_constants, only: len_statname, len_souname, len_sessname, len_gribepochname
    

    ! type for storing options (parameters) for the ray-tracing calculations
    type, public :: parameters_type
        
        ! variables for
        
        ! epolog mode
        character(len=:), allocatable :: epolog_mode
        
        ! vertical interpolation method
        character(len=:), allocatable  :: interpolation_method
        
        ! ray-tracing method
        character(len=:), allocatable  :: raytr_method
        
        ! option of creating trp-files
        logical :: create_trp
        
        ! option of saving error log
        logical :: save_errorlog
        
        ! option for cleaning up of created session index- and epolog-files
        logical :: cleanup
        
        ! option for reading the azel file or creating it internally
        logical :: readAzel
        
        ! specifying path to file with station data
        character(len=:), allocatable :: load_path_stat_info
        
        ! specifying filename of file with station data
        character(len=:), allocatable :: load_filename_stat_info
        
        ! specifying path to AZEL-files
        character(len=:), allocatable :: load_path_azel
        
        ! specifying filename of input AZEL-file
        character(len=:), allocatable :: load_filename_azel
        
        ! specifying path to grib data
        character(len=:), allocatable :: load_path_grib
        
        ! specifying path to undulation data
        character(len=:), allocatable :: load_path_undulation
        
        ! setting the version info of the RADIATE program
        character(len=:), allocatable :: RADIATE_version
        
        ! setting the sub-version info of the RADIATE program
        character(len=:), allocatable :: RADIATE_subversion
        
        ! specifying prefix of session_index-files
        character(len=:), allocatable :: prefix_session_index
        
        ! specifying save path to session index file
        character(len=:), allocatable :: save_path_session_index
        
        ! specifying save path to epolog-files
        character(len=:), allocatable :: save_path_epologs
        
        ! specifying load path to session index file
        character(len=:), allocatable :: load_path_session_index
        
        ! specifying load path to epolog-files
        character(len=:), allocatable :: load_path_epologs
        
        ! specifying save path to error logs
        character(len=:), allocatable :: save_path_errorlog
        
        ! specifying save path to .radiate files
        character(len=:), allocatable :: save_path_radiate
        
        ! specifying save path to .radiate files
        character(len=:), allocatable :: save_path_trp
        
        ! specifying path to file with data for the individual azel creation
        character(len=:), allocatable :: load_path_filename_indAzelSpec
        
    end type parameters_type
    
    !----------------------------------------------------------------------------------------------
    
    ! type for storing azel-file data
    type, public :: observations_type
        
        ! variable for storing the observation data
        
        ! variable for storing the scannumber
        integer :: scannr
        
        ! variable for storing the modified julian date
        double precision :: mjd
        
        ! variable for storing the year
        integer :: year
        
        ! variable for storing the day of year
        integer :: doy
        
        ! variable for storing the hour
        integer :: hour
        
        ! variable for storing the minute
        integer :: min
        
        ! variable for storing the seconds
        double precision :: sec
        
        ! variable for storing the station name
        ! note: see "module_constants" for length of station names
        character(len= len_statname) :: station
        
        ! variable for storing the azimuth
        double precision :: az
        
        ! variable for storing the elevation
        double precision :: elev
        
        ! variable for storing the source name
        ! note: see "module_constants" for length of source names
        character(len= len_souname) :: source
        
        ! variable for storing the temperature
        double precision :: temp
        
        ! variable for storing the pressure
        double precision :: pres
        
        ! variable for storing the water vapour pressure
        double precision :: wvpr
                
    end type observations_type
    
    !----------------------------------------------------------------------------------------------
    
    ! type for storing the data of the observing stations
    type, public :: observing_stations_type
        
        ! variable for storing the station name
        ! note: see "module_constants" for length of station name
        character(len= len_statname) :: station
        
        ! variable for storing the number of observations per station
        integer :: nr_obs_per_obsstat
        
        ! variable for checking if ray-tracing should (needs to) be suspended for the observing station
        logical :: suspend_raytr = .FALSE. ! default initialization to .FALSE.
        
        ! variable for storing the ellipsoidal latitude coordinate
        double precision :: lat_ell
        
        ! variable for storing the ellipsoidal longitude coordinate
        double precision :: lon_ell
        
        ! variable for storing the ellipsoidal height coordinate
        double precision :: h_ell
           
    end type observing_stations_type
    
    !----------------------------------------------------------------------------------------------
    
    ! type for storing only the grid-resolution of the grib-file
    ! note: This type is needed although the resolution is also stored in the "grid_type" as the whole "profile"
    !       structure is deallocated after its use, but the resolution data is needed when creating the output files
    !       containing the results and errors.
    type, public :: gridres_type
        
        ! define variable for storing the grid resolution in latitude
        double precision :: res_lat
        
        ! define variable for storing the grid resolution in longitude
        double precision :: res_lon
        
    end type gridres_type
        
    !----------------------------------------------------------------------------------------------
    
    ! type for storing the epoch time information of a grib-file
    ! note: this is an extended structure of the date_type
    type, public, extends(date_type) :: epoch_type
        
        ! add mjd variable to date_type
        double precision :: mjd
        
    end type epoch_type
    
    !----------------------------------------------------------------------------------------------
    
    ! type for storing information about a grid
    type, public :: grid_type
        
        ! define variable for storing grid of latitude values
        ! note: dimension = nlat, nlon
        double precision, dimension(:, :), allocatable :: grid_lat
        
        ! define variable for storing grid of longitude values
        ! note: dimension = nlat, nlon
        double precision, dimension(:, :), allocatable :: grid_lon
        
        ! define variable for storing the grid size in latitude and longitude
        integer, dimension(2) :: grid_size
        
        ! define variable for storing the grid resolution in latitude
        double precision :: res_lat
        
        ! define variable for storing the grid resolution in longitude
        double precision :: res_lon
        
        ! define control variable for check if grid coverage in latitude is global
        logical ::  global_lat_check = .FALSE. ! default initialization to .FALSE.
        
        ! define control variable for check if grid coverage in longitude is global
        logical ::  global_lon_check = .FALSE. ! default initialization to .FALSE.
        
        ! define control variable for check if latitude of first grid point is equal to requirement
        logical ::  start_lat_check = .FALSE. ! default initialization to .FALSE.
        
        ! define control variable for check if longitude of first grid point is equal to requirement
        logical ::  start_lon_check = .FALSE. ! default initialization to .FALSE.
        
        ! define control variable for check if starting values of first grid point (latitude and longitude) apply to requirements
        ! as well as the grid coverage overall (latitude and longitude) is global
        logical ::  start_and_global_check = .FALSE. ! default initialization to .FALSE.
        
    end type grid_type
    
    !----------------------------------------------------------------------------------------------
    
    ! type for storing original grib-data and some derived parameters
    type, public :: meteo_type
        
        ! variable for storing the Z-records
        ! note dimension= (nr_pres_lev, nlat, nlon)
        double precision, dimension(:, :, :), allocatable :: Z
        
        ! variable for storing the Q-records
        ! note dimension= (nr_pres_lev, nlat, nlon)
        double precision, dimension(:, :, :), allocatable :: Q
        
        ! variable for storing the T-records
        ! note dimension= (nr_pres_lev, nlat, nlon)
        double precision, dimension(:, :, :), allocatable :: T
        
        ! variable for storing the pressure level values
        double precision, dimension(:), allocatable :: p
        
        ! variable for storing the geopotential height
        ! note dimension= (nr_pres_lev, nlat, nlon)
        double precision, dimension(:, :, :), allocatable :: gph
        
        ! variable for storing the water vapour pressure
        ! note dimension= (nr_pres_lev, nlat, nlon)
        double precision, dimension(:, :, :), allocatable :: wvpr
        
    end type meteo_type
    
    !----------------------------------------------------------------------------------------------
    
    ! type for storing vertically interpolated meteorological data and some additional data
    type, public :: meteo_int_type
        
        ! variable for storing the interpolation/extrapolation height levels
        double precision, dimension(:), allocatable :: h
        
        ! variable for storing the number of interpolation/extrapolation height levels
        integer :: nr_h_lev
    
        ! variable for storing the uppermost height level with data from the ECMWF grib-file
        double precision :: upper_limit_ECMWF
    
        ! variable for storing the upper limit of the atmosphere for the extrapolation of meteorologic values
        double precision :: upper_limit_atm
    
        ! variable for storing the interpolated/extrapolated gridded pressure values
        double precision, dimension(:, :, :), allocatable :: p
    
        ! variable for storing the interpolated/extrapolated gridded temperature values
        double precision, dimension(:, :, :), allocatable :: T
    
        ! variable for storing the interpolated/extrapolated gridded water vapour pressure values
        double precision, dimension(:, :, :), allocatable :: wvpr
    
        ! variable for storing the interpolated/extrapolated gridded hydrostatic refractive index
        double precision, dimension(:, :, :), allocatable :: n_h
    
        ! variable for storing the interpolated/extrapolated gridded wet refractive index
        double precision, dimension(:, :, :), allocatable :: n_w
    
        ! variable for storing the interpolated/extrapolated gridded total refractive index
        double precision, dimension(:, :, :), allocatable :: n_total
    
    end type meteo_int_type
    
    !----------------------------------------------------------------------------------------------
    
    ! type for storing interpolated meteorological data at the station position
    type, public :: meteo_stat_type
        
        ! variable for storing the indices of the grid point lat1lon1 next to the station position
        integer, dimension(2) :: ind_lat1lon1
        
        ! variable for storing the indices of the grid point lat1lon2 next to the station position
        integer, dimension(2) :: ind_lat1lon2
        
        ! variable for storing the indices of the grid point lat2lon2 next to the station position
        integer, dimension(2) :: ind_lat2lon2
        
        ! variable for storing the indices of the grid point lat2lon1 next to the station position
        integer, dimension(2) :: ind_lat2lon1
        
        ! variable for storing the interpolated pressure value at the station position
        double precision :: p
    
        ! variable for storing the interpolated temperature value at the station position
        double precision :: T
    
        ! variable for storing the interpolated water vapour pressure value at the station position
        double precision :: wvpr
    
        ! variable for storing the interpolated hydrostatic refractive index at the station position
        double precision :: n_h
    
        ! variable for storing the interpolated wet refractive index at the station position
        double precision :: n_w
    
        ! variable for storing the interpolated total refractive index at the station position
        double precision :: n_total
        
        ! variable for storing the index of the first height level in meteo_int above the station height (needed for
        ! ray-tracing start above station)
        integer:: start_lev
    
    end type meteo_stat_type
    
    !----------------------------------------------------------------------------------------------
    
    ! type for storing meteorological data interpolated at the station position for the output to the .radiate file
    ! Note: This type is used for storing temperature, pressure and water vapour pressure at the station position from the NWM
    !       in order to be output to the resulting .radiate file. According to the ray-tracing program options the values might
    !       be interpolated in time for the output.
    type, public :: meteo_stat_out_type
        
        ! variable for storing the interpolated pressure value at the station position
        double precision :: p
    
        ! variable for storing the interpolated temperature value at the station position
        double precision :: T
    
        ! variable for storing the interpolated water vapour pressure value at the station position
        double precision :: wvpr
    
    end type meteo_stat_out_type
    
    !----------------------------------------------------------------------------------------------
    
    ! type for storing grid info, the meteorologic profiles (from grib-file, refined (interpolated/extrapolated)) and the geoid undulations
    type, public :: profile_type
        
        ! substructure (declared above) for storing information about the grid
        type(grid_type) :: grid
        
        ! variable for storing the undulations fitting the grid
        double precision, dimension(:, :), allocatable :: undulation_grid
        
        ! substructure (declared above) for storing the meteorological data as reveived from the grib-file
        ! including some further parameters derived from the original values
        ! note: no vertical interpolation performed on the data
        type(meteo_type) :: meteo
        
        ! substructure (declared above) for storing vertically interpolated meteorological data and some additional data
        ! note: vertical interpolation performed on the data
        type(meteo_int_type) :: meteo_int
        
        ! type for storing interpolated meteorological data at the station position (dimension = number of observing stations)
        type(meteo_stat_type), dimension(:), allocatable :: meteo_stat
        
    end type profile_type
    
    !----------------------------------------------------------------------------------------------
    
    ! type for storing results from ray-tracing concerning the delays, the elevation angles, mapping factors and possible break in ray-tracing iteration
    type, public :: delay_type

        ! variable for storing the zenith total delay
        double precision :: dz_total
        
        ! variable for storing the zenith hydrostatic delay
        double precision :: dz_h
        
        ! variable for storing the zenith wet delay
        double precision :: dz_w
        
        ! variable for storing the slant total delay including geometric bending effect
        double precision :: ds_total_geom
        
        ! variable for storing the slant total delay
        double precision :: ds_total
        
        ! variable for storing the slant hydrostatic delay including geometric bending effect
        double precision :: ds_h_geom
        
        ! variable for storing the slant hydrostatic delay
        double precision :: ds_h
        
        ! variable for storing the slant wet delay
        double precision :: ds_w
        
        ! variable for storing the elevation angle at the station
        double precision :: e_stat
        
        ! variable for storing the (iteratively) ray-traced outgoing elevation angle
        double precision :: e_outgoing_rt
        
        ! variable for storing the geometric bending effect
        double precision :: dgeo
        
        ! variable for storing the value for total mapping factor (includes treatment of geometric bending effect)
        double precision :: mf_total_geom
        
        ! variable for storing the value for hydrostatic mapping factor (includes treatment of geometric bending effect)
        double precision :: mf_h_geom
        
        ! variable for storing the value for wet mapping factor
        double precision :: mf_w
        
        ! variable for storing the difference: outgoing (theoretical) - outgoing (ray-traced) elevation angle
        double precision :: diff_e
        
        ! variable for storing the logical signaling if a break in the while loop for calculating the outgoing elevation angle has occured (.FALSE. = no break, .TRUE. = break)
        logical :: break_elev = .FALSE. ! default initialization to .FALSE.
        
        ! variable for storing the logical signaling if a break in the while loop for calculating the next intersection point has occured (at least for one intersection point) ( .FALSE. = no break, .TRUE. = break)
        logical :: break_layer = .FALSE. ! default initialization to .FALSE.
        
    end type delay_type
    
    !----------------------------------------------------------------------------------------------
    
    ! type for storing results from ray-tracing concerning the ray path
    type, public :: raytrace_type
        
        ! note: most of the variables contained in the type are per height level, which means that their dimension
        !       will be allocated with [size(h_lev) = number of levels from station to end of ray path] in the subroutine
        
        ! variable for storing the earth radius calculated by gaussian curvature radius
        double precision :: R_e
        
        ! variable for storing the azimuth
        double precision :: az
        
        ! variable for storing the height levels from station height up to upper limit of the atmosphere
        double precision, dimension(:), allocatable :: h_lev
        
        ! variable for storing the number of height levels from station height up to upper limit of the atmosphere
        integer :: nr_h_lev
        
        ! variable for storing the position-dependend elevation angle in each height level
        double precision, dimension(:), allocatable :: theta
        
        ! variable for storing the elevation angle referenced to station position in each height level
        double precision, dimension(:), allocatable :: e
        
        ! variable for storing the geocentric angle between station and ray point in each height level
        double precision, dimension(:), allocatable :: anggeo
        
        ! variable for storing the latitudes of points along the ray path in each height level
        double precision, dimension(:), allocatable :: ray_lat
        
        ! variable for storing the longitudes of points along the ray path in each height level
        double precision, dimension(:), allocatable :: ray_lon
        
        ! variable for storing the index of point lat1lon1 in the global grid for bilinear interpolation of refractive index along ray-trace
        ! note: second dimension will be allocated with size=2 as for storing the two indices of lat and lon
        integer, dimension(:, :), allocatable :: ind_lat1lon1_trace
        
        ! variable for storing the index of point lat1lon2 in the global grid for bilinear interpolation of refractive index along ray-trace
        ! note: second dimension will be allocated with size=2 as for storing the two indices of lat and lon
        integer, dimension(:, :), allocatable :: ind_lat1lon2_trace
        
        ! variable for storing the index of point lat2lon2 in the global grid for bilinear interpolation of refractive index along ray-trace
        ! note: second dimension will be allocated with size=2 as for storing the two indices of lat and lon
        integer, dimension(:, :), allocatable :: ind_lat2lon2_trace
        
        ! variable for storing the index of point lat2lon1 in the global grid for bilinear interpolation of refractive index along ray-trace
        ! note: second dimension will be allocated with size=2 as for storing the two indices of lat and lon
        integer, dimension(:, :), allocatable :: ind_lat2lon1_trace
        
        ! variable for storing the slant total refractive index in each (intermediate for pwl) height level
        double precision, dimension(:), allocatable :: n_total
        
        ! variable for storing the slant hydrostatic refractive index in each (intermediate for pwl) height level
        double precision, dimension(:), allocatable :: n_h
        
        ! variable for storing the slant wet refractive index in each (intermediate for pwl) height level
        double precision, dimension(:), allocatable :: n_w
        
        ! variable for storing the zenith total refractive index in each (intermediate for pwl) height level
        double precision, dimension(:), allocatable :: n_total_z
        
        ! variable for storing the zenith hydrostatic refractive index in each (intermediate for pwl) height level
        double precision, dimension(:), allocatable :: n_h_z
        
        ! variable for storing the zenith wet refractive index in each (intermediate for pwl) height level
        double precision, dimension(:), allocatable :: n_w_z
        
    end type raytrace_type
    
    !----------------------------------------------------------------------------------------------
    
    ! type for storing epoch specific ray-tracing data and calculations (comprehension of substructures)
    type, public :: epolog_type
        
        ! variable for storing the session name
        ! note: see "module_constants" for length of session names
        character(len= len_sessname) :: session_name
        
        ! variable for storing the epoch name
        ! note: see "module_constants" for length of grib epoch name
        character(len= len_gribepochname) :: epoch_name
        
        ! structure (type declared above) for storing the epoch time data
        type(epoch_type) :: epoch
        
        ! variable for storing the grib-file name of the specific epolog (epoch)
        ! note: Grib-file filename without file extension. Should be equal to value of "epoch_name".
        character(len=:), allocatable :: grib_filename
        
        ! variable for storing the total number of observations in the epolog
        integer :: nr_obs
        
        ! variable for storing the number of observing stations (unique stations) in the epolog
        integer :: nr_observing_stations
        
        ! structure (type declared above) for storing the azel-file data
        ! note: each observation is one own entry in observations_type
        type(observations_type), dimension(:), allocatable :: observations
        
        ! structure (type declared above) for storing the data of the observing stations
        ! note: each observing station is one own entry in observing_stations_type
        type(observing_stations_type), dimension(:), allocatable :: observing_stations
        
        ! structure (type declared above) for storing grid info, the meteorologic profiles (from grib-file, refined (interpolated/extrapolated)) and the geoid undulations
        ! note: profile must be declared allocatable in order to be able to deallocate the structure later in the program when its usage is not needed any more in order to free RAM
        type(profile_type), allocatable :: profile
        
        ! structure (type declared above) for storing grid resolution
        ! note: this information is also part of the profile structure, but this is deallocated at a specific point and therefore the information needs to be stored again
        type(gridres_type) :: gridres
        
        ! structure (type declared above) for storing results from ray-tracing concerning the delays, the elevation angles, mapping factors and possible break in ray-tracing iteration
        type(delay_type), dimension(:), allocatable :: delay
        
        ! structure (type declared above) for storing meteorological data interpolated at the station position for the output to the .radiate file
        ! Note: The data might be interpolated in time for the output according to the settings.
        type(meteo_stat_out_type), dimension(:), allocatable :: meteo_stat_out
             
    end type epolog_type
    
    !----------------------------------------------------------------------------------------------
    
    ! type for storing ray-tracing data including the observation data
    type, public :: rdlog_type
        
        ! variable for storing the session name
        ! note: see "module_constants" for length of session names
        character(len= len_sessname) :: session_name
        
        ! structure (type declared above) for storing the azel-file data
        ! note: each observation is one own entry in observations_type
        type(observations_type), dimension(:), allocatable :: observations
        
        ! structure (type declared above) for storing results from ray-tracing concerning the delays, the elevation angles, mapping factors and possible break in ray-tracing iteration
        type(delay_type), dimension(:), allocatable :: delay
        
        ! structure (type declared above) for storing meteorological data interpolated at the station position for the output to the .radiate file
        ! Note: The data might be interpolated in time for the output according to the settings.
        type(meteo_stat_out_type), dimension(:), allocatable :: meteo_stat_out
        
    end type rdlog_type
    
    !----------------------------------------------------------------------------------------------
    
    ! type for storing the epoch time information of a grib-file
    ! note: this is an extended structure of the date_type
    type, public, extends(rdlog_type) :: rdlog_epoch_type
        
        ! add variable to rdlog_type for storing the mjd of the grib epoch at which each delay for each observation has been calculated
        double precision, dimension(:), allocatable :: grib_epoch_mjd
        
    end type rdlog_epoch_type
    
    !----------------------------------------------------------------------------------------------
    
    ! type for storing data from duplicate observations needed for time domain interpolation
    type, public :: duplicate_obs_data_type
        
        ! structure (type declared above) for storing results from ray-tracing concerning the delays, the elevation angles, mapping factors and possible break in ray-tracing iteration
        type(delay_type), dimension(:), allocatable :: delay
        
        ! structure (type declared above) for storing meteorological data interpolated at the station position for the output to the .radiate file
        ! Note: The data is interpolated in time for the output.
        type(meteo_stat_out_type), dimension(:), allocatable :: meteo_stat_out
        
        ! variable for storing the mjd of the grib epoch at which each delay for each observation has been calculated
        double precision, dimension(:), allocatable :: grib_epoch_mjd
        
        ! variable for storing the mjd of each observation
        double precision :: mjd_obs
        
    end type duplicate_obs_data_type
    
    !----------------------------------------------------------------------------------------------
    
    ! type for storing errors that occur during the ray tracing processing
    type, public :: errorlog_type
        
        ! variable for storing the error number
        integer :: error_nr
        
        ! variable for storing the error type definition
        character(len=:), allocatable :: error_type
        
        ! variable for storing the error description
        character(len=:), allocatable :: error_descr
        
    end type errorlog_type
    
    !----------------------------------------------------------------------------------------------
    
    ! type for storing grib-data (from grib-textfile or directly from grib-file)
    type, public :: gribdata_type
        
        ! variables for storing the number of grid points in latitude and longitude
        integer :: nlat, nlon
    
        ! variables for storing the latitude of the first and last grid point
        double precision :: lat_first, lat_last
    
        ! variables for storing the longitude of the first and last grid point
        double precision :: lon_first, lon_last
    
        ! variables for storing the grid resolution in latitude and longitude
        double precision :: dint_lat, dint_lon
    
        ! variable for storing the number of pressure levels
        integer :: nr_pres_lev
    
        ! variable for storing the pressure level values
        double precision, dimension(:), allocatable :: p
    
        ! variable for storing the number of parameters used from the grib-file (grib-textfile)
        integer :: nr_param
        
        ! variable for storing the Z-records
        ! note dimension= (nr_pres_lev, nlat * nlon)
        double precision, dimension(:, :), allocatable :: Z
        
        ! variable for storing the Q-records
        ! note dimension= (nr_pres_lev, nlat * nlon)
        double precision, dimension(:, :), allocatable :: Q
        
        ! variable for storing the T-records
        ! note dimension= (nr_pres_lev, nlat * nlon)
        double precision, dimension(:, :), allocatable :: T
    
        ! variables for storing the time of grib epoch
        type(epoch_type) :: epoch
        
    end type gribdata_type
    
    
end module module_type_definitions