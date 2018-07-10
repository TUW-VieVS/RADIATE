! module_gribtxt_get_meteo_undulation_refr_glob.f90

!****************************************************************************
!
! PURPOSE:  Module containing subroutine to read grib-file, derive meterological data, add undulation, interpolate meteorological values vertically.
!           
!           This function reads the ECMWF grib-file data (daily pressure-level data) from a textfile and delivers the
!           meteorological data of the 3 parameters Z, Q, T derived from the grib-file at one epoch. The delivered profile
!           contains the meteorologic values from the whole grib-file grid.
!           Undulation information is loaded according to the grid of the grib-file and used to calculate ellipsoidal heights.
!           The meteorological profiles will be interpolated vertically in order to increase height level resolution.
!
!           Attention: This version of the module "module_get_meteo_undulation_refr_glob" uses the grib-data stored in a textfile for loading!
!                      The grib data in th textfile may only contain one epoch!
!
!
! INPUT:
!       input_path_grib......... path to the grib-file
!       input_gribfilename...... filename of the grib-file
!       POIs_lat................ vector of ellips. latitudes of POIs in [°], interval [-90°,90°]
!       POIs_lon................ vector of ellips. longitudes of POIs in [°], interval [0°,360°[
!       POIs_h_ell.............. vector of ellips. heigths of POIs in [m]
!       POIs_name............... vector of names of POIs (station name)
!       input_path_undulation... path to the directory with the global geoid undulation files with
!                                different grid resolution
!       interpolation_method.... variable that defines the selected interpolation method
!       suspend_raytr........... variable which reports, if ray-tracing is not possible, because POI is not
!                                supported by the grid or has not been found in the catalogue
!                                .FALSE. ... ray-tracing will be done for this specific station
!                                .TRUE. ... ray-tracing will be suspended for this specific station
!       wavelength.............. wavelength range of observations (microwave or optical)
! 
! 
! OUTPUT:
!       suspend_raytr... variable which reports, if ray-tracing is not possible, because POI is not
!                        supported by the grid or has not been found in the catalogue
!       profile......... structure array that contains the grid info, the meteorologic profiles (from
!                        grib-file, refined (interpolated/extrapolated)) and the geoid undulations
!       epoch........... structure containing datum, time and mjd of grib-file epoch
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
! 02.12.2014: create the Fortran-file based on the Matlab-file "get_meteo_undulation_refr_ind_global.m"
! 03.12.2014: programming
! 08.01.2015: programming
! 12.01.2015: programming
! 13.01.2015: programming
! 14.01.2015: programming
! 15.01.2015: programming
! 19.01.2015: programming
! 20.01.2015: programming
! 21.01.2015: programming
! 22.01.2015: programming
! 26.01.2015: comments
! 28.01.2015: programming
! 29.01.2015: change to module as to save explicit interface block
!             use modules for subroutines with need of explicit interface (module_date2mjd, module_gridwise_refrHD_ECMWFmin,
!             module_profilewise_refrHD_ECMWFmin)
!             rename module and subroutine
!             programming
!             transfer constant "gn" to "module_constants"
! 03.02.2015: add output of "nr_h_lev" to subroutines for vertical interpolation
! 09.02.2015: programming
! 20.04.2015: correct comments and messages
! 05.05.2015: add report of grib-file resolution, declare dimension of input variables for POIs depending on "POIs_lat"
! 11.05.2015: add comments
! 25.08.2015: rename only the filename of this module to differentiate the version from the version of direct grib-file reading
!             add comments
! 31.08.2015: change == to .eqv. if logical values are compared as gfortran does not support == for logical comparison
!             set parenthesis in creation of vlat vector to avoid warning in gfortran
! 02.09.2015: comments
!             remove unecessary variables and the input variable "ind_epoch_grib"
! 03.09.2015: comments
! 06.10.2015: changes due to the input of the minimum station height to the subroutines for the vertical interpolation
! 21.12.2015: changes due to adding of the file extension of the grib-file in this subroutine instead of adding it in the main subroutine
! 11.01.2016: remove comments
! 14.01.2016: catch problem for vertical interpolation in case of no found station data for any observing station in the current
!             epolog,
!             avoid unnecessary calculations of ind_NPOI a further checks of valid station coordinates in case of stations with suspend_raytr == .TRUE.
!   
! Changes by Daniel Landskron:
! 05.02.2018: for Q, not only negative values are corrected but also those below 8*10^-7
!
! Changes by Janina Boisits:
! 27.02.2018: add input parameter 'wavelength',
!             add "module_gridwise_refrHD_optical_ECMWFmin" and "module_profilewise_refrHD_optical_ECMWFmin" to used modules,
!             add case selection when calling modules for computing refractivity (module_gridwise_refrHD_ECMWFmin/module_profilewise_refrHD_ECMWFmin or
!             module_gridwise_refrHD_optical_ECMWFmin/module_profilewise_refrHD_optical_ECMWFmin)
!
!****************************************************************************

module module_get_meteo_undulation_refr_glob

contains
    
    subroutine get_meteo_undulation_refr_glob( input_path_grib, &
                                               input_gribfilename, &
                                               POIs_lat, &
                                               POIs_lon, &
                                               POIs_h_ell, &
                                               POIs_name, &
                                               input_path_undulation, &
                                               interpolation_method, &
                                               suspend_raytr, &
                                               profile, &
                                               epoch, & 
                                               wavelength )
    
                                                 
        ! Define modules to be used
        use module_constants, only: gn
        use module_global_var
        use module_type_definitions, only: gribdata_type, profile_type, epoch_type
        use module_date_type_definition
        use module_meshgrid2D
        use module_date2mjd
        use module_gridwise_refrHD_ECMWFmin
        use module_profilewise_refrHD_ECMWFmin
        use module_gridwise_refrHD_optical_ECMWFmin
        use module_profilewise_refrHD_optical_ECMWFmin
        use module_calc_refr_ind_at_stations
    
        !----------------------------------------------------------------------------
    
        ! Variable definitions
        implicit none
    
    
        ! dummy arguments
        !----------------
    
        ! INPUT
    
        ! variable for input of path to the grib-file directory
        character(len=*), intent(in) :: input_path_grib
    
        ! variable for input of grib-filename
        ! note: grib-file filename without file extension
        character(len=*), intent(in) :: input_gribfilename
    
        ! vector for input vector of ellips. latitudes of POIs
        double precision, dimension(:), intent(in) :: POIs_lat
    
        ! vector for input vector of ellips. longitudes of POIs
        ! note: size is the number of input "POIs_lat"
        double precision, dimension(size(POIs_lat)), intent(in) :: POIs_lon
    
        ! vector for input vector of ellips. heights of POIs
        ! note: size is the number of input "POIs_lat"
        double precision, dimension(size(POIs_lat)), intent(in) :: POIs_h_ell
    
        ! vector for input vector of POIs names
        ! note: size is the number of input "POIs_lat"
        character(len=*), dimension(size(POIs_lat)), intent(in) :: POIs_name
    
        ! variable for input of path to the undulation-file directory
        character(len=*), intent(in) :: input_path_undulation
    
        ! variable for input of vertical interpolation method
        character(len=*), intent(in) :: interpolation_method
    
        ! variable for input of wavelength of observations
        character(len=*), intent(in) :: wavelength
    
    
        ! INPUT and OUTPUT
    
        ! define variable that reports if ray-tracing for a specific observing station needs to be suspended (default: .FALSE. --> no suspension == ray-tracing will be done)
        logical, dimension(:), intent(in out) :: suspend_raytr
    
    
        ! OUTPUT
    
        ! structure for storing grid info, the meteorologic profiles (from grib-file, refined (interpolated/extrapolated)) and the geoid undulations
        type(profile_type), intent(out) :: profile
    
        ! structure for storing datum, time and mjd of grib-file epoch
        type(epoch_type), intent(out) :: epoch
    
    
        ! local variables
        !----------------
    
        ! define variable for storing the number of pressure levels in the grib-file
        integer :: nr_pres_lev
    
        ! define variable for the grib-data text-filename
        ! note: grib-file filename with file extension
        character(len=:), allocatable :: input_gribfilename_txt
    
        ! define structure for storing the grib-data retrieved from the textfile
        type(gribdata_type) :: gribdata
    
        ! variables for storing the number of grid points in latitude and longitude
        integer :: nlat, nlon
    
        ! variables for storing the latitude of the first and last grid point
        double precision :: lat_first, lat_last
    
        ! variables for storing the longitude of the first and last grid point
        double precision :: lon_first, lon_last
    
        ! variables for storing the grid resolution in latitude and longitude
        double precision :: dint_lat, dint_lon
    
        ! define variables for storing the reference values of the desired starting position of the grid in latitude and longitude
        double precision :: reference_start_lat, reference_start_lon
    
        ! loop variable
        integer :: i
    
        ! variable for a vector storing the latitudes covered by the grid from the grib-file (lat_first:lat_last)
        double precision, dimension(:), allocatable :: vlat
    
        ! variable for a vector storing the longitudes covered by the grid from the grib-file (lon_first:lon_last)
        double precision, dimension(:), allocatable :: vlon
    
        ! define variable for index of pressure level
        integer :: ind_pres_lev
    
        ! variable for temporal storage of grib-data after reshape
        double precision, dimension (:, :), allocatable :: temp_2D_grid
    
        ! variable for temporal storage of grib-data after flipping pressure level dimension
        double precision, dimension (:, :, :), allocatable :: temp_3D_grid
    
        ! variable for temporal storage of pressure levels after flipping pressure levels
        double precision, dimension (:), allocatable :: temp_vec
    
        ! variable for storing the epoch civil date for calculating the epoch in mjd
        type(date_type) :: date_of_epoch
    
        ! variable for storing the undulation grid-size of a possible previously loaded undulation grid
        integer, dimension(2) :: und_size
    
        ! logical variable to determine if reloading of geoid undulations is necessary
        logical :: noreload
    
        ! variable for storing the number of points of interest, i.e. number of stations
        integer :: nr_POIs
    
        ! variables for indices of NPOIs, i.e. nearest grid point indices for latitude and longitude of the POIs
        integer, dimension(:), allocatable :: ind_NPOI_lat, ind_NPOI_lon
        
        ! variable for storing the minimum ellipsoidal height of all stations for which ray-tracing will be done in a specific epoch
        double precision :: min_stat_h
    
        
        ! CONSTANTS
        
        ! see "module_constants" for "gn": constant for normal gravity in [m/s^2]
        
        
        !----------------------------------------------------------------------------
    
        !============================================================================
        ! Body of the subroutine
        !============================================================================
    
        ! report process
        write(unit= *, fmt= '(3a)') 'Start reading the grib-file "', input_gribfilename,'" ...'
    
    
        ! Set global variable
        ! declared in module_global_var:
        !   undulation_grid_global
    
    
        ! Parameter definition
    
        !***************************************************************************************************
        !**************************** Define PARAMETER *****************************************************
    
        ! Note: Parameter table for grib-file is supposed to be be ECMWF128.
    
        ! define number of pressure levels in the grib-file
        ! !!! setting the correct value is important for getting the desired records!!!
        nr_pres_lev=25 ! Standard for ECMWF 128 grib-files: 25 pressure levels)
    
        ! Note: Only the 3 parameters Z,T,Q (Records nr. 129, 130, 133) will be read in by the subroutine
        !       "get_gribdata_txt".
    
        !***************************************************************************************************
        !***************************************************************************************************
    
    
        ! define the grib-data textfile name by adding the file extension
        input_gribfilename_txt= input_gribfilename // '.txt'
    
        ! call subroutine to read in the grib data that have been externally stored in a textfile
        call get_gribdata_txt( input_path_grib, &
                               input_gribfilename_txt, &
                               gribdata)
    
    
        ! report elapsed time
        call get_elapsed_time( icount_start, count_rate, icount_interm, elapsed_time )
        write(unit= *, fmt= '(a, f0.3, a, /)') 'Total elapsed time after reading the grib-file: ', elapsed_time, ' s'
        
        
        ! Get the information about the grid from the grib-file
        ! Note: In principle the following data could always be directly taken from the gridata-structure for later expressions,
        !       but this would lead to long not reader-friendly expressions.
    
        ! report process
        write(unit= *, fmt= '(a)') 'Start deriving the meteorological values ...'
    
        ! number of points along a latitude circle
        nlat= gribdata % nlat
    
        ! number of points along a (longitude) meridian
        nlon= gribdata % nlon
        
        ! latitude of first grid point
        lat_first= gribdata % lat_first ! possible values [90°,-90°]; first means max. latitude
    
        ! latitude of last grid point
        lat_last= gribdata % lat_last ! possible values [90°,-90°]; last means min. latitude
    
        ! longitude of first grid point
        lon_first= gribdata % lon_first ! possible values [0°,360°[; first means min. longitude
    
        ! longitude of last grid point
        lon_last= gribdata % lon_last ! possible values [0°,360°[; last means max. longitude
    
        ! grid intervall for latitude
        dint_lat = gribdata % dint_lat ! in [°]
    
        ! grid intervall for longitude
        dint_lon = gribdata % dint_lon ! in [°]
    
        ! report grib-file resolution
        ! note: fmt: ss --> no plus sign, f5.3 --> form x.yyy, 1 leading zero in case of shorter value
        write(unit= *, fmt= '(tr4, ss, a, f5.3, a, f5.3, a)') 'Grib-file resolution: res_lat= ', dint_lat, ' deg, res_lon= ', dint_lon, ' deg'
    
        ! check if the input grid is a global grid and starting with the desired values in latitude and
        ! longitude
        ! starting point of the grid must be 90° latitude and 0° longitude in order to apply undulations
        ! correctly in a later processing step (undulations start at 90° latitude and 0° longitude)
    
        ! set reference starting values for latitude and longitude used in the check
        reference_start_lat=90
        reference_start_lon=0
    
        ! call the subroutine to do the global-checks
        call test_start_values_and_global_coverage( nlat, &
                                                    nlon, &
                                                    dint_lat, &
                                                    dint_lon, &
                                                    lat_first, &
                                                    lon_first, &
                                                    reference_start_lat, &
                                                    reference_start_lon, &
                                                    profile % grid % global_lat_check, &
                                                    profile % grid % global_lon_check, &
                                                    profile % grid % start_lat_check, &
                                                    profile % grid % start_lon_check, &
                                                    profile % grid % start_and_global_check)
    
        ! error message in case of a grid that is not global or starting with a different latitude or
        ! longitude value
        if (.NOT. profile % grid % start_and_global_check) then
            ! report error in case that grib-file is not global or starting with wrong latitude or longitude
            write(unit= *, fmt= '(3a)') 'Error: Grib-file "', input_gribfilename, &
                                        '" does not deliver data for global coverage or grid boundaries are shifted to other values than 0° in longitude and 90° in latitude. Thus program is stopped (undulations can not be assigned correctly and other problems would occur)!'
            ! stop the program
            stop
        end if
    
    
        ! Create meshgrids which contain the latitude and longitude values of the global grid nodes
    
        ! allocate vectors for storing the latitude and longitude values
        allocate(vlat(nlat), vlon(nlon))
    
        ! create latitude vector using implied do loop
        ! note: loop index must start at 0 and end at nlat-1 to get vector from lat_first to lat_last
        vlat= [(lat_first + i*(-dint_lat), i= 0, nlat-1)] ! in [°], -dint_lat, because latitude is decreasing from first latitude entry per convention of grib-file, note: third value in implied loop would be the stride and is default set to 1
    
        ! create longitude vector using implied do loop
        ! note: loop index must start at 0 and end at nlon-1 to get vector from lon_first to lon_last
        vlon= [(lon_first + i*dint_lon, i= 0, nlon-1)] ! in [°], note: third value in implied loop would be the stride and is default set to 1
    
        ! allocate variables for storing the grids
        allocate(profile % grid % grid_lat(nlat, nlon), profile % grid % grid_lon(nlat, nlon))
    
        ! create meshgrids using the latitude and longitude vectors
        call meshgrid2D(vlon, &
                        vlat, &
                        profile % grid % grid_lon, &
                        profile % grid % grid_lat )
    
        ! determine size of global grid
        profile % grid % grid_size= [nlat, nlon]
    
        ! store grid resolution
        profile % grid % res_lat= dint_lat
        profile % grid % res_lon= dint_lon
    
        ! store global grid check values
        ! already done when calling the subroutine "test_start_values_and_global_coverage"
    
    
        ! Get values of recorded parameters (Z, Q, T) from the grib-file
    
        ! allocate arrays for the variables
        ! Z: geopotential in [m^2/s^2]
        ! Q: specific humidity in [kg/kg]
        ! T: temperature in [K]
        ! p: pressure in [hPa]
        allocate( profile % meteo % Z(nr_pres_lev,nlat,nlon), & 
                  profile % meteo % Q(nr_pres_lev,nlat,nlon), &
                  profile % meteo % T(nr_pres_lev,nlat,nlon), &
                  profile % meteo % p(nr_pres_lev) )
        
        ! Q-values:
        ! catch values below 8*10^-7 and set them to 8*10^-7
        where (gribdata % Q < 8.0d-007)
            gribdata % Q = 8.0d-007
        end where
    
        ! allocate temporal variable
        allocate(temp_2D_grid(nlon,nlat))
    
        ! loop over all pressure levels
        do ind_pres_lev= 1, nr_pres_lev
        
            ! note: reshape fills array columnwise
            ! note: temporal variable is necessary to split step of reshape and transpose to avoid stack overflow
        
            ! Z-values:
            ! to temporal variable: reshape the grib-data from from the textfile from a vector to a rectangular grid
            temp_2D_grid= reshape(gribdata % Z(ind_pres_lev, :), [nlon, nlat])
            ! transpose temporal variable and assign values to receive output in dimension(nlat, nlon)
            profile % meteo % Z(ind_pres_lev, :, :)= transpose(temp_2D_grid)
        
            ! Q-values:
            ! catch negative humidity values and set them to 0
            ! note: This step is already done prior to the loop over all pressure levels to speed up the processing.
        
            ! to temporal variable: reshape the grib-data from from the textfile from a vector to a rectangular grid
            temp_2D_grid= reshape(gribdata % Q(ind_pres_lev, :), [nlon, nlat])
            ! transpose temporal variable and assign values to receive output in dimension(nlat, nlon)
            profile % meteo % Q(ind_pres_lev, :, :)= transpose(temp_2D_grid)
        
            ! T-values:
            ! to temporal variable: reshape the grib-data from from the textfile from a vector to a rectangular grid
            temp_2D_grid= reshape(gribdata % T(ind_pres_lev, :), [nlon, nlat])
            ! transpose temporal variable and assign values to receive output in dimension(nlat, nlon)
            profile % meteo % T(ind_pres_lev, :, :)= transpose(temp_2D_grid)
        
        end do
    
        ! deallocate temporal variable
        deallocate(temp_2D_grid)
        
        
        ! Get pressure values at each pressure level from the grib-file
        profile % meteo % p = gribdata % p
    
    
        ! Flip meteorological values
    
        ! flip the entries in order to have the first pressure level (at 1000 hPa) and the according values as
        ! first entries in the matrices
        ! note: indexing: start:end:stride
    
        ! note: temporal variable is necessary to split to avoid stack overflow
    
        ! allocate temporal variables
        allocate(temp_3D_grid(nr_pres_lev, nlat, nlon))
        allocate(temp_vec(nr_pres_lev))
    
        ! flip Z
        temp_3D_grid = profile % meteo % Z(nr_pres_lev:1:-1, :, :)
        profile % meteo % Z = temp_3D_grid
    
        ! flip Q
        temp_3D_grid = profile % meteo % Q(nr_pres_lev:1:-1, :, :)
        profile % meteo % Q = temp_3D_grid
    
        ! flip T
        temp_3D_grid = profile % meteo % T(nr_pres_lev:1:-1, :, :)
        profile % meteo % T = temp_3D_grid
    
        ! flip p
        temp_vec = profile % meteo % p(nr_pres_lev:1:-1)
        profile % meteo % p = temp_vec
    
        ! deallocate temporal variables
        deallocate(temp_3D_grid, temp_vec)
    
    
        ! Determine time of epoch of the grib-file records for the selected epoch
    
        !epoch % year = gribdata % epoch % year
        !epoch % month = gribdata % epoch % month
        !epoch % day = gribdata % epoch % day
        !epoch % hour = gribdata % epoch % hour
        !epoch % min = gribdata % epoch % min
    
        epoch = gribdata % epoch
        epoch % sec = 0d0 ! set to 0 as usually not provided by the original grib-file
        
        ! calculate mjd for epoch
    
        ! assign civil date of epoch to date structure for calculating the mjd
        ! note: The structure "epoch" can not be used directly, because it is a different data type
        !       with more fields than required in the subroutine "date2mjd".
        date_of_epoch % year= epoch % year
        date_of_epoch % month= epoch % month
        date_of_epoch % day= epoch % day
        date_of_epoch % hour= epoch % hour
        date_of_epoch % min= epoch % min
        date_of_epoch % sec= epoch % sec
    
        ! call subroutine to calculate the mjd
        call date2mjd(date_of_epoch, epoch % mjd)
    
        !----------------------------------------------------------------
        
        
        ! clear grib variable to prevent overwriting and speed up calculation (prevent RAM overflow)
    
        ! deallocate gribdata-variables, which need a lot of memory
        deallocate(gribdata % p, gribdata % Z, gribdata % Q, gribdata % T)
    
        ! report elapsed time
        call get_elapsed_time( icount_start, count_rate, icount_interm, elapsed_time )
        write(unit= *, fmt= '(a, f0.3, a, /)') 'Total elapsed time after deriving the meteorological values: ', elapsed_time, ' s'
    
        !----------------------------------------------------------------
        
        
        ! Load global geoid undulation data
        ! This section loads the global geoid undulation data retrieved from the EGM2008. The resolution of
        ! the geoid undulation must be the same as the resolution of the ECMWF grib-file. This is recognized
        ! through the filename of the geoid undulations.
    
        ! report process
        write(unit= *, fmt= '(a)') 'Start importing the global geoid undulations ...'
        
        ! check if undulations have already been loaded in a previous run with the same size in lat and lon
        ! (--> same resolution)
    
        ! get actual size of global variable undulation_grid_global
        ! note: und_size is 0 if there hasn't been a prevous run
        ! initialize variable "und_size" with values 0 to be sure that size-determination of an not initialized "undulation_grid_global" works fine
        ! uninitialized variable-size should deliver size=0
        und_size = 0
    
        und_size(1) = size(undulation_grid_global, 1) ! size for dim= 1
        und_size(2) = size(undulation_grid_global, 2) ! size for dim= 2
    
        ! check size
        ! --> size in latitude and longitude should fit according to the current grib-file grid
        ! (then also the resolution fits)
        ! note: result of .AND.-expression is only .TRUE. in case all connected items are .TRUE.
        noreload = (und_size(1) == nlat .AND. und_size(2) == nlon) ! logical variable for defining necessity of reloading undulations (.TRUE. ... no reloading, .FALSE. ... reloading)
    
        ! check if reloading is necessary
        if (.NOT. noreload) then ! --> if noreload = .FALSE.
        
            ! allocate
            allocate(undulation_grid_global(nlat, nlon))
        
            ! call the subroutine "get_global_undulation" in order to load the global geoid undulations with the
            ! corresponding resolution compared to the grib-file
            call get_global_undulation( input_path_undulation, &
                                        dint_lat, &
                                        dint_lon, &
                                        nlat, &
                                        nlon, &
                                        undulation_grid_global )
         
        else
        
            ! report skipping of loading of the geoid undulations
            write(unit= *, fmt= '(tr4, a)') 'Skipping import of the global geoid undulations as needed undulations have already been loaded for a previous epolog!'
        
        end if

    
        ! store geoid undulation values from the global grid
        ! allocate the variable in the profile structure
        allocate(profile % undulation_grid(nlat, nlon))
        profile % undulation_grid = undulation_grid_global
    
    
        ! report elapsed time
        call get_elapsed_time( icount_start, count_rate, icount_interm, elapsed_time )
        write(unit= *, fmt= '(a, f0.3, a, /)') 'Total elapsed time after importing the global geoid undulations: ', elapsed_time, ' s'
        
        !----------------------------------------------------------------
        
        
        ! Main section to derive meteorologic profiles (global and at exact station position)
    
        ! report process
        write(unit= *, fmt= '(a)') 'Start deriving meteorological profiles ...'
    
        ! get number of POIs
        nr_POIs=size(POIs_lat, 1) ! dim= 1
    
        ! allocate variables for indices of NPOI
        allocate(ind_NPOI_lat(nr_POIs), ind_NPOI_lon(nr_POIs))
    
        ! intitialize variables for indices of NPOI
        ! note: for this variables the initialized value 0 is used to recognize "no data"!
        ind_NPOI_lat= 0
        ind_NPOI_lon= 0
    
        ! loop over all all desired POIs
        POI_loop: do i= 1, nr_POIs
            
            ! note: The check if NPOI has no valid data (NaN), i.e. POI is not found in catalogue and therefore has NaN values
            !       has already been done by checking for missing stations in "RayTrace_main_global".
            !       Nevertheless a check of "suspend_raytr" is introduced here since the following calculations of ind_NPOI
            !       with nint(NaN) result in Fortran in a number, but the calculationd can be avoided by the check anyway.
            if ( .NOT. suspend_raytr(i) ) then
                
                ! Calculate index of nearest point of interest to the actual POI and determine if POI lies in the grid area of the grib-file
                
                ! get index of the nearest point (NPOI) to the POI
                ! note: nint() returns nearest integer and returns an integer value, but input must be a kind of real
                ! usage of nint() determines which grid point is the nearest to the POI
                ! attention: for a positive latitude value of a POI that lies in the exact middle (.5) between
                !            two grid points, round delivers the higher (and integer!) index of the two possible, which
                !            results in an apparently lower latitude than expected when round would be used on the latitude
                !            value because of the data format begins with higher latitudes and goes to lower latitudes
                !            (for positive latitudes) delivering higher indices for lower (positive) latitudes
                ! note: nint() needs real input and delivers integer output
                ind_NPOI_lat(i)= nint( (lat_first - POIs_lat(i)) / dint_lat ) + 1 ! +1 because otherwise only the number of intervals is determined
                ind_NPOI_lon(i)= nint( (POIs_lon(i) - lon_first) / dint_lon ) + 1 ! +1 because otherwise only the number of intervals is determined
                
                ! determine if NPOI lies outside the grid
                ! catch also POIs in global grids between 360-dint_lon and 360
                ! note: result of .AND.-expression is only .TRUE. in case all connected items are .TRUE.
                ! note: Use .eqv. to check for == in case of logical values as gfortran does not support == for logicals.
                if ( (ind_NPOI_lat(i) < 1) .OR. (ind_NPOI_lat(i) > nlat) .OR. (ind_NPOI_lon(i) < 1) .OR. ( (ind_NPOI_lon(i) > nlon) .AND. (profile % grid % global_lon_check .eqv. .FALSE.) ) .OR. ( (ind_NPOI_lon(i) > nlon+1) .AND. (profile % grid % global_lon_check .eqv. .TRUE.) ) ) then
                    
                    ! report message
                    ! convert latitude and longitude of the affected station to a string of 8 characters in the form xxx.yyyy (sp --> force sign for all until sign specifier changed again, ss --> suppress plus sign, f8.4 --> leading zero in case of shorter value including possible sign)
                    write(unit= *, fmt= '(tr4, a, a, a, sp, f8.4, a, ss, f8.4, a)') 'Warning: POI does not lie in the grid supported by the selected grib-file: Station ', POIs_name(i), ', lat: ', POIs_lat(i),', lon: ', POIs_lon(i), '. Ray-tracing will be suspended for this station!'
                    ! set variable for suspending ray-tracing for the current POI to 1
                    suspend_raytr(i)= .TRUE.
                    
                end if
                
            end if
            
        end do POI_loop! end of loop over all POIs
    
    
        ! allocate arrays for the upcoming variables
        ! gph: geopotential in [m]
        ! wvpr: specific humidity in [hPa]
        allocate( profile % meteo % gph(nr_pres_lev,nlat,nlon), & 
                  profile % meteo % wvpr(nr_pres_lev,nlat,nlon) )
    
        ! calculate geopotential height in [m] from geopotential in [m^2/s^2] at all pressure levels
        ! note: since Fortran 90: two arrays that are conformable (same shape) are processed element-wise in a common expression
        !                         array and scalar in a statement lead to element-wise calculation of array elements
        profile % meteo % gph= profile % meteo % Z / gn ! in [m], see scriptum Atmospheric Effects in Geodesy, equation (1.62)
    
        ! calculate water vapour pressure
        ! see scriptum Atmospheric Effects in Geodesy, equation (1.11) for calculation of specific
        ! humidity and express the formula to e=...
        ! note: since Fortran 90: two arrays that are conformable (same shape) are processed element-wise in a common expression
        !                         array and scalar in a statement lead to element-wise calculation of array elements
        do i= 1, nr_pres_lev
            profile % meteo % wvpr(i, :, :)= profile % meteo % Q(i, :, :) * profile % meteo % p(i) / ( 0.622d0 + 0.378d0 * profile % meteo % Q(i, :, :) ) ! in [hPa]
        end do
    
        !----------------------------------------------------------------
        
        
        ! report process
        write(unit= *, fmt= '(tr4, a)') 'Start vertical interpolation ...'
        ! report elapsed time
        call get_elapsed_time( icount_start, count_rate, icount_interm, elapsed_time )
        write(unit= *, fmt= '(tr4, a, f0.3, a)') 'Total elapsed time at the start of the vertical interpolation: ', elapsed_time, ' s'
        
        
        ! determine the minimum station height used for the start of the vertical interpolation
        
        ! check if ray-tracing is suspended for all observing stations in the current epolog
        ! note: This test is necessary since all stations have in this case their height set to NaN, which would lead
        !       to problems in case minval(POIs_h_ell) is used for the minimum station height as it would deliver an infinity value.
        !       As a remedy a dummy value for the minimum station height is set in order to be able to proceed with the program
        !       since later epologs of the session may still work out fine if their station data have been found and a program stop
        !       would lead to possible loss of results from these epologs.
        ! note: all() determines if all elements of an array have the value .TRUE.
        if ( all(suspend_raytr) ) then
            ! set a dummy value for the minimum station height
            ! note: Usage of the uppermost possible height for the vertical interpolation, i.e. 84000 m is sufficient to reduce the processing time.
            min_stat_h= 84000d0 ! in [m]
            
        else
            ! note: minval(POIs_h_ell) returns the minimum station height.
            min_stat_h= minval(POIs_h_ell)
            
        end if
        
        
        ! calculate the refractive index profiles (grids for each desired height level)
        ! determine the method which is used for interpolation
        select case (interpolation_method)
            ! case of profilewise interpolation approach
            case ('profilewise')
                ! determine wavelength range of observations
                select case (wavelength)
                    ! case of microwave observations
                    case ('microwave')
                        ! call subroutine for profilewise interpolation (microwave)
                        call profilewise_refrHD_ECMWFmin( min_stat_h, &
                                                          profile % grid % grid_lat, &
                                                          profile % grid % grid_size, &
                                                          nr_pres_lev, &
                                                          profile % meteo % gph, &
                                                          profile % meteo % p, &
                                                          profile % meteo % wvpr, &
                                                          profile % meteo % T, &
                                                          profile % undulation_grid, &
                                                          profile % meteo_int % h, &
                                                          profile % meteo_int % nr_h_lev, &
                                                          profile % meteo_int % upper_limit_ECMWF, &
                                                          profile % meteo_int % upper_limit_atm, &
                                                          profile % meteo_int % p, &
                                                          profile % meteo_int % T, &
                                                          profile % meteo_int % wvpr, &
                                                          profile % meteo_int % n_h, &
                                                          profile % meteo_int % n_w, &
                                                          profile % meteo_int % n_total )
                    ! case of optical observations
                    case ('optical')
                        ! call subroutine for profilewise interpolation (optical)
                        call profilewise_refrHD_optical_ECMWFmin( min_stat_h, &
                                                                  profile % grid % grid_lat, &
                                                                  profile % grid % grid_size, &
                                                                  nr_pres_lev, &
                                                                  profile % meteo % gph, &
                                                                  profile % meteo % p, &
                                                                  profile % meteo % wvpr, &
                                                                  profile % meteo % T, &
                                                                  profile % undulation_grid, &
                                                                  profile % meteo_int % h, &
                                                                  profile % meteo_int % nr_h_lev, &
                                                                  profile % meteo_int % upper_limit_ECMWF, &
                                                                  profile % meteo_int % upper_limit_atm, &
                                                                  profile % meteo_int % p, &
                                                                  profile % meteo_int % T, &
                                                                  profile % meteo_int % wvpr, &
                                                                  profile % meteo_int % n_h, &
                                                                  profile % meteo_int % n_w, &
                                                                  profile % meteo_int % n_total )
                end select

    
            ! case of gridwise interpolation approach
            case ('gridwise')
                ! determine wavelength range of observations
                select case (wavelength)
                    ! case of microwave observations
                    case ('microwave')
                        ! call subroutine for gridwise interpolation (microwave)
                        call gridwise_refrHD_ECMWFmin( min_stat_h, &
                                                       profile % grid % grid_lat, &
                                                       profile % grid % grid_size, &
                                                       nr_pres_lev, &
                                                       profile % meteo % gph, &
                                                       profile % meteo % p, &
                                                       profile % meteo % wvpr, &
                                                       profile % meteo % T, &
                                                       profile % undulation_grid, &
                                                       profile % meteo_int % h, &
                                                       profile % meteo_int % nr_h_lev, &
                                                       profile % meteo_int % upper_limit_ECMWF, &
                                                       profile % meteo_int % upper_limit_atm, &
                                                       profile % meteo_int % p, &
                                                       profile % meteo_int % T, &
                                                       profile % meteo_int % wvpr, &
                                                       profile % meteo_int % n_h, &
                                                       profile % meteo_int % n_w, &
                                                       profile % meteo_int % n_total )
                    ! case of optical observations
                    case ('optical')
                        ! call subroutine for gridwise interpolation (optical)
                        call gridwise_refrHD_optical_ECMWFmin( min_stat_h, &
                                                               profile % grid % grid_lat, &
                                                               profile % grid % grid_size, &
                                                               nr_pres_lev, &
                                                               profile % meteo % gph, &
                                                               profile % meteo % p, &
                                                               profile % meteo % wvpr, &
                                                               profile % meteo % T, &
                                                               profile % undulation_grid, &
                                                               profile % meteo_int % h, &
                                                               profile % meteo_int % nr_h_lev, &
                                                               profile % meteo_int % upper_limit_ECMWF, &
                                                               profile % meteo_int % upper_limit_atm, &
                                                               profile % meteo_int % p, &
                                                               profile % meteo_int % T, &
                                                               profile % meteo_int % wvpr, &
                                                               profile % meteo_int % n_h, &
                                                               profile % meteo_int % n_w, &
                                                               profile % meteo_int % n_total )
                end select
            
            ! case if other cases fail
            case default
                ! report error message
                print '(tr8, a)', 'Error: Unexpected interpolation method! Please set correct interpolation method! Program stopped!'
                ! stop the program
                stop
            
            end select
    
        ! clean up (deallocate) the variables in the "meteo" structure of the structure "profile", which contains grib-file values that are not needed any further for calculations
        deallocate( profile % meteo % Z, &
                    profile % meteo % Q, &
                    profile % meteo % T, &
                    profile % meteo % p, &
                    profile % meteo % gph, &
                    profile % meteo % wvpr ) 
    
        ! report elapsed time
        call get_elapsed_time( icount_start, count_rate, icount_interm, elapsed_time )
        write(unit= *, fmt= '(tr4, a, f0.3, a, /)') 'Total elapsed time after the vertical interpolation: ', elapsed_time, ' s'
        
        !----------------------------------------------------------------
        
        
        ! Get interpolated meteorological values at station position and height
        ! This section is necessary to interpolate the meteorological values at the station
        ! positions and their heights using the previously interpolated height levels.
        ! For each station the meteorlogical values will be determined and the refractive index will be
        ! calculated.
    
        ! report process
        write(unit= *, fmt= '(tr4, a)') 'Start deriving meteorological values at the stations ...'
        
        ! allocate sub-structure in structure "profile" for storing the meteorological values at the station positions
        ! that will be determined in "calc_refr_ind_at_stations"
        ! note: size = number of observing stations
        allocate( profile % meteo_stat(nr_POIs))
        
        ! call subroutine to determine the refractive index at the station position and height
        call calc_refr_ind_at_stations( POIs_lat, &
                                        POIs_lon, &
                                        POIs_h_ell, &
                                        profile % meteo_int, &
                                        dint_lat, &
                                        dint_lon, &
                                        profile % grid % grid_lat, &
                                        profile % grid % grid_lon, &
                                        profile % grid % grid_size, &
                                        profile % grid % start_and_global_check, &
                                        suspend_raytr, &
                                        profile % meteo_stat )
        
        ! report elapsed time
        call get_elapsed_time( icount_start, count_rate, icount_interm, elapsed_time )
        write(unit= *, fmt= '(tr4, a, f0.3, a)') 'Total elapsed time after deriving meteorological values at the stations: ', elapsed_time, ' s'
        
        ! report elapsed time
        call get_elapsed_time( icount_start, count_rate, icount_interm, elapsed_time )
        write(unit= *, fmt= '(a, f0.3, a, /)') 'Total elapsed time after deriving meteorological profiles: ', elapsed_time, ' s'
        
        !----------------------------------------------------------------
        
        
        ! report process
        write(unit= *, fmt= '(a)') 'Start cleaning up unneeded values ...'
        
        ! Free memory
        ! deallocate variables in "profile % meteo_int" that are not needed any further for calculations
        deallocate( profile % meteo_int % p, &
                    profile % meteo_int % T, &
                    profile % meteo_int % wvpr )
        
        ! deallocate undulation values stored in "profile % undulation_grid"
        ! note: only deallocate undulations stored in profile structure not the global variable "undulation_grid_global" as it may be needed for next epolog
        deallocate( profile % undulation_grid )
        
        ! report elapsed time
        call get_elapsed_time( icount_start, count_rate, icount_interm, elapsed_time )
        write(unit= *, fmt= '(a, f0.3, a, /)') 'Total elapsed time after cleaning up unneeded values: ', elapsed_time, ' s'
        
     
    end subroutine get_meteo_undulation_refr_glob
    
end module module_get_meteo_undulation_refr_glob