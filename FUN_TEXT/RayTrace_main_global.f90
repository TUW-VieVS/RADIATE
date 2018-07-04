! RayTrace_main_global.f90 

!****************************************************************************
!
! PURPOSE:  Main subroutine in RADIATE program for processing tasks
!            
!           This subroutine is the main program subroutine for ray-tracing.
!           It contains the main steps to call several subroutines in order to load, process and save data used
!           to determine the final outputs like total path delay and zenith delay for observations of
!           different VLBI stations in specific sessions.
! 
! 
! INPUT:
!         parameters... structure containing processing modes and options, paths and filenames
!                       needed to carry out the ray-tracing program
! 
! 
! OUTPUT:
!         parameters... some parameters have been added during the processing
!
!----------------------------------------------------------------------------
! 
! Fortran-file created by Armin Hofmeister
!   based on Matlab scripts created by Armin Hofmeister
! 
!----------------------------------------------------------------------------
! History:
! 
! 03.11.2014: create the Fortran-file based on the Matlab-file "RayTrace_main_global.m"
! 04.11.2014: proceed with Fortran translation of Matlab code
! 05.11.2014: proceed with Fortran translation of Matlab code
! 06.11.2014: programming
! 13.11.2014: comments
! 18.11.2014: programming
! 20.11.2014: programming
! 24.11.2014: programming
! 25.11.2014: programming
! 26.11.2014: programming
! 27.11.2014: programming
! 01.12.2014: programming
! 02.12.2014: programming
! 03.12.2014: programming
! 08.01.2015: programming
! 12.01.2015: programming
! 20.01.2015: formatting
! 21.01.2015: comments
! 28.01.2015: comments
! 29.01.2015: use modules for subroutines with need of explicit interface (module_import_station_data, module_load_session_index,
!             module_get_coord_observing_stations, module_get_coord_observing_stations, module_get_meteo_undulation_refr_glob)
! 09.02.2015: programming
! 10.02.2015: programming
! 11.02.2015: programming
! 12.02.2015: programming
!             use constant to define length of station names
!             use constant to define length of session names
! 16.02.2015: programming
! 18.02.2015: programming
! 25.02.2015: programming
! 26.02.2015: programming
! 02.03.2015: programming
! 09.04.2015: programming
! 20.04.2015: programming
!             correct comments and messages
! 21.04.2015: programming
! 22.04.2015: programming
! 23.04.2015: programming
! 27.04.2015: programming
! 28.04.2015: programming
! 30.04.2015: programming
! 04.05.2015: programming
! 05.05.2015: programming
! 06.05.2015: programming
! 07.05.2015: programming
! 11.05.2015: programming
! 12.05.2015: programming
! 13.05.2015: programming
! 02.06.2015: programming
! 03.06.2015: programming
! 09.06.2015: comments
! 11.06.2015: comments
! 31.08.2015: change \ to / in path definitions for Windows and Linux compatibility
!             change == to .eqv. if logical values are compared as gfortran does not support == for logical comparison
! 02.09.2015: remove variable "ind_epoch_grib" for input to subroutine "get_meteo_undulation_refr_glob"
! 10.09.2015: change directory names to upper case only
!             change usage of "epolog" to "epolog1" and "session_index" to "session_index1"
! 14.09.2015: changes of input variables for subroutines "get_observing_stations" and "get_coord_observing_stations"
! 01.10.2015: correct prefix for session index file for mode "one_epoch_per_obs" to "session_index1_"
! 19.11.2015: correct comments
! 17.12.2015: changes due to the output of the total mapping factor and the meteorological data at the station position from the NWM to the .radiate-file
! 21.12.2015: move adding of grib-file file extension to the subroutine where the file is really loaded
! 11.01.2016: correct indentings at do loop for assigning values to substructure "meteo_stat_out"
! 13.01.2016: correct comments
! 14.01.2016: correct message
    
! Changes by Daniel Landskron:
! 19.11.2017: also enabled for individual internal creation of uniform azel-files
! 05.02.2018: epolog % raytrace is not globally stored anymore because it requires huge disk space and is not
!             needed outside of the get_RayTrace2D_* subroutines
!
!****************************************************************************
    
    
subroutine RayTrace_main_global( parameters )
    
    ! Define modules to be used
    
    ! Module for type definitions
    use, intrinsic :: ieee_arithmetic ! intrinsic module to e.g. specify NaN and Inf
    use module_type_definitions, only: parameters_type, epolog_type, errorlog_type, rdlog_epoch_type, observing_stations_type
    use module_constants, only: len_sessname, len_statname, len_epolog_filename
    use module_import_station_data
    use module_load_session_index
    use module_get_observing_stations
    use module_get_coord_observing_stations
    use module_resize_errorlog
    use module_get_meteo_undulation_refr_glob
    use module_global_var, only: icount_start, count_rate, icount_interm, elapsed_time ! not needed here: "undulation_grid_global"
    use module_get_bilint_value
    use module_get_RayTrace2D_pwl_global
    use module_get_RayTrace2D_ref_pwl_global
    use module_get_RayTrace2D_Thayer_global
    use module_combine_and_sort_rd
    use module_time_interpolation
    use module_create_radiate_global
    use module_get_unique_stations_with_coord
    use module_create_trp_global
    use module_create_errorlog
    
    !----------------------------------------------------------------------------
    
    ! Variable definitions
    implicit none
    
    ! dummy argument
    !---------------
    
    ! define structure for parameter settings for the calculation process
    type(parameters_type), intent(in out) :: parameters
    
    
    ! local variables
    !----------------
    
    ! define variable for storing the "NaN" value
    double precision :: NaN
    
    ! define variable for storing the current date
    character(len=8) :: date_value
    character(len=10) :: date_formatted
    
    ! define variable for storing the current time
    character(len=10) :: time_value
    character(len=12) :: time_formatted
    
    ! define variable for storing the session name
    ! note: see "module_constants" for length of session names
    character(len= len_sessname) :: session_name
    
    ! define variable for storing the station names
    ! note: see "module_constants" for length of station names
    character(len= len_statname), dimension(:), allocatable :: name_list ! attention: dimension will be specified through counting the number of lines in the file
    
    ! define variable for storing the ellipsoidal latitude values in [°]
    double precision, dimension(:), allocatable :: lat_ell ! attention: dimension will be specified through counting the number of lines in the file
    
    ! define variable for storing the ellipsoidal longitude values in [°]
    double precision, dimension(:), allocatable :: lon_ell ! attention: dimension will be specified through counting the number of lines in the file
    
    ! define variable for storing the ellipsoidal height values in [m]
    double precision, dimension(:), allocatable :: h_ell ! attention: dimension will be specified through counting the number of lines in the file
    
    ! define character string array for storing the filenames of the epologs
    ! note: see "module_constants" for "len_epolog_filename"
    character(len=len_epolog_filename), dimension(:), allocatable :: epolog_filenames
    
    ! define variable for storing the number of epolog-files that make up the session
    integer :: nr_epologs
    
    ! define loop variable for epolog index
    integer :: curr_epolog_ind
    
    ! define structure for storing the epolog data and other data necessary for ray-tracing
    ! note: allocatable as number of epologs for the input session will be determined during run-time
    type(epolog_type), dimension(:), allocatable :: epolog
    
    ! define variable for the unique station names
    ! note: see "module_constants" for length of station names
    character(len= len_statname), dimension(:), allocatable :: obs_stations
    
    ! define variable for the number of observations per (unique) station
    integer, dimension(:), allocatable :: nr_obs_per_obsstat
    
    ! define structure for errors of missing station ccordinates in one epolog
    ! note: allocatable as number of errors will be determined during run-time
    type(errorlog_type), dimension(:), allocatable :: errorlog_curr_epolog_obsstat_coord
    
    ! define structure for gathering all errors
    ! note: allocatable as number of errors will be determined during run-time
    type(errorlog_type), dimension(:), allocatable :: errorlog
    
    ! define variable for storing the current number of occured errors (just a part of the total errors)
    integer :: nr_curr_errors
    
    ! define variable for storing the total number of occured errors (prior to new errors)
    integer :: nr_total_errors
    
    ! define variable for storing the total number of occured errors (with added new errors)
    integer :: new_nr_total_errors
    
    ! define variable for storing the start index of the observations for a specific station
    integer :: obs_begin
    
    ! define variable for storing the end index of the observations for a specific station
    integer :: obs_end
    
    ! define loop variable for current station
    integer :: ind_stat
    
    ! define variable to store values of refractive indices of the height level directly above the station level at the station's horizontal position
    ! note: these variables are only needed in case of using "pwl" for ray-tracing
    double precision :: n_h_stat_upper
    double precision :: n_w_stat_upper
    double precision :: n_total_stat_upper
    
    ! define variable for storing the total number of height levels
    integer :: nr_lev
    
    ! define loop variable
    integer :: i
    
    ! define variable for the name of the currently processed station
    ! note: see "module_constants" for length of station names
    character(len= len_statname) :: curr_stat
    
    ! define loop variable
    integer :: r
    
    ! define variable for index of current observation of the current station
    integer :: obs_ind_curr_stat
    
    ! define variable for storing the index of the current error
    integer :: ind_error
    
    ! define variable for storing error messages
    ! note: variables are initialized with specific length that should be long enough to hold the description
    !       This is necessary as the write statement needs an allocated variable.
    !       Assigning a new string to the variable that is shorter than the old one is no problem as the remaining
    !       space in the string is automatically filled with spaces
    character(len=1000) :: error_description
    
    ! define variable for storing the combination of all observations and their delays (including duplicate observations with delays calculated at other epochs)
    type(rdlog_epoch_type) :: rdlog_epoch
    
    ! define structure for storing the combined observing station data of each epolog (possibly including duplicates)
    type(observing_stations_type), dimension(:), allocatable :: observing_stations_cumul_epolog
    
    ! define variable for storing the index of the beginning position where to add epolog data to the combined structure "observing_stations_cumul_epolog"
    integer :: ind_begin
        
    ! define variable for storing the index of the end position where to add epolog data to the combined structure "observing_stations_cumul_epolog"
    integer :: ind_end
    
    ! define structure for storing the session-unique observing station data
    type(observing_stations_type), dimension(:), allocatable :: observing_stations_session
    
        
    ! CONSTANTS
    
    ! see "module_constants" for "len_sessname" "len_statname" and "len_epolog_filename"
    
    
    !----------------------------------------------------------------------------
    
    !============================================================================
    ! Body of the subroutine
    !============================================================================

    ! Print informations about the program start and the chosen parameter options
    
    ! get the current date and time
    call DATE_AND_TIME(date= date_value, time= time_value)
    
    ! start timing of the program
    call system_clock(icount_start, count_rate)
    
    ! format the date and time
    date_formatted= date_value(1:4) // '-' // date_value(5:6) // '-' // date_value(7:8)
    time_formatted= time_value(1:2) // ':' // time_value(3:4) // ':' // time_value(5:)
    
    ! display start of the program
    write(unit= *, fmt= '(a)') 'RADIATE Ray-tracing program'
    write(unit= *, fmt= '(tr4, a, /)') 'created by Armin Hofmeister'
    
    ! write current date and time
    write(unit= *, fmt= '(a)') 'Start time:'
    write(unit= *, fmt= '(2( tr3, 2a, /))') 'Date: ', date_formatted, 'Time: ', time_formatted
    
    
    !---------------------------------------------------------------------------------------------------
    
    ! Define program version of the RADIATE program
    
    ! define the version info
    parameters % RADIATE_version= '2.0_Fortran'
    
    ! define subversion info
    parameters % RADIATE_subversion= 'global_limit'
    
    
    !---------------------------------------------------------------------------------------------------
    
    ! display program version info
    write(unit= *, fmt= '(a, a)') 'Version: ', parameters % RADIATE_version
    write(unit= *, fmt= '(a, a, /, /)') 'Sub-version: ', parameters % RADIATE_subversion
    
    
    !---------------------------------------------------------------------------------------------------
    
    ! get session name
    ! find first "_" after "azel" and take filename until exclusive "."
    session_name=parameters % load_filename_azel(index(parameters % load_filename_azel, '_')+1:index(parameters % load_filename_azel, '.')-1)
    
    
    !---------------------------------------------------------------------------------------------------
    
    ! display session name and starting time of processing
    write(unit= *, fmt= '(a)') '**************************************'
    write(unit= *, fmt= '(3a)') '* Processed session: ', session_name, '  *'
    write(unit= *, fmt= '(a, /)') '**************************************'
    
    ! write header
    write(unit= *, fmt= '(a)') '==================='
    write(unit= *, fmt= '(a)') '= Parametrization ='
    write(unit= *, fmt= '(a, /)') '==================='
    
    ! write info about option how station data is imported
    write(unit= *, fmt= '(a, a, a, /)') 'Station data are loaded from textfile: "', parameters % load_path_stat_info // parameters % load_filename_stat_info, '"'
    
    ! write info about chosen mode of assigning the observations to epochs
    write(unit= *, fmt= '(a, a, /)') 'Mode of assigning the observations to epochs: ', parameters % epolog_mode
    
    ! write info about chosen extension in lat and lon of the subgrid around a station
    write(unit= *, fmt= '(a, /)') 'Global meteorological grid is used for ray-tracing at all stations!'
    
    ! write info about chosen mode of interpolation
    write(unit= *, fmt= '(a, a, /)') 'Interpolation mode: ', parameters % interpolation_method
    
    ! write info about chosen ray-tracing method
    write(unit= *, fmt= '(a, a, /)') 'Ray-tracing approach chosen: ', parameters % raytr_method
    
    ! report if .trp-file will be created
    if (parameters % create_trp) then
        ! report that .trp-file will be created
        write(unit= *, fmt= '(a, /)') '.trp-file will be created.'
    else
        ! report that no .trp-file will be created
        write(unit= *, fmt= '(a, /)') 'No .trp-file will be created.'
    end if
    
    ! write info about the chosen option for saving an error log
    if (parameters % save_errorlog) then
        write(unit= *, fmt= '(a, /)') 'Option for saving an error log is active. An .err-file containing errors encountered during the processing will be saved.'
    else
        write(unit= *, fmt= '(a, /)') 'Option for saving an error log is inactive. No .err-file containing errors encountered during the processing will be saved.'
    end if
    
    ! write info about the chosen option for cleaning up created session index- and epolog-files
    if (parameters % cleanup) then
        write(unit= *, fmt= '(a, /)') 'Option for cleaning up created session index- and epolog-files is active. Session index- and epolog-files will be deleted after their use.'
    else
        write(unit= *, fmt= '(a, /)') 'Option for cleaning up created session index- and epolog-files is inactive. Session index- and epolog-files will be kept.'
    end if
    
    ! write info about the chosen option for reading or internally creating the uniform azel file
    if (parameters % readAzel) then
        write(unit= *, fmt= '(a, a, a /)') 'The azel-file ', parameters % load_path_azel // parameters % load_filename_azel , ' is read.'
    else
        write(unit= *, fmt= '(a, a, a /)') 'A uniform azel-file is created internally for the individual settings from ', parameters % load_path_filename_indAzelSpec , '. No azel file is read.'
    end if
    
    
    !---------------------------------------------------------------------------------------------------
    
    ! Define paths and filenames if necessary for saving and loading of variables and results
    !
    ! epologs and session indices
    !
    ! define prefix for session index name (used for defining the loading filename) depending on the
    ! chosen epolog creation mode
    ! define path for saving the new session index and epolog textfiles depending on the chosen epolog
    ! creation mode
    select case (parameters % epolog_mode)
    
        case ('one_epoch_per_obs')
        
            ! define prefix_session_index used for defining the loading filename of the session index
            parameters % prefix_session_index= 'session_index1_'
        
            ! define path for saving the new session index
            parameters % save_path_session_index= '../DATA/SESSION_INDEX1/'
        
            ! define path for saving the new epolog textfiles
            parameters % save_path_epologs= '../DATA/EPOLOG1/'
        
        case ('two_epochs_per_obs')
        
            ! define prefix_session_index used for defining the loading filename of the session index
            parameters % prefix_session_index= 'session_index2_'
        
            ! define path for saving the new session index2
            parameters % save_path_session_index= '../DATA/SESSION_INDEX2/'
        
            ! define path for saving the new epolog textfiles
            parameters % save_path_epologs= '../DATA/EPOLOG2/'
            
        case default ! case if other cases fail
            ! report error message
            write(unit= *, fmt= '(a)') 'Error: Set correct epolog mode! Program stopped!'
            ! stop the program
            stop
            
        end select
    
    
    ! define loading path to the stored session indices
    parameters % load_path_session_index= parameters % save_path_session_index
    
    ! define loading path to the stored epologs
    parameters % load_path_epologs= parameters % save_path_epologs
    
    ! results .radiate- and trp-files
    !
    ! define path to the output directory, where the .radiate-files should be stored
    parameters % save_path_radiate= '../RESULTS/RADIATE/'
    
    ! define path to the output directory, where the .trp-files should be stored
    parameters % save_path_trp= '../RESULTS/TRP/'
    
    
    ! optional parameter
    !
    ! depending on the parameter define path for saving the error log data
    if (parameters % save_errorlog) then
        ! define path for saving the log entries for each epolog
        parameters % save_path_errorlog= '../DATA/ERROR_LOG/'
    end if
    
    
    !---------------------------------------------------------------------------------------------------
    
    ! report the beginning of data handling and processing
    write(unit= *, fmt= '(/, a)') '=============='
    write(unit= *, fmt= '(a)') '= Processing ='
    write(unit= *, fmt= '(a, /)') '=============='
    
    
    !---------------------------------------------------------------------------------------------------
        
    ! Stations - read in station dependent values (station catalog)
    
    ! report process
    write(unit= *, fmt= '(a)') 'Reading the station catalog ...'
    
    ! loading from textfile
    call import_station_data ( parameters % load_path_stat_info, &
                               parameters % load_filename_stat_info, &
                               name_list, &
                               lat_ell, &
                               lon_ell, &
                               h_ell )
    
    ! report elapsed time
    call get_elapsed_time( icount_start, count_rate, icount_interm, elapsed_time )
    write(unit= *, fmt= '(a, f0.3, a, /)') 'Total elapsed time after reading the station catalog: ', elapsed_time, ' s'
    
    
    !---------------------------------------------------------------------------------------------------
    
    ! Start processing the input session
    
    ! Create grib-file-epoch-dependent observation-log-files (epologs) and session index
    
    ! report process
    write(unit= *, fmt= '(a)') 'Creating epologs and session index ...'
    
    ! This section splits the input azel-file (=session) with all observations for each station and
    ! creates new grib-file-epoch-dependent observation-log-files called "epolog".
    
    ! determine the mode of assigning the observations to epochs (select method of creating epologs)
    select case (parameters % epolog_mode)
        ! case of one_epoch_per_obs mode
        case ('one_epoch_per_obs')
            ! run the subroutine to create the epologs for the input azel-file
            call create_session_index_and_epologs_single( parameters % load_path_azel, &
                                                          parameters % load_filename_azel, &
                                                          parameters % save_path_session_index, &
                                                          parameters % save_path_epologs, &
                                                          parameters % readAzel, &
                                                          parameters % load_path_filename_indAzelSpec, &
                                                          size(name_list), &
                                                          name_list ) 
        
        ! case of two_epochs_per_obs mode
        case ('two_epochs_per_obs')
            ! run the subroutine to create the epologs for the input azel-file
            call create_session_index_and_epologs_dual( parameters % load_path_azel, &
                                                        parameters % load_filename_azel, &
                                                        parameters % save_path_session_index, &
                                                        parameters % save_path_epologs, &
                                                        parameters % readAzel, &
                                                        parameters % load_path_filename_indAzelSpec, &
                                                        size(name_list), &
                                                        name_list ) 
    
        case default ! Case if other cases fail. This should not happen as this case has already been checked in a prior step.
            ! report error message
            write(unit= *, fmt= '(a)') 'Error: Unexpected epolog mode encountered! Program stopped!'
            ! stop the program
            stop
            
    end select
    
    ! report elapsed time
    call get_elapsed_time( icount_start, count_rate, icount_interm, elapsed_time )
    write(unit= *, fmt= '(a, f0.3, a, /)') 'Total elapsed time after creating epologs and session index: ', elapsed_time, ' s'
    
    
    !---------------------------------------------------------------------------------------------------
         
    ! Determine which epologs should be loaded for the session using the session index textfile
    
    ! report process
    write(unit= *, fmt= '(a)') 'Loading session index ...'
    
    ! get all epolog-filenames for the session (=AZEL-file)
    
    ! load session index
    ! get all epolog-filenames necessary for the current session (azel-file)
    call load_session_index( parameters % load_path_session_index, &
                             parameters % prefix_session_index // session_name // '.txt', &
                             parameters % cleanup, &
                             epolog_filenames, &
                             nr_epologs )
    
    ! report elapsed time
    call get_elapsed_time( icount_start, count_rate, icount_interm, elapsed_time )
    write(unit= *, fmt= '(a, f0.3, a, /)') 'Total elapsed time after loading session index: ', elapsed_time, ' s'
    
    
    !---------------------------------------------------------------------------------------------------
    
    ! Load all epologs found for the current session sequentially and process the data in loops over the single epologs
    ! A loop over all epologs is used to process (ray-trace) all observations for one session.
    ! In the next loop the next epolog is processed.

    
    ! allocate the structure type for storing the epolog data (and all upcoming epolog-based data)
    ! note: allocation is done for the total number of epologs
    allocate(epolog(nr_epologs))
    
    ! allocate the errorlog for the session with size 0 to prevent errors, e.g. deallocation of unallocated variable or size() of unallocated variable
    allocate( errorlog(0) )
    
    ! loop over all epologs associated to one specific azel-file (session)
    epolog_loop: do curr_epolog_ind= 1, nr_epologs

        ! report info of processing status
        ! notes on output format:
        ! i0.1 outputs an integer with the minimal field width and at least one digit printed
        write(unit= *, fmt= '(/, a)') '============================================================================================'
        write(unit= *, fmt= '(a, a, a, i0.1, a, i0.1, a, a, a )') 'Processing "', epolog_filenames(curr_epolog_ind), '" [nr. ', curr_epolog_ind, ' of ', nr_epologs, '] for session "', session_name, '"'
        write(unit= *, fmt= '(a, /)')  '============================================================================================'
        
        ! report new processing section
        write(unit= *, fmt= '(a)') '------------------------------'
        write(unit= *, fmt= '(a)') 'Preparing the observation data'
        write(unit= *, fmt= '(a, /)') '------------------------------'
        
        
        ! report process
        write(unit= *, fmt= '(a)') 'Loading epolog...'
        
        ! load the epolog contents
        call load_epologs( parameters % load_path_epologs, &
                           epolog_filenames(curr_epolog_ind), &
                           parameters % cleanup, &
                           epolog(curr_epolog_ind) )
        
        ! report elapsed time
        call get_elapsed_time( icount_start, count_rate, icount_interm, elapsed_time )
        write(unit= *, fmt= '(a, f0.3, a, /)') 'Total elapsed time after loading epolog: ', elapsed_time, ' s'
        
        !---------------------------------------------------------------------------------------------------
    
        ! Determine the observing stations and their number of observations in the current epolog
        ! This section gets the unique names of each station observing in the current epolog.
        ! Additionally the number of observations by each station is counted. So each station with
        ! at least one appearance in the observations and the number of observations per station are
        ! determined.
        ! note: The observations in the epolog are already sorted alphabetically (as required for the subroutine "get_observing_stations")
        !       after the station name as this has been done during the epolog-file creation.
    
        ! report process
        write(unit= *, fmt= '(a)') 'Determining names of observing stations and their number of observations ...'
        
        ! call subroutine to get observing stations and number of observations per station
        ! note: assigning of the observing station names and their observation number directly to the
        !       "observing_stations" substructure is not possible, because automatic allocation is not supported
        !       in an assignment of struct(i) % field= vector(i)
        !       Therefore the subroutine output has to be stored in temporary vectors and the substructure needs
        !       to be allocated before the final assignment of the values to the substructure!
        call get_observing_stations( epolog(curr_epolog_ind) % nr_obs, &
                                     epolog(curr_epolog_ind) % observations % station, &
                                     obs_stations, &
                                     nr_obs_per_obsstat, &
                                     epolog(curr_epolog_ind) % nr_observing_stations )
        
        ! allocate the substructure "observing_stations"
        allocate(epolog(curr_epolog_ind) % observing_stations(epolog(curr_epolog_ind) % nr_observing_stations))
        
        ! assign the values of the vectors "obs_stations" and "nr_obs_per_obsstat" to the substructure "observing_stations"
        epolog(curr_epolog_ind) % observing_stations % station= obs_stations
        epolog(curr_epolog_ind) % observing_stations % nr_obs_per_obsstat= nr_obs_per_obsstat
        
        ! deallocate the temporary vectors
        deallocate(obs_stations, nr_obs_per_obsstat)
        
        ! report elapsed time
        call get_elapsed_time( icount_start, count_rate, icount_interm, elapsed_time )
        write(unit= *, fmt= '(a, f0.3, a, /)') 'Total elapsed time after determining names of observing stations and their number of observations: ', elapsed_time, ' s'
        
        
        !---------------------------------------------------------------------------------------------------
        
        ! Initialize variable for each observing station for checking if ray-tracing should be suspended
        
        ! note: This is already done in the default initialization in the type definition (see module for type definitions)
        
        ! initialize variable for checking if ray-tracing should be suspended for a specific station
        ! due to station coordinates (lat, lon, height) that are not supported through the grid or the
        ! maximum interpolated height level
        ! .FALSE. ... ray-tracing will be done for this specific station
        ! .TRUE. ... ray-tracing will be suspended for this specific station
        
        ! note: This is already done in the default initialization in the type definition (see module for type definitions)
        !epolog(curr_epolog_ind) % observing_stations % suspend_raytr= .FALSE.
        
        
        !---------------------------------------------------------------------------------------------------
        
        ! Get coordinates of the observing stations from the station catalogue
        ! This section determines the coordinates of the observing stations in the current epolog.
        ! In case that not all stations listed with observations in the epolog can be found in the station
        ! catalogue in order to derive the position information, a message is shown in the command window
        ! explaining which station is missing. All encontered problems are saved to an error log. If
        ! more than one matching station is found in the catalogue also a warning message is
        ! created.
        
        ! report process
        write(unit= *, fmt= '(a)') 'Getting coordinates for the observing stations ...'
        
        ! call subroutine to get station coordinates from station catalogue by comparing station names
        call get_coord_observing_stations( epolog(curr_epolog_ind) % nr_observing_stations, &
                                           epolog(curr_epolog_ind) % observing_stations % station, &
                                           size(name_list), &
                                           name_list, &
                                           lat_ell, &
                                           lon_ell, &
                                           h_ell, &
                                           epolog(curr_epolog_ind) % observing_stations % lat_ell, &
                                           epolog(curr_epolog_ind) % observing_stations % lon_ell, &
                                           epolog(curr_epolog_ind) % observing_stations % h_ell, &
                                           errorlog_curr_epolog_obsstat_coord, &
                                           epolog(curr_epolog_ind) % observing_stations % suspend_raytr)
        
        ! determine the number of errors in "errorlog_curr_epolog_obsstat_coord"
        ! note: size() without dimension specification returns the total number total elements in the array
        nr_curr_errors= size(errorlog_curr_epolog_obsstat_coord)
        
        ! check if any error has occurred
        if (nr_curr_errors /= 0) then
            
            ! determine the current size of "errorlog"
            ! note: if " errorlog" has not been allocated before size() will returns 0 as desired
            ! note: size() without dimension specification returns the total number total elements in the array
            nr_total_errors= size(errorlog)
            
            ! determine the new size to which "errorlog" will be resized
            new_nr_total_errors= nr_total_errors + nr_curr_errors
        
            ! call subroutine to resize "errorlog"
            call resize_errorlog( errorlog, new_nr_total_errors )
            
            ! assign the errors for missing station coordinates in the current epolog to the general errorlog
            errorlog(nr_total_errors + 1 : new_nr_total_errors)= errorlog_curr_epolog_obsstat_coord
        
        end if
        
        ! note: as possible errors of missing station coordinates are different for each epolog the epolog-specific errorlog-variable has to be deallocated
        !       after its usage
        deallocate( errorlog_curr_epolog_obsstat_coord )
        
        ! report elapsed time
        call get_elapsed_time( icount_start, count_rate, icount_interm, elapsed_time )
        write(unit= *, fmt= '(a, f0.3, a, /)') 'Total elapsed time after getting coordinates for the observing stations: ', elapsed_time, ' s'
        
        
        !---------------------------------------------------------------------------------------------------
                
        ! Get meteorological values - get global grid and improve vertical resolution
        ! This section is necessary to retrieve profiles with high (interpolated) vertical
        ! resolution for the whole global grid.
        
        ! report new processing section
        write(unit= *, fmt= '(/, a)') '---------------------------------'
        write(unit= *, fmt= '(a)') 'Preparing the meteorological data'
        write(unit= *, fmt= '(a, /)') '---------------------------------'
        
        ! get grib-file name without extension
        ! note: The filename without extension is equal to the epoch name.
        epolog(curr_epolog_ind) % grib_filename= epolog(curr_epolog_ind) % epoch_name
        
        ! allocate the substructure "profile", which stores the grid info, the meteorologic profiles (from grib-file, refined (interpolated/extrapolated)) and the geoid undulations
        ! note: "profile" is needed to be allocatable, because at the end of processing the current epolog its "profile" part is deallocated to free RAM
        allocate(epolog(curr_epolog_ind) % profile)
        
        ! get the global meteorologic profiles
        ! Note: Read in grib-file must always contain only one epoch.
        call get_meteo_undulation_refr_glob( parameters % load_path_grib, &
                                             epolog(curr_epolog_ind) % grib_filename, &
                                             epolog(curr_epolog_ind) % observing_stations % lat_ell, &
                                             epolog(curr_epolog_ind) % observing_stations % lon_ell, &
                                             epolog(curr_epolog_ind) % observing_stations % h_ell, &
                                             epolog(curr_epolog_ind) % observing_stations % station, &
                                             parameters % load_path_undulation, &
                                             parameters % interpolation_method, &
                                             epolog(curr_epolog_ind) % observing_stations % suspend_raytr, &
                                             epolog(curr_epolog_ind) % profile, &
                                             epolog(curr_epolog_ind) % epoch )
        
        ! report elapsed time
        call get_elapsed_time( icount_start, count_rate, icount_interm, elapsed_time )
        write(unit= *, fmt= '(a, f0.3, a, /)') 'Total elapsed time after preparing the meteorological data: ', elapsed_time, ' s'
        
        
        !---------------------------------------------------------------------------------------------------
        
        ! Define ray path and calculate refractivity, elevation and delay
        ! This section is necessary to define the intersection points of the ray with the interpolated
        ! height levels.
        ! In order to determine the intersection points one of the following ray-tracing approaches is
        ! used together with the observed outgoing elevation angle:
        ! piece-wise linear, refined piece-wise linear or Thayer approximation 
        
        ! At the intersection points the total refractive indices will be be used to determine the
        ! delay along the path.
    
        ! report new processing section
        write(unit= *, fmt= '(/, a)') '-----------'
        write(unit= *, fmt= '(a)') 'Ray-tracing'
        write(unit= *, fmt= '(a, /)') '-----------'
        
        ! store for each observation the mjd of the grib-file epoch at which epoch the delay has been calculated
        !logs.(curr_epolog).observations.grib_epoch_mjd(1:logs.(curr_epolog).observations.nr_obs,1)=logs.(curr_epolog).epoch.mjd;
        
        ! allocate sub-structure "delay" for storing the calculated delays, elevation angles, mapping factors and possible iteration breaks
        allocate( epolog(curr_epolog_ind) % delay(epolog(curr_epolog_ind) % nr_obs) )
        
        ! allocate sub-structure "meteo_stat_out" for storing meteorological data interpolated at the station position for the output to the .radiate file
        allocate( epolog(curr_epolog_ind) % meteo_stat_out(epolog(curr_epolog_ind) % nr_obs) )
        
        ! determine the NaN-value as signaling NaN (producing exceptions when used in calculation)
        NaN= ieee_value(0.0d0, ieee_signaling_nan) ! 0.0d0 specifies double precision kind
        
        ! initialize variables in the substructure "delay"
        ! note: This step is needed only for observations from suspended stations that would otherwise have a random value written to the variables.
        !       The setting to NaN will be done in a later part in a positive case of suspending ray-tracing.
        !       For all other really ray-traced variables the initial values will be overwritten anyway by the results.
        ! note: Initilization of the logical values for break-determination in iteration to .FALSE. is already done in the type definition.
        !epolog(curr_epolog_ind) % delay % break_elev= .FALSE. ! already done in type definition, this is necessary in case that specific observations are not treated in ray-tracing and therefore the value would remain unset
        !epolog(curr_epolog_ind) % delay % break_layer= .FALSE. ! already done in type definition, this is necessary in case that specific observations are not treated in ray-tracing with Thayer method and therefore the value would remain unset
                
        
        ! report elapsed time
        call get_elapsed_time( icount_start, count_rate, icount_interm, elapsed_time )
        write(unit= *, fmt= '(a, f0.3, a, /)') 'Total elapsed time at the start of the ray-tracing for all stations and their observations: ', elapsed_time, ' s'
        
        ! check if method "pwl" should be used for ray-tracing (for whole sessions and all stations)
        select case (parameters % raytr_method)
            
            ! case of piece-wise linear approach
            case ('pwl')
                
                ! report process
                write(unit= *, fmt= '(tr4, a)') 'Ray-tracing method "pwl" requires mean refractive indices built from values of two consecutive levels used for ray-tracing in the intermediate space between the levels.'
                write(unit= *, fmt= '(tr4, a)') 'Start calculating mean refractive indices for the first ray path (first intermediate space above the station) ...'
                
                ! calculate the mean refractive index just for the level between station height
                ! and first interpolated height level above the station is calculated for all stations
                ! (loop over the stations)
                ! attention: this step has to be done prior to calculating the mean refractive indices for
                ! all interpolated height levels, because their values will be overwritten due to RAM
                ! reasons
                
                ! loop over all observing stations
                do ind_stat= 1, epolog(curr_epolog_ind) % nr_observing_stations
                    
                    ! check if ray-tracing can be done for the current station (suspend_raytr= .FALSE. ... ray-tracing, .TRUE. ... no ray-tracing)
                    ! e.g. POI has valid coordinates (lat, lon, height)
                    if ( .NOT. epolog(curr_epolog_ind) % observing_stations(ind_stat) % suspend_raytr ) then
                        
                        ! calculate the first mean refractive index between station level and the 
                        ! first height level above the station with the index "start_lev"
        
                        ! get value of refractive index at the station's horizontal position of the
                        ! height level above the station
                        ! hydrostatic, wet and total refractive index
                        ! note: see subroutine for optional arguments, use keywords for all arguments as specified in the subroutine to get correct call if some optional arguments are missing in between! 
                        call get_bilint_value( v1= epolog(curr_epolog_ind) % profile % meteo_int % n_h, &
                                               v2= epolog(curr_epolog_ind) % profile % meteo_int % n_w, &
                                               v3= epolog(curr_epolog_ind) % profile % meteo_int % n_total, &
                                               level= epolog(curr_epolog_ind) % profile % meteo_stat(ind_stat) % start_lev, &
                                               POI_lat= epolog(curr_epolog_ind) % observing_stations(ind_stat) % lat_ell, &
                                               POI_lon= epolog(curr_epolog_ind) % observing_stations(ind_stat) % lon_ell, &
                                               dint_lat= epolog(curr_epolog_ind) % profile % grid % res_lat, &
                                               dint_lon= epolog(curr_epolog_ind) % profile % grid % res_lon, &
                                               ind_lat1lon1_in= epolog(curr_epolog_ind) % profile % meteo_stat(ind_stat) % ind_lat1lon1, & ! input of indices to avoid the call of subroutine to determine them as indices have already been determined in the subroutine "calc_refr_ind_at_stations"
                                               ind_lat1lon2_in= epolog(curr_epolog_ind) % profile % meteo_stat(ind_stat) % ind_lat1lon2, & ! grid_lat, grid_lon, grid_size and start_and_global_check are therefore not needed
                                               ind_lat2lon2_in= epolog(curr_epolog_ind) % profile % meteo_stat(ind_stat) % ind_lat2lon2, &
                                               ind_lat2lon1_in= epolog(curr_epolog_ind) % profile % meteo_stat(ind_stat) % ind_lat2lon1, &
                                               v1_bilint= n_h_stat_upper, &
                                               v2_bilint= n_w_stat_upper, &
                                               v3_bilint= n_total_stat_upper )
                                               
                                               ! note: usage of optional input of indices lat1lon1 etc. as values have already been determined in the subroutine "calc_refr_ind_at_stations"
                                               ! repeated output of these indices is therefore not necessary
                        
                        ! overwrite the original refractivity values at the station height and position with mean values
                        ! calculated from the original values and the ones from the upper level at the station position
                        ! hydrostatic refractive index
                        epolog(curr_epolog_ind) % profile % meteo_stat(ind_stat) % n_h= (epolog(curr_epolog_ind) % profile % meteo_stat(ind_stat) % n_h + n_h_stat_upper) / 2
                        ! wet refractive index
                        epolog(curr_epolog_ind) % profile % meteo_stat(ind_stat) % n_w= (epolog(curr_epolog_ind) % profile % meteo_stat(ind_stat) % n_w + n_w_stat_upper) / 2
                        ! total refractive index
                        epolog(curr_epolog_ind) % profile % meteo_stat(ind_stat) % n_total= (epolog(curr_epolog_ind) % profile % meteo_stat(ind_stat) % n_total + n_total_stat_upper) / 2
                        
                    end if ! end of if-expressing for determinining if ray-tracing is possible
                    
                end do! end of for-loop over all stations
                
                ! report elapsed time
                call get_elapsed_time( icount_start, count_rate, icount_interm, elapsed_time )
                write(unit= *, fmt= '(tr4, a, f0.3, a, /)') 'Total elapsed time after calculating mean refractive indices for the first intermediate space above the station: ', elapsed_time, ' s'
                
                ! report process
                write(unit= *, fmt= '(tr4, a)') 'Start calculating mean refractive indices for the other intermediate spaces between the levels ...'
                
                ! mean refractive index calculation is done for all height levels interpolated
                ! from the original grib-file.
                ! note: the new mean refractive indices at the station level from the original station level value and the
                !       above interpolated height level has to be done before (as done above), because the original values
                !       will be overwritten in the following following steps!
        
                ! get number of height levels
                nr_lev= epolog(curr_epolog_ind) % profile % meteo_int % nr_h_lev
                
                ! calculate mean refractive indices for all height levels from two consecutive
                ! height levels, which will be used for ray tracing in "pwl" approach
                
                ! note: No temporal variable for temporal storage of calculated mean refractive indices needed
                !       as mean values will always be saved to an index level that will not be used again during the mean calculation!
                !       The correct output has been tested.
                
                ! calculate mean refractive indices for one level per loop cycle
                ! loop from level 1 to nr_lev-1 as i+1 is needed inside the loop and the mean building delivers nr_lev-1 levels with mean values
                do i= 1, nr_lev-1
                    
                    ! calculate mean for hydrostatic refractive index
                    ! note: the current height level and the next height level build the mean and this will be assigned to the index of the current level index,
                    !       so no overwriting of orignal data that would be needed again happens
                    epolog(curr_epolog_ind) % profile % meteo_int % n_h(i, :, :)= (epolog(curr_epolog_ind) % profile % meteo_int % n_h(i, :, :) + epolog(curr_epolog_ind) % profile % meteo_int % n_h(i + 1, :, :)) / 2
                    
                    ! calculate mean for wet refractive index
                    ! note: the current height level and the next height level build the mean and this will be assigned to the index of the current level index,
                    !       so no overwriting of orignal data that would be needed again happens
                    epolog(curr_epolog_ind) % profile % meteo_int % n_w(i, :, :)= (epolog(curr_epolog_ind) % profile % meteo_int % n_w(i, :, :) + epolog(curr_epolog_ind) % profile % meteo_int % n_w(i + 1, :, :)) / 2
                    
                    ! calculate mean for total refractive index
                    ! note: the current height level and the next height level build the mean and this will be assigned to the index of the current level index,
                    !       so no overwriting of orignal data that would be needed again happens
                    epolog(curr_epolog_ind) % profile % meteo_int % n_total(i, :, :)= (epolog(curr_epolog_ind) % profile % meteo_int % n_total(i, :, :) + epolog(curr_epolog_ind) % profile % meteo_int % n_total(i + 1, :, :)) / 2
                    
                end do
                
                ! change last level values (not changed by mean calculation), which stand for the value of
                ! space to n= 1. This is necessary because the loop of calculating the ray trace via "pwl" needs the index
                ! i+1 (in the last step this means the value of nr_lev).
                epolog(curr_epolog_ind) % profile % meteo_int % n_h(nr_lev, :, :)= 1
                epolog(curr_epolog_ind) % profile % meteo_int % n_w(nr_lev, :, :)= 1
                epolog(curr_epolog_ind) % profile % meteo_int % n_total(nr_lev, :, :)= 1
        
                ! report elapsed time
                call get_elapsed_time( icount_start, count_rate, icount_interm, elapsed_time )
                write(unit= *, fmt= '(tr4, a, f0.3, a, /)') 'Total elapsed time after calculating mean refractive indices for the other intermediate spaces between the levels: ', elapsed_time, ' s'
                
        end select ! end of "pwl"-query
        
        
        ! set index for end of observations
        ! note: setting "obs_end" to 0 at the beginning of station_loop is necessary as it is used when determining "obs_begin"
        obs_end= 0
        
        ! for the current epoch a loop over all observing stations is used and another loop over each
        ! observation of the specific station
        station_loop: do ind_stat= 1, epolog(curr_epolog_ind) % nr_observing_stations
            
            ! determine index for the beginning of the station's observations
            obs_begin= obs_end + 1
            
            ! get index for end of observations of the current observing station
            obs_end= obs_end + epolog(curr_epolog_ind) % observing_stations(ind_stat) % nr_obs_per_obsstat
            
            ! get name of current station
            curr_stat= epolog(curr_epolog_ind) % observing_stations(ind_stat) % station
            
            ! write current status of ray-tracing (which station is processed)
            write(unit= *, fmt= '(/, tr4, a)') '.................................................................................'
            write(unit= *, fmt= '(tr4, 5a)') 'Ray-tracing for station "', trim(curr_stat), '" within "', epolog_filenames(curr_epolog_ind), '"'
            write(unit= *, fmt= '(tr4, a)') '.................................................................................'
            
            ! check if ray-tracing can be done (suspend_raytr= .FALSE. ... ray-tracing, .TRUE. ... no ray-tracing)
            ! e.g. POI has valid coordinates (lat, lon, height)
            if ( .NOT. epolog(curr_epolog_ind) % observing_stations(ind_stat) % suspend_raytr ) then
                
                ! loop over each observation of the current station
                obs_loop: do r= obs_begin, obs_end
                    
                    ! index of current observation of the current station
                    obs_ind_curr_stat= r - obs_begin + 1
                    
                    ! write current status of ray-tracing (which observation is processed) after every 10
                    ! observations
                    ! notes on output format:
                    ! i0.1 outputs an integer with the minimal field width and at least one digit printed
                    ! note: modulo() in Fortran works like mod() in MATLAB (if non-zero result then the sign of the second argument is retained)
                    if ( modulo(obs_ind_curr_stat,10) == 0 ) then
                        write(unit= *, fmt= '(tr8, a, i0.1, a, i0.1, a)') 'Delays of ', obs_ind_curr_stat, ' of ', epolog(curr_epolog_ind) % observing_stations(ind_stat) % nr_obs_per_obsstat, ' observations processed.'
                    end if
                    
                    
                    ! determine the ray-tracing method which is selected for estimating the intersection
                    ! points and the delays
                    select case (parameters % raytr_method)
                        ! case of piece-wise linear approach
                        case ('pwl')
                            
                            call get_RayTrace2D_pwl_global( epolog(curr_epolog_ind) % observations(r) % az, &
                                                            epolog(curr_epolog_ind) % observations(r) % elev, &
                                                            epolog(curr_epolog_ind) % observing_stations(ind_stat) % lat_ell, &
                                                            epolog(curr_epolog_ind) % observing_stations(ind_stat) % lon_ell, &
                                                            epolog(curr_epolog_ind) % observing_stations(ind_stat) % h_ell, &
                                                            epolog(curr_epolog_ind) % profile % meteo_stat(ind_stat) % ind_lat1lon1, &
                                                            epolog(curr_epolog_ind) % profile % meteo_stat(ind_stat) % ind_lat1lon2, &
                                                            epolog(curr_epolog_ind) % profile % meteo_stat(ind_stat) % ind_lat2lon2, &
                                                            epolog(curr_epolog_ind) % profile % meteo_stat(ind_stat) % ind_lat2lon1, &
                                                            epolog(curr_epolog_ind) % profile % meteo_int % nr_h_lev, &
                                                            epolog(curr_epolog_ind) % profile % meteo_int % h, &
                                                            epolog(curr_epolog_ind) % profile % meteo_stat(ind_stat) % start_lev, &
                                                            epolog(curr_epolog_ind) % profile % grid % grid_lat, &
                                                            epolog(curr_epolog_ind) % profile % grid % grid_lon, &
                                                            epolog(curr_epolog_ind) % profile % grid % grid_size, &
                                                            epolog(curr_epolog_ind) % profile % grid % res_lat, &
                                                            epolog(curr_epolog_ind) % profile % grid % res_lon, &
                                                            epolog(curr_epolog_ind) % profile % meteo_stat(ind_stat) % n_total, &
                                                            epolog(curr_epolog_ind) % profile % meteo_stat(ind_stat) % n_h, &
                                                            epolog(curr_epolog_ind) % profile % meteo_stat(ind_stat) % n_w, &
                                                            epolog(curr_epolog_ind) % profile % meteo_int % n_total, &
                                                            epolog(curr_epolog_ind) % profile % meteo_int % n_h, &
                                                            epolog(curr_epolog_ind) % profile % meteo_int % n_w, &
                                                            epolog(curr_epolog_ind) % profile % grid % start_and_global_check, &
                                                            epolog(curr_epolog_ind) % delay(r) )
                            
                        ! case of refined piece-wise linear approach
                        case ('ref_pwl')
                            
                            call get_RayTrace2D_ref_pwl_global( epolog(curr_epolog_ind) % observations(r) % az, &
                                                                epolog(curr_epolog_ind) % observations(r) % elev, &
                                                                epolog(curr_epolog_ind) % observing_stations(ind_stat) % lat_ell, &
                                                                epolog(curr_epolog_ind) % observing_stations(ind_stat) % lon_ell, &
                                                                epolog(curr_epolog_ind) % observing_stations(ind_stat) % h_ell, &
                                                                epolog(curr_epolog_ind) % profile % meteo_stat(ind_stat) % ind_lat1lon1, &
                                                                epolog(curr_epolog_ind) % profile % meteo_stat(ind_stat) % ind_lat1lon2, &
                                                                epolog(curr_epolog_ind) % profile % meteo_stat(ind_stat) % ind_lat2lon2, &
                                                                epolog(curr_epolog_ind) % profile % meteo_stat(ind_stat) % ind_lat2lon1, &
                                                                epolog(curr_epolog_ind) % profile % meteo_int % nr_h_lev, &
                                                                epolog(curr_epolog_ind) % profile % meteo_int % h, &
                                                                epolog(curr_epolog_ind) % profile % meteo_stat(ind_stat) % start_lev, &
                                                                epolog(curr_epolog_ind) % profile % grid % grid_lat, &
                                                                epolog(curr_epolog_ind) % profile % grid % grid_lon, &
                                                                epolog(curr_epolog_ind) % profile % grid % grid_size, &
                                                                epolog(curr_epolog_ind) % profile % grid % res_lat, &
                                                                epolog(curr_epolog_ind) % profile % grid % res_lon, &
                                                                epolog(curr_epolog_ind) % profile % meteo_stat(ind_stat) % n_total, &
                                                                epolog(curr_epolog_ind) % profile % meteo_stat(ind_stat) % n_h, &
                                                                epolog(curr_epolog_ind) % profile % meteo_stat(ind_stat) % n_w, &
                                                                epolog(curr_epolog_ind) % profile % meteo_int % n_total, &
                                                                epolog(curr_epolog_ind) % profile % meteo_int % n_h, &
                                                                epolog(curr_epolog_ind) % profile % meteo_int % n_w, &
                                                                epolog(curr_epolog_ind) % profile % grid % start_and_global_check, &
                                                                epolog(curr_epolog_ind) % delay(r) )
                            
                        ! case of Thayer approach    
                        case ('Thayer')
                            
                            call get_RayTrace2D_Thayer_global( epolog(curr_epolog_ind) % observations(r) % az, &
                                                               epolog(curr_epolog_ind) % observations(r) % elev, &
                                                               epolog(curr_epolog_ind) % observing_stations(ind_stat) % lat_ell, &
                                                               epolog(curr_epolog_ind) % observing_stations(ind_stat) % lon_ell, &
                                                               epolog(curr_epolog_ind) % observing_stations(ind_stat) % h_ell, &
                                                               epolog(curr_epolog_ind) % profile % meteo_stat(ind_stat) % ind_lat1lon1, &
                                                               epolog(curr_epolog_ind) % profile % meteo_stat(ind_stat) % ind_lat1lon2, &
                                                               epolog(curr_epolog_ind) % profile % meteo_stat(ind_stat) % ind_lat2lon2, &
                                                               epolog(curr_epolog_ind) % profile % meteo_stat(ind_stat) % ind_lat2lon1, &
                                                               epolog(curr_epolog_ind) % profile % meteo_int % nr_h_lev, &
                                                               epolog(curr_epolog_ind) % profile % meteo_int % h, &
                                                               epolog(curr_epolog_ind) % profile % meteo_stat(ind_stat) % start_lev, &
                                                               epolog(curr_epolog_ind) % profile % grid % grid_lat, &
                                                               epolog(curr_epolog_ind) % profile % grid % grid_lon, &
                                                               epolog(curr_epolog_ind) % profile % grid % grid_size, &
                                                               epolog(curr_epolog_ind) % profile % grid % res_lat, &
                                                               epolog(curr_epolog_ind) % profile % grid % res_lon, &
                                                               epolog(curr_epolog_ind) % profile % meteo_stat(ind_stat) % n_total, &
                                                               epolog(curr_epolog_ind) % profile % meteo_stat(ind_stat) % n_h, &
                                                               epolog(curr_epolog_ind) % profile % meteo_stat(ind_stat) % n_w, &
                                                               epolog(curr_epolog_ind) % profile % meteo_int % n_total, &
                                                               epolog(curr_epolog_ind) % profile % meteo_int % n_h, &
                                                               epolog(curr_epolog_ind) % profile % meteo_int % n_w, &
                                                               epolog(curr_epolog_ind) % profile % grid % start_and_global_check, &
                                                               epolog(curr_epolog_ind) % delay(r) )
                            
                        ! case if other cases fail
                        case default
                            ! report error message
                            write(unit= *, fmt= '(/, a)') 'Error: Unexpected ray-tracing method! No delay calculation executed. Please set correct ray-tracing approach. Program stopped!'
                            ! stop the program
                            stop
                            
                        end select
                    
                    ! display warning if break occured during ray-tracing due to loop
                    ! overflow of elevation angle iteration
                        
                    ! report if desired accuracy for outgoing elevation angle couldn't be
                    ! reached (iteration has been terminated due to reaching max. loop number)
                    ! note: Use .eqv. to check for == in case of logical values as gfortran does not support == for logicals.
                    if (epolog(curr_epolog_ind) % delay(r) % break_elev .eqv. .TRUE.) then
                        
                        ! define error description
                        ! i0.1 outputs an integer with the minimal field width and at least one digit printed
                        ! "es" is the exponential format, "sp" prints plus sign
                        write(unit= error_description, fmt= '(a, i0.1, a, a, a, a, a, sp, es12.5)') 'For scannumber "', &
                                                                                                    epolog(curr_epolog_ind) % observations(r) % scannr, &
                                                                                                    '", station "', &
                                                                                                    curr_stat, &
                                                                                                    '" in "', &
                                                                                                    epolog_filenames(curr_epolog_ind), &
                                                                                                    '" the desired accuracy of the outgoing elevation angle couldn''t be reached! Resulting difference [outgoing (theoretical) - outgoing (ray-traced) elevation angle] in [rad]: ', &
                                                                                                    epolog(curr_epolog_ind) % delay(r) % diff_e
                        
                        ! report message
                        ! note: as "error_description" is initialized with a too long length trim() removes trailing blanks
                        write(unit= *, fmt= '(tr12, a, a)') 'Warning: ', trim(error_description)
                        
                        ! determine the current size of "errorlog"
                        ! note: if " errorlog" has not been allocated before size() will returns 0 as desired
                        ! note: size() without dimension specification returns the total number total elements in the array
                        nr_total_errors= size(errorlog)
            
                        ! determine the new size to which "errorlog" will be resized = index of new error
                        ! former total size + 1 as one error should be added
                        ind_error= nr_total_errors + 1
        
                        ! call subroutine to resize "errorlog"
                        call resize_errorlog( errorlog, ind_error )
                        
                        ! assign the error to the general errorlog as in the current observation a break
                        ! occured during ray-tracing due to loop overflow of elevation angle iteration
                        ! note: as "error_description" is initialized with a too long length trim() removes trailing blanks
                        errorlog(ind_error) % error_nr = ind_error
                        errorlog(ind_error) % error_type= 'Break in outgoing elevation angle iteration'
                        errorlog(ind_error) % error_descr= trim(error_description)
                        
                    end if
                    
                    ! case of Thayer or other approach where possible:
                    ! report if desired accuracy of the position for at least one of the intersection points along the ray
                    ! couldn't be reached (iteration has been terminated due to reaching max. loop number)
                    ! note: In principle this check is only necessary in case of Thayer approach or ray-tracing approaches where this error can occur.
                    !       Pwl or refined pwl can not create a true value for this parameter.
                    ! note: Use .eqv. to check for == in case of logical values as gfortran does not support == for logicals.
                    if (epolog(curr_epolog_ind) % delay(r) % break_layer .eqv. .TRUE.) then
                        
                        ! define error description
                        ! i0.1 outputs an integer with the minimal field width and at least one digit printed
                        ! "es" is the exponential format, "sp" prints plus sign
                        write(unit= error_description, fmt= '(a, i0.1, a, a, a, a, a, sp, es12.5)') 'For scannumber "', &
                                                                                                    epolog(curr_epolog_ind) % observations(r) % scannr, &
                                                                                                    '", station "', &
                                                                                                    curr_stat, &
                                                                                                    '" in "', &
                                                                                                    epolog_filenames(curr_epolog_ind), &
                                                                                                    '" the desired position accuracy for at least one intersection point couldn''t be reached! Difference [outgoing (theoretical) - outgoing (ray-traced) elevation angle] in [rad]: ', &
                                                                                                    epolog(curr_epolog_ind) % delay(r) % diff_e
                        
                        ! report message
                        ! note: as "error_description" is initialized with a too long length trim() removes trailing blanks
                        write(unit= *, fmt= '(tr12, a, a)') 'Warning: ', trim(error_description)
                        
                        ! determine the current size of "errorlog"
                        ! note: if " errorlog" has not been allocated before size() will returns 0 as desired
                        ! note: size() without dimension specification returns the total number total elements in the array
                        nr_total_errors= size(errorlog)
            
                        ! determine the new size to which "errorlog" will be resized = index of new error
                        ! former total size + 1 as one error should be added
                        ind_error= nr_total_errors + 1
        
                        ! call subroutine to resize "errorlog"
                        call resize_errorlog( errorlog, ind_error )
                        
                        ! assign the error to the general errorlog as in the current observation a break
                        ! occured during ray-tracing due to loop overflow of intersection point position iteration
                        ! note: as "error_description" is initialized with a too long length trim() removes trailing blanks
                        errorlog(ind_error) % error_nr = ind_error
                        errorlog(ind_error) % error_type= 'At least one break in intersection point position iteration'
                        errorlog(ind_error) % error_descr= trim(error_description)
                        
                    end if
                    
                    
                    ! display information when ray-tracing for current station has been finished
                    if ( obs_ind_curr_stat == epolog(curr_epolog_ind) % observing_stations(ind_stat) % nr_obs_per_obsstat ) then
                        write(unit= *, fmt= '(tr8, a, i0.1, a)') 'Delays of all ', epolog(curr_epolog_ind) % observing_stations(ind_stat) % nr_obs_per_obsstat, ' observations processed. Finished!'
                    end if
                    
                end do obs_loop
                        
            else ! in case of suspending the ray-tracing for a specific station
        
                ! initialize the "delay"-structure for the suspended observations
                ! note: Setting the variables for not ray-traced observations in "delay" structure to NaN is done here
                !       for all observations of the suspended station instead of initializing the variables for all stations and their observations to NaN!
                ! attention: Treatment of raytraces structure elements for observations that were not processed need to be done here!
                epolog(curr_epolog_ind) % delay(obs_begin: obs_end) % dz_total= NaN
                epolog(curr_epolog_ind) % delay(obs_begin: obs_end) % dz_h= NaN
                epolog(curr_epolog_ind) % delay(obs_begin: obs_end) % dz_w= NaN
                epolog(curr_epolog_ind) % delay(obs_begin: obs_end) % ds_total_geom= NaN
                epolog(curr_epolog_ind) % delay(obs_begin: obs_end) % ds_total= NaN
                epolog(curr_epolog_ind) % delay(obs_begin: obs_end) % ds_h_geom= NaN
                epolog(curr_epolog_ind) % delay(obs_begin: obs_end) % ds_h= NaN
                epolog(curr_epolog_ind) % delay(obs_begin: obs_end) % ds_w= NaN
                epolog(curr_epolog_ind) % delay(obs_begin: obs_end) % e_stat= NaN
                epolog(curr_epolog_ind) % delay(obs_begin: obs_end) % e_outgoing_rt= NaN
                epolog(curr_epolog_ind) % delay(obs_begin: obs_end) % diff_e= NaN
                epolog(curr_epolog_ind) % delay(obs_begin: obs_end) % dgeo= NaN
                epolog(curr_epolog_ind) % delay(obs_begin: obs_end) % mf_total_geom= NaN
                epolog(curr_epolog_ind) % delay(obs_begin: obs_end) % mf_h_geom= NaN
                epolog(curr_epolog_ind) % delay(obs_begin: obs_end) % mf_w= NaN
                ! note: Initilization of the logical values for break-determination in iteration to .FALSE. is already done in the type definition.
                !epolog(curr_epolog_ind) % delay(obs_begin: obs_end) % break_elev= .FALSE. ! already done in type definition, this is necessary in case that specific observations are not treated in ray-tracing and therefore the value would remain unset
                !epolog(curr_epolog_ind) % delay(obs_begin: obs_end) % break_layer= .FALSE. ! already done in type definition, this is necessary in case that specific observations are not treated in ray-tracing with Thayer method and therefore the value would remain unset
                
                
                
                ! display warning if current observation has not been ray-traced due to suspending,
                ! e.g. due to missing station coordinates
                
                ! define error description
                ! i0.1 outputs an integer with the minimal field width and at least one digit printed
                ! "es" is the exponential format, "sp" prints plus sign
                write(unit= error_description, fmt= '(a, a, a, a, a)') 'Ray-tracing for station "', &
                                                                                            curr_stat, &
                                                                                            '" in "', &
                                                                                            epolog_filenames(curr_epolog_ind), &
                                                                                            '" has been skipped!'
                        
                ! report message
                ! note: as "error_description" is initialized with a too long length trim() removes trailing blanks
                write(unit= *, fmt= '(tr8, a, a)') 'Warning: ', trim(error_description)
                        
                ! determine the current size of "errorlog"
                ! note: if " errorlog" has not been allocated before size() will returns 0 as desired
                ! note: size() without dimension specification returns the total number total elements in the array
                nr_total_errors= size(errorlog)
            
                ! determine the new size to which "errorlog" will be resized = index of new error
                ! former total size + 1 as one error should be added
                ind_error= nr_total_errors + 1
        
                ! call subroutine to resize "errorlog"
                call resize_errorlog( errorlog, ind_error )
                        
                ! assign the error for skipped ray-tracing of a specific station to the general errorlog
                ! note: as "error_description" is initialized with a too long length trim() removes trailing blanks
                errorlog(ind_error) % error_nr = ind_error
                errorlog(ind_error) % error_type= 'Ray-tracing skipped'
                errorlog(ind_error) % error_descr= trim(error_description)
                
            end if ! end of if expression checking the suspending of ray-tracing for the current station
            
            
            ! assign the meteorological parameters at the station position to "meteo_stat_out" for the later output to the .radiate file
            ! Note: Depending on the ray-tracing settings a time interpolation of the assigned values will be carried out later.
            
            ! loop over each observation of the current station
            ! Note: In order to always have the assignments independent of wether the ray-tracing for the current station is suspended or not
            !       a separate loop over all observations is needed.
            obs_loop2: do r= obs_begin, obs_end
                
                ! pressure
                epolog(curr_epolog_ind) % meteo_stat_out(r) % p= epolog(curr_epolog_ind) % profile % meteo_stat(ind_stat) % p
                
                ! temperature
                epolog(curr_epolog_ind) % meteo_stat_out(r) % T= epolog(curr_epolog_ind) % profile % meteo_stat(ind_stat) % T
                
                ! water vapour pressure
                epolog(curr_epolog_ind) % meteo_stat_out(r) % wvpr= epolog(curr_epolog_ind) % profile % meteo_stat(ind_stat) % wvpr
                
            end do obs_loop2
            
            
        end do station_loop
        
        
        ! report elapsed time
        call get_elapsed_time( icount_start, count_rate, icount_interm, elapsed_time )
        write(unit= *, fmt= '(/, a, f0.3, a, /)') 'Total elapsed time after the ray-tracing for all stations and their observations in the current epolog: ', elapsed_time, ' s'
        
        
        ! report new processing section
        write(unit= *, fmt= '(/, a)') '--------'
        write(unit= *, fmt= '(a)') 'Clean up'
        write(unit= *, fmt= '(a, /)') '--------'
        
        ! report process
        write(unit= *, fmt= '(a)') 'Cleaning up the meteorological data in the current epolog ...'
        
        ! store resolution of grid (grib-file) to new variable as the information is needed later again, but the original structure will be deallocated
        epolog(curr_epolog_ind) % gridres % res_lat= epolog(curr_epolog_ind) % profile % grid % res_lat
        epolog(curr_epolog_ind) % gridres % res_lon= epolog(curr_epolog_ind) % profile % grid % res_lon
        
        ! delete the substructure "profile" in the "epolog"-structure that contains the grid info, the meteorologic profiles (from grib-file, refined (interpolated/extrapolated)) and the geoid undulations
        ! note: this clean up is necessary due to saving RAM, especially in cases of highly resoluted numerical weather models
        deallocate( epolog(curr_epolog_ind) % profile )
        
        ! report elapsed time
        call get_elapsed_time( icount_start, count_rate, icount_interm, elapsed_time )
        write(unit= *, fmt= '(a, f0.3, a, /)') 'Total elapsed time after cleaning up the grid-, meteorological- and geoid undulation-data in the current epolog: ', elapsed_time, ' s'
    
    end do epolog_loop ! end of epolog-based calculations (calculations for one specific epolog)
    
    
    !---------------------------------------------------------------------------------------------------
    
    ! Switch to session-based calculations
    
    ! This section is needed for cases of an epolog approach, where one observation is assigned to more than one epolog.
    ! Then time interpolation using all the delays for the same observation needs to be done in order to get one final delay at the exact observation time.
    ! Matching observations are found using the scannumber and the station name. These two parameters can only match for true duplicates. 
    
    ! report the beginning of finalizing the ray-traced delays
    write(unit= *, fmt= '(/, a)') '======================='
    write(unit= *, fmt= '(a)') '= Finalize the delays ='
    write(unit= *, fmt= '(a, /)') '======================='
    
    
    !---------------------------------------------------------------------------------------------------
    
    ! combine and sort the observations and delays from all epologs for processed session
    
    ! report process
    write(unit= *, fmt= '(a)') 'Start combining and sorting the epolog-wise ray-traced delays ...'
    
    ! call subroutine to combine all observations, delays and meteorological parameters at the station position from the NWM contained in the
    ! epologs (in any case of epolog mode) from the whole session and sort them.
    ! The sorting is done in a way so that the observations have the order of:
    !   1.) chronological mjd
    !   2.) if same mjd then alphabetical order of station name
    !   3.) if same mjd and same station name then alphabetical order of source names
    call combine_and_sort_rd( epolog, &
                              rdlog_epoch )
    
    ! report elapsed time
    call get_elapsed_time( icount_start, count_rate, icount_interm, elapsed_time )
    write(unit= *, fmt= '(a, f0.3, a, /)') 'Total elapsed time after combining and sorting the epolog-wise ray-traced delays: ', elapsed_time, ' s'
    
    
    !---------------------------------------------------------------------------------------------------
    
    ! time domain interpolation of delays for the same observations in order to get the final delay at the exact observation time
    ! note: Time domain interpolation is only done for epolog modes where more than one delay has been calculated for one specific observation.
    
    ! report process
    write(unit= *, fmt= '(a)') 'Start calculating final delays at actual observation time through time domain interpolation of delays at specific epochs ...'
    
    ! determine the mode of assigning the observations to epochs
    select case (parameters % epolog_mode)
        ! case of one_epoch_per_obs mode
        case ('one_epoch_per_obs')
            
            ! note: No time domain interpolation necessary, because there is only one delay per observation!
        
            ! report if time interpolation is needed
            write(unit= *, fmt= '(tr4, a, a, a)') 'Due to epolog mode "', parameters % epolog_mode, '" no time domain interpolation is necessary to receive the final ray-traced delays.'
            
            ! deallocate the grib epoch vector in the combined structure as this information per observation is not needed any further
            deallocate( rdlog_epoch % grib_epoch_mjd )
    
            
        ! case of two_epochs_per_obs mode
        case ('two_epochs_per_obs')
            
            ! note: Time domain interpolation necessary, because there are two delays per observation at different grib-epochs!
            
            ! report if time interpolation is needed
            write(unit= *, fmt= '(tr4, a, a, a)') 'Due to epolog mode "', parameters % epolog_mode, '" time domain interpolation is necessary to receive the final ray-traced delays.'
            
            ! do time domain interpolation for observations which have been ray-traced for different grib-file epochs
            call time_interpolation( rdlog_epoch )
            
            ! deallocate the grib epoch vector in the combined structure as this information per observation is not needed any further
            ! note: this information is just needed for time domain interpolation using the delays of duplicate observations at different grib-epochs
            deallocate( rdlog_epoch % grib_epoch_mjd )
            
        case default ! Case if other cases fail. This should not happen as this case has already been checked in a prior step at the very begginning of the program.
            
            ! report error message
            write(unit= *, fmt= '(a)') 'Error: Unexpected epolog mode encountered! Program stopped!'
            ! stop the program
            stop
        
    end select
    
    ! report elapsed time
    call get_elapsed_time( icount_start, count_rate, icount_interm, elapsed_time )
    write(unit= *, fmt= '(a, f0.3, a, /)') 'Total elapsed time after calculating final delays at actual observation time through time domain interpolation of delays at specific epochs: ', elapsed_time, ' s'
    
    
    !---------------------------------------------------------------------------------------------------
    
    ! Output of the results
    
    ! report the beginning of finalizing the ray-traced delays
    write(unit= *, fmt= '(/, a)') '================='
    write(unit= *, fmt= '(a)') '= Create output ='
    write(unit= *, fmt= '(a, /)') '================='
    
    
    !---------------------------------------------------------------------------------------------------
    
    ! create .radiate-file as output of the ray-tracing program
    
    ! report process
    write(unit= *, fmt= '(a)') 'Start creating the .radiate-file from all final ray-tracing results for the session ...'
    
    ! call subroutine to write output file .radiate for the session
    call create_radiate_global( rdlog_epoch, &
                                parameters % save_path_radiate, &
                                parameters % RADIATE_version, &
                                parameters % RADIATE_subversion, &
                                parameters % epolog_mode, &
                                parameters % interpolation_method, &
                                parameters % raytr_method, &
                                parameters % load_filename_stat_info, &
                                epolog % epoch_name, &
                                epolog % gridres )
    
    ! report elapsed time
    call get_elapsed_time( icount_start, count_rate, icount_interm, elapsed_time )
    write(unit= *, fmt= '(a, f0.3, a, /)') 'Total elapsed time after creating the .radiate-file from all final ray-tracing results for the session: ', elapsed_time, ' s'
    
    
    !---------------------------------------------------------------------------------------------------
    
    ! create .trp-file as output of the ray tracing program
    ! note: .trp-file will only be created in case that this is specified in the parameter "create_trp"
    
    ! check if .trp-file should be created (in case "parameters % create_trp" is .TRUE.)
    if ( parameters % create_trp ) then
        
        ! report process
        write(unit= *, fmt= '(a)') 'Start creating the .trp-file for the session ...'
    
        ! determine unique station names (alphabetically sorted) and coordinates for the whole session
        ! note: Up to now each epolog contains only the station data of stations really observing at
        !       the specific epoch of the epolog. Therefore a common station dataset covering all
        !       epologs must be formed as the unique observing stations for the whole session are
        !       needed in the trp-file.
        ! attention: The final unique station data needs to be sorted alphabetically after the station name
        !            to create a sorted output in the trp-file.
        
        ! get all observing stations from all epologs of the current session
        
        ! allocate structure "observing_stations_cumul_epolog" for storing the observing stations and their
        ! coordinates and all other data contained in the "observing_stations"-structure from each epolog (with duplicates)
        ! note: the size is determined by the sum of all numbers of observing stations per each epolog
        !       sum() returns the sum of total elements in an array if no dimension is specified
        allocate( observing_stations_cumul_epolog(sum( epolog(:) % nr_observing_stations )) )
        
        ! combine observing stations from all epologs
        
        ! initialize variable for index of end for added data
        ! note: setting "ind_end" to 0 at begin of loop is necessary as it is used when determining "ind_begin"
        ind_end= 0
        
        ! loop over all epologs
        ! combine all observing stations
        do curr_epolog_ind= 1, nr_epologs
            
            ! determine beginning index for adding data
            ind_begin= ind_end + 1
            
            ! determine index of end of added data
            ! note: use the number of observing stations in the current epolog to define the end index
            ind_end= ind_end + epolog(curr_epolog_ind) % nr_observing_stations
            
            ! add all observing stations of the current epolog to the cumulating structure "observing_stations_cumul_epolog"
            observing_stations_cumul_epolog(ind_begin:ind_end)= epolog(curr_epolog_ind) % observing_stations
        
        end do
        
        ! determine the unique observing stations, their coordinates and additional data from all epologs (session unique stations)
        ! note: additional data as supported by the "observing_station_type":
        !           % nr_obs_per_obsstat... number in how many epologs the station has been present
        !           % suspend_raytr........ .TRUE. if in any epolog this value was .TRUE.
        ! call subroutine to determine the unique observing stations for the whole session
        call get_unique_stations_with_coord( observing_stations_cumul_epolog, &
                                             observing_stations_session )
        
        
        ! call subroutine to create .trp-file
        call create_trp_global( observing_stations_session, &
                                 rdlog_epoch, &
                                 parameters % save_path_trp, &
                                 parameters % RADIATE_version, &
                                 parameters % RADIATE_subversion, &
                                 parameters % epolog_mode, &
                                 parameters % interpolation_method, &
                                 parameters % raytr_method, &
                                 parameters % load_filename_stat_info, &
                                 epolog % epoch_name, &
                                 epolog % gridres )
        
        ! report elapsed time
        call get_elapsed_time( icount_start, count_rate, icount_interm, elapsed_time )
        write(unit= *, fmt= '(a, f0.3, a, /)') 'Total elapsed time after creating the .trp-file for the session: ', elapsed_time, ' s'
    
    ! in case that no trp-file should be created
    else
        
        ! report process
        write(unit= *, fmt= '(a, /)') 'Due to parameter setting, no .trp-file is created for the session.'
        
    end if
    
    
    !---------------------------------------------------------------------------------------------------
    
    ! save the errorlog data in a file
    
    ! Depending on the parameter setting a file containing information about occurred errors during
    ! the processing is created for the session.
    ! Note: In case an .err-file should be created, but no error occurred, then no file will be created.
    !       This check is done within the subroutine "create_errorlog".
        
    ! check if a file for possible occurred errors should be created
    if ( parameters % save_errorlog ) then
              
        ! report process
        write(unit= *, fmt= '(a)') 'Start creating the .err-file if necessary containing occured errors during ray-tracing of the session ...'
        
        ! call subroutine to create the file containing the errors
        ! note: within this subroutine it is checked wether an .err-file should be created at all as no file is created in case of no occurred errors
        call create_errorlog( errorlog, &
                              parameters % save_path_errorlog, &
                              session_name, &
                              parameters % RADIATE_version, &
                              parameters % RADIATE_subversion, &
                              parameters % epolog_mode, &
                              parameters % interpolation_method, &
                              parameters % raytr_method, &
                              parameters % load_filename_stat_info, &
                              epolog % epoch_name, &
                              epolog % gridres )
            
        ! report elapsed time
        call get_elapsed_time( icount_start, count_rate, icount_interm, elapsed_time )
        write(unit= *, fmt= '(a, f0.3, a, /)') 'Total elapsed time after creating the .err-file if necessary: ', elapsed_time, ' s'
        
    ! if no err-file should be created in any case
    else
        
        ! report process
        write(unit= *, fmt= '(a, /)') 'Due to parameter setting, no .err-file is created for the session.'
        
    end if
    
    
    !---------------------------------------------------------------------------------------------------
    
    ! clean up the "epolog"-structure
    ! note: Deallocation of "epolog"-structure or other variables is not necessary as they will be automatically deallocated at the end of this subroutine.
    !       In principle most substructures of the "epolog"-structure like "% observations" and "% delay" and could be deallocated directly after
    !       the creation of "rdlog_epoch", but as there is no more need of additional memory the occupied memory will not be needed until the end of the program
    !       and the forced deallocation would only lead to an additional processing step.
    
    
    !---------------------------------------------------------------------------------------------------
    
    ! End of the program is reached
    ! display information, elapsed time and current time
    
    ! display finish of processing of the session
    write(unit= *, fmt= '(/, a)') '************************************************'
    write(unit= *, fmt= '(3a)') '* Finished processing session: ', session_name, '  *'
    write(unit= *, fmt= '(a, /)') '************************************************'
    
    ! get the current date and time
    call DATE_AND_TIME(date= date_value, time= time_value)
    
    ! format the date and time
    date_formatted= date_value(1:4) // '-' // date_value(5:6) // '-' // date_value(7:8)
    time_formatted= time_value(1:2) // ':' // time_value(3:4) // ':' // time_value(5:)
    
    ! report elapsed time
    call get_elapsed_time( icount_start, count_rate, icount_interm, elapsed_time )
    write(unit= *, fmt= '(a, f0.3, a, /)') 'Total elapsed time at the end of the program: ', elapsed_time, ' s'
    
    ! write current date and time
    write(unit= *, fmt= '(a)') 'End time:'
    write(unit= *, fmt= '(2( tr3, 2a, /))') 'Date: ', date_formatted, 'Time: ', time_formatted

end subroutine RayTrace_main_global