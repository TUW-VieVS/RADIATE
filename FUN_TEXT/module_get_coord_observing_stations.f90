! module_get_coord_observing_stations.f90

!****************************************************************************
!
! PURPOSE:  Module containing subroutine to import the filenames of all epolog-files for the specific session created with "create_session_index_and_epologs" (dual or single)
! 
!           "get_coord_observing_stations" is a subroutine to determine the latitude, longitude and height of an
!           observing station by searching in data from a catalogue file using the station's name.
!           If more than one matching entry or no match has been found an appropriate warning message is
!           created and the encountered error is reported to an error log.
!           For not found stations ray-tracing will be suspended. Therefore a marker is set.
! 
! 
! INPUT:
!         nr_observing_stations... number of observing stations (unique)
!         obs_names............... vector of names of all observing stations, that should be found in the catalogue
!         nr_cat_names............ number of catalogue stations
!         cat_names............... vector of names of all stations in the catalogue
!         cat_lat_ell............. vector of (ellipsoidal) latitude in [°] according to the catalogue stations
!         cat_lon_ell............. vector of (ellipsoidal) longitude in [°] according to the catalogue stations
!         cat_h_ell............... vector of (ellipsoidal) height in [m] according to the catalogue stations
! 
! 
! OUTPUT:
!         obs_lat_ell............. vector of (ellipsoidal) latitude in [°] according to the observing stations
!         obs_lon_ell............. vector of (ellipsoidal) longitude in [°] according to the observing stations
!         obs_h_ell............... vector of (ellipsoidal) height in [m] according to the observing stations
!         errorlog................ structure to store possible errors during the processing
!         suspend_raytr........... marker to decide if ray-tracing will be done
!                                  0... ray-tracing will be done for this specific station
!                                  1... ray-tracing will be suspended for this specific station
!
!----------------------------------------------------------------------------
! 
! Fortran-file created by Armin Hofmeister
!   based on Matlab scripts created by Armin Hofmeister
! 
!----------------------------------------------------------------------------
! History:
! 
! 27.11.2014: create the Fortran-file based on the Matlab-file "get_coordinates_of_observing_stations.m"
! 01.12.2014: programming
! 02.12.2014: programming
! 03.12.2014: comments
! 20.01.2015: comments, enhance warning message for missing station
! 29.01.2015: change to module as to save explicit interface block
!             rename mdule and subroutine
! 12.02.2015: use assumed length operator for input of station names
! 05.05.2015: programming: declare dimension of input variables for catalogue station position data depending on size of catalogue station names, comments
! 06.05.2015: improve the code by using pack() instead of a do loop to determine the indices of matching catalog data with searched station name
! 11.05.2015: add comments
! 14.09.2015: enhance program code by using nr_observing_stations and nr_cat_names as input to determine dimensions of variables
!             in the declaration part
!
!****************************************************************************

module module_get_coord_observing_stations

contains
    
    subroutine get_coord_observing_stations( nr_observing_stations, &
                                             obs_names , &
                                             nr_cat_names, &
                                             cat_names, &
                                             cat_lat_ell, &
                                             cat_lon_ell, &
                                             cat_h_ell, &
                                             obs_lat_ell, &
                                             obs_lon_ell, &
                                             obs_h_ell, &
                                             errorlog, &
                                             suspend_raytr)

        ! Define modules to be used
    
        use, intrinsic :: ieee_arithmetic ! intrinsic module to e.g. specify NaN and Inf
    
        use module_type_definitions, only: errorlog_type
        
        use module_constants, only: len_statname
    
    
        !----------------------------------------------------------------------------
    
        ! Variable definitions
        implicit none
    
    
        ! dummy arguments
        !----------------
    
        ! INPUT
        
        ! define variable for the input of the number of observing stations (unique)
        integer, intent(in) :: nr_observing_stations
        
        ! define variable for the input of (unique) observing station names
        character(len=*), dimension(:), intent(in) :: obs_names
    
        ! define variable for the input of the number of catalogue stations
        integer, intent(in) :: nr_cat_names
    
        ! define variable for input of station names from the catalogue
        character(len=*), dimension(:), intent(in) :: cat_names
    
        ! vector for input of the latitudes from the station catalogue
        ! note: size is the number of input catalogue stations
        double precision, dimension(nr_cat_names), intent(in) :: cat_lat_ell
    
        ! vector for input of the longitudes from the station catalogue
        ! note: size is the number of input catalogue stations
        double precision, dimension(nr_cat_names), intent(in) :: cat_lon_ell
    
        ! vector for input of the heights from the station catalogue
        ! note: size is the number of input catalogue stations
        double precision, dimension(nr_cat_names), intent(in) :: cat_h_ell
    
        
        ! OUTPUT
    
        ! vector for storing the latitudes of the found observing stations from the catalogue
        ! note: size is the number of observing stations
        double precision, dimension(nr_observing_stations), intent(out) :: obs_lat_ell
    
        ! vector for storing the longitudes of the found observing stations from the catalogue
        ! note: size is the number of observing stations
        double precision, dimension(nr_observing_stations), intent(out) :: obs_lon_ell
    
        ! vector for storing the heights of the found observing stations from the catalogue
        ! note: size is the number of observing stations
        double precision, dimension(nr_observing_stations), intent(out) :: obs_h_ell
    
    
        ! define variable for output of occurred errors
        type(errorlog_type), dimension (:), allocatable, intent(out) :: errorlog
    
    
        ! define variable that reports if ray-tracing for a specific observing station needs to be suspended (default: .FALSE. --> no suspension == ray-tracing will be done)
        logical, dimension(nr_observing_stations), intent(out) :: suspend_raytr
    
    
        ! local variables
        !----------------
    
        ! define logical vector for checking if observing station has been found in the catalog
        ! note: size is the number of observing stations
        logical, dimension(nr_observing_stations) :: station_found
    
        ! define variable for storing the "NaN" value
        double precision :: NaN
    
        ! loop variable
        integer :: i
        
        ! define index vector for the indices of the catalogue data
        ! note: size is the number of station names in the catalogue
        integer, dimension(nr_cat_names) :: ind_vec
    
        ! define logical vector for storing comparison results of station names
        ! note: size is the number of station names in the catalogue
        logical, dimension(nr_cat_names) :: station_log_ind
    
        ! define variable for storing the number of found matching station names
        integer :: nr_match
    
        ! define variable for temporarily storing occurred errors
        ! note: size is the number of observing stations as the maximum error number possible is the same as the number of observing stations
        type(errorlog_type), dimension(nr_observing_stations) :: temp_errorlog
    
        ! define variable for storing the index of the current error
        integer :: ind_error
    
        ! define variables for storing error messages
        ! note: variables are allocatable to fit for the specific message
        character(len=:), allocatable :: error_description_duplicate_station_entry
        character(len=:), allocatable :: error_description_station_not_found
    
        ! define variable for total number of found observing stations in the catalogue
        integer :: nr_found_stations
    
        ! define variable for storing the (first) index of the station in the catalogue which matches the observing station
        integer, dimension(:), allocatable :: ind_match
    
        
        ! CONSTANTS
        
        ! see "module_constants" for "len_statname"
        
    
        !----------------------------------------------------------------------------
    
        !============================================================================
        ! Body of the subroutine
        !============================================================================
        

        ! initialize indicator variable for checking if station has been found in catalog
        ! note: size is determined in definition by number of observing stations
        station_found= .FALSE.
    
        ! determine the NaN-value as signaling NaN (producing exceptions when used in calculation)
        NaN= ieee_value(0.0d0, ieee_signaling_nan) ! 0.0d0 specifies double precision kind
    
        ! initialize variables for observing stations to which matching catalogue data will be
        ! assigned as NaN
        obs_lat_ell= NaN ! ellips. latitude of station in [°], interval [-90°,90°]
        obs_lon_ell= NaN ! ellips. longitude of station in [°], interval [0°,360°[
        obs_h_ell= NaN ! ellips. height of station in [m]
    
        ! initialize variable for suspension of ray-tracing to default of .FALSE.
        suspend_raytr= .FALSE.
        
        ! create an index vector for the indices of the catalogue data
        ! note: size is the number of station names in the catalogue
        ind_vec= [ (i, i= 1, nr_cat_names) ]
    
        ! initialize index variable for errors
        ind_error= 0
    
        ! find the according station coordinates and assign them to the variables
        do i= 1, nr_observing_stations
    
            ! compare name of one of the observing stations with all names in the station catalogue
            station_log_ind= obs_names(i) == cat_names
        
            ! determine number of matches found
            ! note: this should normally be 1 or 0
            !       In case the value is more than 1 the catalogue file should be revised!
            ! note: count determines the number of true elements
            nr_match= count(station_log_ind)
        
        
            ! in case one or more matches have been found
            if (nr_match > 0) then
                
                ! determine the index (or indices in case of duplicate entries in the station catalogue)
                ! of the found station match
                
                ! allocate the index variable for the matching station
                ! note: size is determined by the number of found matches
                allocate( ind_match(nr_match) )
                
                ! assign the index or indices in case of duplicates to the variable
                ind_match= pack(ind_vec, station_log_ind)
                
                ! assign the values of the first matching station in the catalogue to the observing station
                ! note: (in case there is more than one match) only the first match is considered to get the station data from the catalogue
            
                ! assign values of found match
                obs_lat_ell(i)= cat_lat_ell(ind_match(1)) ! ellips. latitude of station in [°], interval [-90°,90°]
                obs_lon_ell(i)= cat_lon_ell(ind_match(1)) ! ellips. longitude of station in [°], interval [0°,360°[
                obs_h_ell(i)= cat_h_ell(ind_match(1)) ! ellips. height of station in [m]
            
                ! assign an indicator in order to see which station has been found in the station catalog
                station_found(i)= .TRUE.
                
                ! deallocate the index variable for the matching station
                deallocate( ind_match )
                    
                ! check for and treat possible error in case there is more than one matching station in the catalogue
                ! attention: if variable station_log_ind has more than one entry that is .TRUE., this means that the station catalogue contains the same station twice    
                if (nr_match > 1) then
            
                    ! raise index variable for errors
                    ind_error= ind_error + 1
            
                    ! define error description
                    error_description_duplicate_station_entry= 'For station "' // obs_names(i) // '" more than one matching station has been found in the catalog! Therefore the first matching entry is taken for receiving the station coordinates!'
            
                    ! report message
                    write(unit= *, fmt= '(tr4, a)') error_description_duplicate_station_entry
            
                    ! assign error to the "errorlog" variable
                    temp_errorlog(ind_error) % error_nr = ind_error
                    temp_errorlog(ind_error) % error_type= 'Duplicate station entry in catalogue'
                    temp_errorlog(ind_error) % error_descr= error_description_duplicate_station_entry
            
                    ! deallocate the current error description
                    deallocate(error_description_duplicate_station_entry)
            
                end if
            
            ! if no match has been found
            else if (nr_match == 0) then ! note: in principle == 0 checking is not necessary, just an else would be sufficient
            
                ! raise index variable for errors
                ind_error= ind_error + 1
            
                ! define error description
                error_description_station_not_found= 'Station "' // obs_names(i) // '" has not been found in the catalog! Please add station data to the file! Ray tracing for this station will be suspended!'
            
                ! report message
                write(unit= *, fmt= '(tr4, a, a)') 'Warning: ', error_description_station_not_found
            
                ! assign error to the "errorlog" variable
                temp_errorlog(ind_error) % error_nr = ind_error
                temp_errorlog(ind_error) % error_type= 'Station missing in catalogue'
                temp_errorlog(ind_error) % error_descr= error_description_station_not_found
            
                ! deallocate the current error description
                deallocate(error_description_station_not_found)
            
                ! set station data for missing station to NaN (not necessary here, already done via
                ! initializing storage variables as NaN)
                ! --> due to invalid coordinates of station --> no meteorological data will be
                ! extracted-->  no raytracing will be done
            
                ! set variable for suspending ray-tracing to .TRUE.
                ! note: initializiation to .FALSE. should have been done in the beginning
                suspend_raytr(i)= .TRUE.
            
            end if
        
        end do
    
        ! determine number of found stations in the catalog
        ! note: count determines the number of true elements
        nr_found_stations=count(station_found)
    
        ! determine if all stations have been found
        if (nr_observing_stations == nr_found_stations) then
             write(unit= *, fmt= '(tr4, a)') 'All station data has been found in the catalog for all observing stations in the current epolog!'
        end if
    
        ! allocate output structure of errorlog with necessary size as defined by the number of errors (ind_error value)
        allocate(errorlog(ind_error))
    
        ! assign errors from temporary structure to output structure
        errorlog= temp_errorlog(1:ind_error)
    
    end subroutine get_coord_observing_stations
    
end module module_get_coord_observing_stations