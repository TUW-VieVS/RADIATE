! create_session_index_and_epologs_dual.f90

!****************************************************************************
!
! PURPOSE:  Subroutine to create epolog-files and session index file fo a specific azel-file (session)
!           Assignment mode to epologs: dual
!            
!           "create_session_index_and_epologs_dual" is a subroutine to create
!           epolog-files from the VLBI-observations in one azel-file.
!           Additionally a session index file is created that shows all epologs
!           necessary for one session (azel-file) when ray tracing should be done.
!           For the input azel-file epologs following the grib-file epochs will be
!           created in such a way that these created epologs contain combined each
!           observation of the specific azel-file twice.
!           The epologs contain the observations that fall into the interval of
!           the specific grib-file epoch.
!           
!           The observations are written to the "epolog" files station-wise, so that
!           each station's observations are consecutive entries.
!           
!
!           Each created "epolog" is saved as an ascii-file.
!
!           NOTE: The interval of the grib-file epochs must be set to 6 (hours) in this subroutine
!                 in order that the subroutine works correctly and the assignments of the observations
!                 to the epochs is done correctly!
! 
!
! INPUT:
!         load_path_azel............ specifies path to the loading directory of the azel-files
!         load_filename_azel........ filename of the azel-textfile for which epologs will be created
!         save_path_session_index... specifies path to the saving directory for the session index-files
!         save_path_epolog.......... specifies path to the saving directory for the epolog-files
!
!----------------------------------------------------------------------------
! 
! Fortran-file created by Armin Hofmeister
!   based on Matlab scripts created by Armin Hofmeister
! 
!----------------------------------------------------------------------------
! History:
! 
! 10.11.2014: create the Fortran-file based on the Matlab-file "create_session_index_and_epologs_dual.m"
! 11.11.2014: programming
! 12.11.2014: programming
! 13.11.2014: programming
! 17.11.2014: programming
! 18.11.2014: programming
! 19.11.2014: programming
! 20.11.2014: programming
! 25.11.2014: programming
! 26.11.2014: programming, change the dimensional architecture of the observations structure type
! 27.11.2014: add usage of interface module only for mjd2date subroutine
! 02.12.2014: correct error message, comments
! 22.01.2015: correct error that single interfaces can not be used by an only statement
! 29.01.2015: use modules for subroutines with need of explicit interface (module_mjd2date)
! 12.02.2015: use constant "len_sessname" to specify the length of a session name
!             use constant "len_gribepochname" to specify the length of a grib epoch name
! 19.02.2015: comments
! 20.04.2015: correct comments and messages
! 21.04.2015: correct comments
! 07.05.2015: enhance error message on opening AZEL-file
! 11.05.2015: correct format for writing the station and source names in the epolog with the full provided length, change to modulo()
! 03.06.2015: change to implied loop for writing the epolog and session index entries to the file
! 31.08.2015: change loop variables "first_mjd" and "last_mjd" to integer type
!             add special integer variable "curr_mjd_integer" for usage as loop index
! 09.09.2015: correct comments
! 10.09.2015: use constants for string length determination of epolog names, epolog filenames and session index filenames
! 14.09.2015: solve stack overflow problem at applying sorting order through temporal variable usage
! 15.09.2015: comments
! 17.09.2015: comments
!             add usage of constants "len_year", "len_month", "len_day", "len_hour"
! 30.09.2015: add check if all observations are assigned correctly to the epochs
! 17.12.2015: changes due to azel file extension by water vapour pressure
! 21.12.2015: add strict formatting of the field widths for the strings of year, month, day and epoch hour according to the defined constants

! Changes by Daniel Landskron:
! 19.11.2017: also enabled for individual internal creation of uniform azel-files
!
!****************************************************************************
    
    
subroutine create_session_index_and_epologs_dual( load_path_azel, &
                                                  load_filename_azel, &
                                                  save_path_session_index, &
                                                  save_path_epolog, &
                                                  readAzel, & 
                                                  load_path_filename_indAzelSpec, &
                                                  n_name_list, & 
                                                  name_list )
    
    ! Define modules to be used
    
    use module_type_definitions, only : observations_type
    
    use module_constants, only: len_sessname, len_gribepochname, len_epologname, len_session_index_filename, len_year, len_month, len_day, len_hour, deg2rad
    
    use module_char_case
    
    ! note: usage of "module_sort_type_definitions" not necessary as use declared in module for sorting and types are public
    !use module_sort_type_definitions ! module for type definitions used for sorting with module_qsort and module_msort
    
    use module_msort ! module for sorting character strings or double or integer values ascending (alphabetical)
    
    use module_date_type_definition
    
    use module_mjd2date
    
    use module_date2mjd
    
    use module_mjd2doy
    
    !----------------------------------------------------------------------------
    
    ! Variable definitions
    implicit none
    
    
    ! dummy arguments
    !----------------
    
    ! define variable for the path to the azel-file
    character(len=*), intent(in) :: load_path_azel
    
    ! define variable for the filename
    character(len=*), intent(in) :: load_filename_azel
    
    ! define variable for the saving path of the session index file
    character(len=*), intent(in) :: save_path_session_index
    
    ! define variable for the saving path of the epolog files
    character(len=*), intent(in) :: save_path_epolog
    
    ! define structure for parameter settings for the calculation process
    logical, intent(in) :: readAzel
    
    ! define variable for the azel_spec filename
    character(len=*), intent(in) :: load_path_filename_indAzelSpec
    
    ! define variable for the length of the station list
    integer, intent(in) :: n_name_list
    
    ! define variable for the list of stations
    character(len=8), dimension(n_name_list), intent(in) :: name_list
    
    
    ! local variables
    !----------------
    
    ! variable with the unit of the file
    integer, parameter :: file_unit=1
    
    ! variable for storing the opening status of the file
    integer :: open_status
    
    ! variable for storing the reading status of the file
    integer :: read_status
    
    ! variable for storing the total number of lines present in the file
    integer :: nr_lines
    
    ! variable for storing the number of comment lines present in the file
    integer :: nr_comment_lines
    
    ! variable for storing the number of lines containing data present in the file
    integer :: nr_data_lines
    
    ! variable for storing the first character of a line of the file
    character :: first_char
    
    ! define structure for temporal storage of all azel-file data
    type(observations_type), dimension(:), allocatable :: all_obs
    
    ! define structure for temporal storage of all azel-file data for sorting the data
    type(observations_type), dimension(:), allocatable :: temp_all_obs
    
    ! loop indices
    integer :: i, j, k
    
    ! data index
    integer :: ind_data
    
    ! variable for storing the session name
    ! note: see "module_constants" for length of session name
    character(len= len_sessname) :: session_name
    
    ! define structure that is used for sorting (see module_sort_type_definitions for definition)
    type(sort_char_type), dimension(:), allocatable :: sort_struct
    
    ! define structure that is used for sorting (see module_sort_type_definitions for definition)
    type(sort_char_type), dimension(:), allocatable :: sort_struct_temp
    
    ! define variable for storing the interval of the epochs (grib-file epoch)
    ! note: this interval is needed in full hour values and 24/interval must be an integer value
    !       first epoch must always be at 0 h!
    integer :: interval_epoch
    
    ! define variables for storing the field widths for the writing format of year, month, day and epoch hour
    ! note: In order to avoid a character string that is to short for string representation of the number,
    !       the string length is set to the size of the actual parameter.
    !       Later trim() will be used to get rid of the trailing blanks, so that the field width is correctly written as character without blanks.
    !       In this way complicated calculation of the needed string length is avoided by the use of an acceptable setting of a too large length.
    character(len=len_year) :: len_year_str
    character(len=len_month) :: len_month_str
    character(len=len_day) :: len_day_str
    character(len=len_hour) :: len_hour_str
    
    ! define variable for storing the epoch hours
    integer, dimension(:), allocatable :: epoch_hour
    
    ! define variable for storing the number of epochs per day
    integer :: nr_epochs_per_day
    
    ! define variable for storing the epoch names
    ! note: see "module_constants" for "len_hour"
    character(len=len_hour), dimension(:), allocatable :: epoch
    
    ! define variables for storing first and last mjd (= at at 00:00:00 h) of the observations in the azel-file
    ! note: use integer type as variables are used for loops
    integer :: first_mjd_integer, last_mjd_integer
    
    ! define variable for storing the number of days covered by the observations
    integer :: delta_days
    
    ! define variable for index for the epoch-assignment indices
    integer :: index_mjd
    
    ! define variable for storing the index of the current epoch
    !integer :: ind_epoch
    
    ! define counter for epochs
    ! note: 1st dimension: day (mjd), 2nd dimension: epoch
    integer, dimension(:, :), allocatable :: count_epoch
    
    ! define variable for storing the epolog file name for each epoch
    ! note: 1st dimension: day (mjd), 2nd dimension: epoch
    ! note: see "module_constants" for "len_epologname"
    character(len=len_epologname), dimension(:, :), allocatable :: name_epoch
    
    ! define variable for storing the indices of the observations that will be assigned to specific epochs
    ! note: 1st dimension: day (mjd), 2nd dimension: epoch, 3rd dimension: observation
    ! note: as each observation can be assigned at most once to a specific epoch, the later allocated size of the third dimension of the arrays is specified by the number of observations
    integer, dimension(:, :, :), allocatable :: asgm_epoch
    
    ! define variable for temporarily storing the content of observation assignments of a specific epoch (for assigning sorting order of epoch 0h)
    integer, dimension(:), allocatable :: asgm_epoch_temp
    
    ! define loop index variable holding the current mjd  (= at at 00:00:00 h) 
    ! note: integer type is needed to avoid warning with gfortran as it is a loop variable
    integer :: curr_mjd_integer
    
    ! define variable holding the current mjd  (= at at 00:00:00 h) 
    double precision :: curr_mjd
    
    ! define struture for storing date from conversion of mjd
    type(date_type) :: date
    
    ! define character strings for storing date values as strings
    ! note: see "module_constants" for "len_year", "len_month", "len_day"
    character(len=len_year) :: year_str
    character(len=len_month) :: month_str
    character(len=len_day) :: day_str
    character(len= len_year + len_month + len_day) :: curr_date_str
    
    ! define variable for storing the filename of the session index file
    ! note: see "module_constants" for "len_session_index_filename"
    character(len=len_session_index_filename) :: save_filename_session_index
    
    ! define variable for storing the file unit of the session index
    integer, parameter :: file_unit_si= 2
    
    ! variable for storing the number of epologs that are created for the session (azel-file)
    integer :: nr_epologs
    
    ! index variable for existing epologs
    integer :: ind_epolog
    
    ! define variable for storing the epolog names (that actually exist)
    ! note: same length as "name_epoch"
    ! note: see "module_constants" for "len_epologname"
    character(len=len_epologname), dimension(:), allocatable  :: epolog_names
    
    ! define variable for storing the grib-epoch name (epolog name without "epolog" and session name)
    ! note: see "module_constants" for length of grib epoch name
    character(len= len_gribepochname) :: grib_epoch
    
    ! define variable for storing the file unit of the epolog files (unit is always te same)
    integer, parameter :: file_unit_epolog= 3
    
    ! define format for writing the epolog textfiles
    ! attention: The format for station and source names is "a" without specifiying a field width.
    !            In this way the names are always written with the full length of each character variable as
    !            specified by the constants "len_statname" and "len_souname"
    ! notes on output format:
    ! 1x outputs one blank
    ! ss supresses plus sign output for all following formats
    character(len=*), parameter :: data_fmt='(a1, 1x, ss, i5, 1x, f11.5, 1x, i4, 1x, i3, 1x, i2, 1x, i2, 1x, f5.2, 1x, a, 1x, f20.15, 1x, f20.15, 1x, a, 1x, f6.2, 1x, f7.2, 1x, f6.2)'
    
    ! define variable for storing the file unit of the azel specifications file (unit is always the same)
    integer, parameter :: file_unit_azel_spec = 4
    
    ! variable for storing the skipped lines from azel_spec
    character(len=1000) :: skip_line
    
    ! variable for storing the azimuths from azel_spec
    integer :: i_az
    integer :: n_az
    
    ! array for all azimuths
    double precision, dimension(:), allocatable :: azimuths
    
    ! variable for storing the elevations from azel_spec as a string
    character(len=1000) :: elevations_str
    character(len=:), allocatable :: elevations_str_trim
    
    ! variable for the assigment of azimuths
    
    ! variable for the number of elevations
    double precision, dimension(:), allocatable :: elevations
    
    ! variables for the assigment of elevations
    integer :: i_el
    integer :: n_el
    
    ! variable for counting the scannumbers
    integer :: scannumber
    
    ! variables for the assigment of stations
    integer :: i_name_list
    
    ! define variables for the time information
    double precision :: mjd
    integer :: doy
    
    ! variable for storing the epoch civil date for calculating the epoch in mjd
    type(date_type) :: date_of_epoch
    
    ! string which is 'NaN'; when converted to num, one gets an NaN value
    character(len=3) :: NaN_str = 'NaN'
    
    
    ! CONSTANTS
        
    ! see "module_constants" for "len_sessname", "len_gribepochname", "len_epologname", "len_session_index_filename", "len_year", "len_month", "len_day", "len_hour"
    
    
    !----------------------------------------------------------------------------
    
    !============================================================================
    ! Body of the subroutine
    !============================================================================
    
    ! Load data from azel-file or create it
    !
    ! check, if the azel file shall be read or internally created
    if (readAzel) then
    
        ! create command window output
        print '(a, a, a, a, a/)', 'The azel file ',load_filename_azel,' is read from the directory ', load_path_azel ,'.'
        
        ! Open the file
        open(unit= file_unit, file= load_path_azel // load_filename_azel, action= 'read', status= 'old', iostat= open_status)
        ! Check if textfile can be opened (iostat == 0 --> no error)
        if (open_status == 0) then
        
            ! Read in the file once in order to determine the number of lines = number of station entries
            ! This is necessary for allocating the variables for storing the data.
        
            ! import the data using automatic detection of different entries with the correct data format
            ! The comments (e.g. header lines) are ignored.
            ! Columns:
            !  1     .... scannumber
            !  2     .... mjd
            !  3     .... year
            !  4     .... day of year
            !  5     .... hour
            !  6     .... minute
            !  7     .... sec
            !  8     .... station
            !  9     .... azimuth [rad]
            !  10    .... elevation [rad]
            !  11    .... source
            !  12    .... temperature [deg. C]
            !  13    .... pressure [hPa]
            !  14    .... water vapour pressure [hPa]
        
        
            ! Initialize the counter of total number of lines in the file
            nr_lines= 0
        
            ! Initialize the counter of number of comment lines in the file
            nr_comment_lines= 0
        
            ! loop do-if-exit
            do
                ! read next line
                read(unit= file_unit, fmt= '(a)', iostat= read_status) first_char
            
                ! test if end of file is reached and exit loop in this case
                if (is_iostat_end(read_status)) exit
            
                ! test if the first character of the line is a comment sign ("%" or "!")
                if (first_char == '%' .or. first_char == '!') then
                    nr_comment_lines= nr_comment_lines + 1
                end if
            
                ! raise counter for number of lines
                nr_lines= nr_lines + 1
            end do
        
            ! calculate number of lines containing data (number of total lines - number of comment lines)
            nr_data_lines= nr_lines - nr_comment_lines
        
            ! allocate the structure-variable that will store the read in observation data using the determined number of data lines
            allocate(all_obs(nr_data_lines))
        
            ! rewind the file
            rewind(unit= file_unit)
        
            ! initialize index for storage of the data
            ind_data= 0
        
            ! loop over all lines in the file
            do i=1, nr_lines
                ! import the azel-file data
            
                ! read next line
                read(unit= file_unit, fmt= '(a)', iostat= read_status) first_char
            
                ! test if the first character of the line is not a comment sign ("%" or "!")
                if (first_char /= '%' .and. first_char /= '!') then
                
                    ! change the file position backwards by one record (line) in order to read the line again
                    backspace(unit= file_unit)
                
                    ! raise the counter of the index for the data
                    ind_data= ind_data + 1
                
                    ! read the line again
                    ! use automatic format as recommended to reduce risk of errors
                    read(unit= file_unit, fmt=*) all_obs(ind_data) % scannr, &
                                                 all_obs(ind_data) % mjd, &
                                                 all_obs(ind_data) % year, &
                                                 all_obs(ind_data) % doy, &
                                                 all_obs(ind_data) % hour, &
                                                 all_obs(ind_data) % min, &
                                                 all_obs(ind_data) % sec, &
                                                 all_obs(ind_data) % station, &
                                                 all_obs(ind_data) % az, &
                                                 all_obs(ind_data) % elev, &
                                                 all_obs(ind_data) % source, &
                                                 all_obs(ind_data) % temp, &
                                                 all_obs(ind_data) % pres, &
                                                 all_obs(ind_data) % wvpr
            
                end if
            
            end do
        
        else
        
            ! report error message
            print '(a, a, a, a /)', 'Error: Problem with opening azel-file "', load_path_azel, load_filename_azel, '" for observation information! Program stopped!'
            ! stop the program
            stop
        
        end if
    
        ! close the file
        close(unit= file_unit)
    
    else
        
        ! create a command window output
        print '(a, a, a, a, a/)', 'The uniform azel file ',load_filename_azel,' is created internally using the specifications in ', load_path_filename_indAzelSpec , '.'
        
        ! Open the file
        open(unit= file_unit_azel_spec, file= load_path_filename_indAzelSpec, action= 'read', status= 'old', iostat= open_status)
        ! Check if textfile can be opened (iostat == 0 --> no error)
        if (open_status == 0) then
            
            ! import the data using automatic detection of different entries with the correct data format
            ! The comments (e.g. header lines) are ignored.
            ! Lines: 
            ! 1-6... comment
            ! 7..... start date (not needed here)
            ! 8..... end date (not needed here)
            ! 9..... number of epochs per day (not needed here)
            ! 10.... number of azimuths
            ! 11.... elevations list
        
            ! read first character of the next line, and if a comment, then jump over it to the next line
            do
                read(unit= file_unit_azel_spec, fmt= '(a)', iostat= read_status) first_char
                if (first_char /= '%' .and. first_char /= '!' .and. first_char /= '#') then
                    backspace(unit= file_unit_azel_spec)
                    exit
                end if
            end do
                    
            ! skip this line, as it contains mjd_start
            read(unit= file_unit_azel_spec, fmt=*) skip_line
                    
            ! read first character of the next line, and if a comment, then jump over it to the next line
            do
                read(unit= file_unit_azel_spec, fmt= '(a)', iostat= read_status) first_char
                if (first_char /= '%' .and. first_char /= '!' .and. first_char /= '#') then
                    backspace(unit= file_unit_azel_spec)
                    exit
                end if
            end do
            
            ! skip this line, as it contains mjd_end
            read(unit= file_unit_azel_spec, fmt=*) skip_line
            
            ! read first character of the next line, and if a comment, then jump over it to the next line
            do
                read(unit= file_unit_azel_spec, fmt= '(a)', iostat= read_status) first_char
                if (first_char /= '%' .and. first_char /= '!' .and. first_char /= '#') then
                    backspace(unit= file_unit_azel_spec)
                    exit
                end if
            end do
            
            ! skip this line, as it contains the number of epochs
            read(unit= file_unit_azel_spec, fmt=*) skip_line
            
            ! read first character of the next line, and if a comment, then jump over it to the next line
            do
                read(unit= file_unit_azel_spec, fmt= '(a)', iostat= read_status) first_char
                if (first_char /= '%' .and. first_char /= '!' .and. first_char /= '#') then
                    backspace(unit= file_unit_azel_spec)
                    exit
                end if
            end do
            
            ! read the number of azimuths
            read(unit= file_unit_azel_spec, fmt=*) n_az 
            
            ! read first character of the next line, and if a comment, then jump over it to the next line
            do
                read(unit= file_unit_azel_spec, fmt= '(a)', iostat= read_status) first_char
                if (first_char /= '%' .and. first_char /= '!' .and. first_char /= '#') then
                    backspace(unit= file_unit_azel_spec)
                    exit
                end if
            end do
            
            ! read the elevations as a string
            read(unit= file_unit_azel_spec, fmt='(A)') elevations_str
        
            
		    ! close the file
            close(unit= file_unit_azel_spec)
        
            ! allocate the size of azimuths and then add the values to it
            allocate(azimuths(n_az))
            azimuths(1) = 0
            do i_az= 2,n_az
                azimuths(i_az) = azimuths(i_az-1) + 360./n_az
            end do
            azimuths = azimuths*deg2rad   ! convert to radians
            
            ! allocate the size of elevations
            elevations_str_trim = trim(elevations_str)
            n_el = 1
            do i_el = 1,len(elevations_str_trim)
                if (elevations_str_trim(i_el:i_el) == ' ') then
                    n_el = n_el + 1
                 end if
            end do
                    
            allocate(elevations(1:n_el))
            read(elevations_str_trim , *) elevations
            elevations = elevations*deg2rad   ! convert to radians
            
            
            ! create the arrays of the azel-file
            
            ! allocate the structure-variable that will store the read in observation data using the determined number of data lines
            nr_data_lines = n_az*n_el*n_name_list
            allocate(all_obs(nr_data_lines))
            
            
            ! convert the remaining time information from the azel-file string to the respective numbers
            read(load_filename_azel(6:9),*) date_of_epoch % year
            read(load_filename_azel(10:11),*) date_of_epoch % month
            read(load_filename_azel(12:13),*) date_of_epoch % day
            read(load_filename_azel(14:15),*) date_of_epoch % hour
            date_of_epoch % min = 0
            date_of_epoch % sec = 0.0d0
            
            ! determine mjd from the subroutine date2mjd
            call date2mjd(date_of_epoch, mjd)
            
            ! determine doy from the subroutine mjd2doy
            call mjd2doy(mjd, doy)
            
            ! initialize index for storage of the data and scannumber
            ind_data= 0
            scannumber = 0
            do i_az = 1,n_az   ! for all azimuths
                do i_el = 1,n_el   ! for all elevations
                    
                    scannumber = scannumber+1
                    do i_name_list = 1,n_name_list   ! for all stations
                        
                      ind_data = ind_data +1
                      
                      ! create the data
                      all_obs(ind_data) % scannr = scannumber
                      all_obs(ind_data) % mjd = mjd
                      all_obs(ind_data) % year = date_of_epoch % year
                      all_obs(ind_data) % doy = doy
                      all_obs(ind_data) % hour = date_of_epoch % hour
                      all_obs(ind_data) % min = date_of_epoch % min
                      all_obs(ind_data) % sec = date_of_epoch % sec
                      all_obs(ind_data) % station = name_list(i_name_list)
                      all_obs(ind_data) % az = azimuths(i_az)
                      all_obs(ind_data) % elev = elevations(i_el)
                      all_obs(ind_data) % source = 'none    '
                      read(NaN_str,*) all_obs(ind_data)% temp
                      read(NaN_str,*) all_obs(ind_data)% pres
                      read(NaN_str,*) all_obs(ind_data)% wvpr
                      
                        
                    end do
                    
                end do
            end do
            
        else
        
            ! report error message
            print '(a, a, a /)', 'Error: Problem with opening the azel specifications file "', load_path_filename_indAzelSpec, '" for observation information! Program stopped!'
            ! stop the program
            stop
        
        end if
            
        
    end if
    
    
    !----------------------------------------------------------------------------
    
    ! get session-name from azel-filename
    ! find first "_" after "azel" and take filename until exclusive "."
    session_name=load_filename_azel(index(load_filename_azel, '_')+1:index(load_filename_azel, '.')-1)
    
    ! get number of all observations in the AZEL-file
    !nr_obs_all= nr_data_lines
    
    !----------------------------------------------------------------------------
    
    ! Sort all observations alphabetically using the station name
    
    ! --> all observations become sorted alphabetically using the station name, original chronological sorting
    ! order (received automatically through input of the azel-file, which is just chronologically sorted)
    ! is preserved for each station as mergesort algorithm is used for the sorting and this is a stable sorting
    ! algorithm
    
    ! make sure that all station-names are written in uppercase letters
    call upper(all_obs(:) % station)
    
    
    ! sort stations alphabetically
    
    ! allocate the structures used for sorting (see module_msort for help)
    allocate(sort_struct(nr_data_lines))
    allocate(sort_struct_temp((nr_data_lines+1)/2))
    
    ! assign station names and indices to the structure for sorting
    do i=1, nr_data_lines
        sort_struct(i) % value= all_obs(i) % station
        sort_struct(i) % order= i
    end do
    
    ! sort the list with the station names
    ! use mergesort for sorting (note: original order is preserved in case of equal entries; mergesort is a stable sorting algorithm)
    ! note: chronological order per stations is preserved with msort (if AZEL-file was sorted chronological)
    call msort(sort_struct, nr_data_lines, sort_struct_temp)
    
    ! sort all fields in the observations structure according to the new alphabetical order
    ! note: also the station field needs to be sorted as it has not been done yet
    ! note: use a temporal variable for the assingment of the new sorting order to avoid a stack overflow in case of a large data set
    
    ! allocate a temporal variable for assigning new sorting orders
    allocate(temp_all_obs(nr_data_lines))
    
    ! assign the data with the new order to the temporal variable
    temp_all_obs= all_obs(sort_struct % order)
    
    ! assign the sorted data to original variable
    all_obs = temp_all_obs
    
    ! deallocate the temporal variable
    deallocate(temp_all_obs)
    
    ! deacllocate the structure used for sorting
    deallocate(sort_struct)
    deallocate(sort_struct_temp)
    
    
    !----------------------------------------------------------------------------
    
    ! define field widths for the writing format of year, month, day and epoch hour,
    ! which will be needed below
    !       + 1 as "|" is printed at the end
    ! notes on output format:
    ! i0.1 outputs the minimal field width with at least one digit printed
    ! ss supresses plus sign output for all following formats
    write(unit= len_year_str, fmt= '(ss, i0.1)') (len_year)
    write(unit= len_month_str, fmt= '(ss, i0.1)') (len_month)
    write(unit= len_day_str, fmt= '(ss, i0.1)') (len_day)
    write(unit= len_hour_str, fmt= '(ss, i0.1)') (len_hour)
    
            
    !----------------------------------------------------------------------------
    
    ! Assign the observations to the corresponding grib-file epochs (all stations mixed)
    
    ! This section splits the observations from the azel-file and assigns them to the designated
    ! epochs of the grib-files.
       
    ! Assigning is defined here in the following way:
    ! An observation will be assigned to the grib-file epoch before and after the observation.
    ! This means that each observation is covered by the resulting epologs twice.
    ! Note: Also observations at the exact epoch time are covered twice, although the second epolog entry will not alter the time interpolation result!
        
    ! The assignment is done using the mjd-day and the observation hour.
    
    !    --> epologs:
    !                    epolog2_0h...... contains observations between [18h, 6h[; attention: observations before 0h have a different date than the grib-file epoch (-1 day)!
    !                    epolog2_6h...... contains observations between [0h, 12h[
    !                    epolog2_12h..... contains observations between [6h, 18h[
    !                    epolog2_18h..... contains observations between [12h, 0h[
    !    
    !                    epolog2+1_0h.... 
    
    ! define interval of the epochs in integer hours (interval of grib-files)
    ! note: this interval must be set to 6 otherwise the subroutine doesn't work correctly
    interval_epoch= 6
    
    ! check if interval of epochs is suitable to a full day
    ! note: modulo() in Fortran works like mod() in MATLAB (if non-zero result then the sign of the second argument is retained)
    if (modulo(24,interval_epoch) /= 0) then
        ! report error message
        print '(a, /)', 'Error: Problem with creating epologs: The specified interval of epochs leads to non-integer number of epochs per day! Note: First epoch is always 0h! Program stopped!'
        ! stop the program
        stop
    end if
    
    ! get the number of epochs per day
    nr_epochs_per_day= 24/interval_epoch
    
    ! allocate the variable for storing the epoch hours and names
    allocate(epoch_hour(nr_epochs_per_day), epoch(nr_epochs_per_day))
    
    ! assign epoch hours and define the names of the epochs
    ! note: first epoch is always 0h
    do i = 1, nr_epochs_per_day
        epoch_hour(i)= interval_epoch*(i-1)
        
        ! note: The field width is set according to the constant "len_hour" as represented here as a string "len_hour_str".
        !       Since "len_hour_str" has trailing blanks trim() is used.
        !       Thus the hour is converted to a string of 2 characters (fmt: ss --> no plus sign, i2.2 --> 2 leading zeros in case of shorter value)
        write (unit= epoch(i), fmt='(ss, i' // trim(len_hour_str) // '.' // trim(len_hour_str) // ')') epoch_hour(i)
    end do
    
    
    ! define first mjd-day and last mjd-day (at 00:00:00 h)
    ! note: int() is not necessary, but clarifies the assignment
    first_mjd_integer= int( floor( minval(all_obs(:) % mjd) ) )
    last_mjd_integer= int( floor( maxval(all_obs(:) % mjd) ) )
    
    ! number of days covered through observations in the azel-file
    ! note: int() is not necessary
    delta_days= last_mjd_integer - first_mjd_integer + 1 ! (+1 to get correct amount of days)
    
    ! allocate counters for epochs
    ! note: 1st dimension: day (mjd), 2nd dimension: epoch
    allocate(count_epoch(delta_days+1, nr_epochs_per_day))! add. +1 necessary for the splitted interval at epoch 0h (and actually only at epoch0): ]18h,24h[ ,[0h,6h[, if last obs. in azel-file are in the interval ]18h,24h[ --> indexing for next day
    
    ! allocate assignment variables for each epoch (1st dimension: number of days spanned by the observations + 1, 2nd dimension: number of epochs per day, 3rd dimension: max. number of assignments through number of observations)
    ! note: 1st dimension: day (mjd), 2nd dimension: epoch, 3rd dimension: observation
    allocate(asgm_epoch(delta_days+1, nr_epochs_per_day, nr_data_lines))! add. +1 necessary for the splitted interval at epoch 0h (and actually only at epoch0): ]18h,24h[ ,[0h,6h[, if last obs. in azel-file are in the interval ]18h,24h[ --> indexing for next day
    
    ! allocate variables with epoch names
    ! note: 1st dimension: day (mjd), 2nd dimension: epoch
    allocate(name_epoch(delta_days+1, nr_epochs_per_day))! add. +1 necessary for the splitted interval at epoch 0h (and actually only at epoch0): ]18h,24h[ ,[0h,6h[, if last obs. in azel-file are in the interval ]18h,24h[ --> indexing for next day
    
    ! initialize counters for epochs (at the beginning 0 for each day!)
    count_epoch= 0
    
    ! initialize assignment variables for each epoch
    asgm_epoch= 0
    
    ! initialize variables with epoch names
    name_epoch= 'none'
    
    ! define index for the epoch-assignment indices
    index_mjd=0
    
    ! Epoch name creation for all possible epochs based on the day span covered by the observations
    
    ! create epoch names for all days covered by the observations (all epochs per day will be given a name)
    epochname_loop: do curr_mjd_integer= first_mjd_integer, last_mjd_integer + 1 ! add. +1 necessary for the splitted interval at epoch 0h (and actually only at epoch0): ]18h,24h[ ,[0h,6h[, if last obs. in azel-file are in the interval ]18h,24h[ --> indexing for next day
        
        ! raise counter for current day to receive according indexing of the observations
        index_mjd= index_mjd + 1
        
        ! store current loop index value of "curr_mjd_integer" to double precision variable for usage in subroutine "mjd2date"
        ! note: For usage in subroutine "mjd2date" the input mjd must be double precision
        curr_mjd= curr_mjd_integer
        
        ! get date from *curr_mjd* in order to set corresponding filenames for the epologs
        call mjd2date(curr_mjd, date)
        
        ! convert year, month and day integer values to character strings
        
        ! convert the year
        ! note: The field width is set according to the constant "len_year" as represented here as a string "len_year_str".
        !       Since "len_year_str" has trailing blanks trim() is used.
        !       Thus the year is converted to a string of 4 characters (fmt: ss --> no plus sign, i4.4 --> 4 leading zeros in case of shorter value)
        write (unit= year_str, fmt='(ss, i' // trim(len_year_str) // '.' // trim(len_year_str) // ')') date % year
        
        ! convert the month
        ! note: The field width is set according to the constant "len_month" as represented here as a string "len_month_str".
        !       Since "len_month_str" has trailing blanks trim() is used.
        !       Thus the month is converted to a string of 2 characters (fmt: ss --> no plus sign, i2.2 --> 2 leading zeros in case of shorter value)
        write(unit= month_str, fmt='(ss, i' // trim(len_month_str) // '.' // trim(len_month_str) // ')') date % month
        
        ! convert the day
        ! note: The field width is set according to the constant "len_day" as represented here as a string "len_day_str".
        !       Since "len_day_str" has trailing blanks trim() is used.
        !       Thus the day is converted to a string of 2 characters (fmt: ss --> no plus sign, i2.2 --> 2 leading zeros in case of shorter value)
        write(unit= day_str, fmt='(ss, i' // trim(len_day_str) // '.' // trim(len_day_str) // ')') date % day
        
        ! concatenate date of current mjd
        curr_date_str= year_str // month_str // day_str
        
        ! loop over all epochs
        do i= 1, nr_epochs_per_day
            ! create epolog-name for current epoch
            name_epoch(index_mjd,i)= 'epolog2_' // session_name // '_'  // curr_date_str // epoch(i)
        end do
            
    end do epochname_loop
    
    
    !----------------------------------------------------------------------------
    
    ! Assignment of observations to the epochs
    
    ! reset index for the epoch-assignment indices
    index_mjd=0
    
    ! set index of the current epoch to 0
    !ind_epoch= 0
    
    ! assign only observations for specified day to the epochs
    curr_mjd_loop: do curr_mjd_integer= first_mjd_integer, last_mjd_integer
        
        ! raise counter for current day to receive according indexing of the observations
        index_mjd= index_mjd + 1
        
        ! store current loop index value of "curr_mjd_integer" to double precision variable for usage in later check to compare the same data types
        curr_mjd= curr_mjd_integer
        
        ! assign the observations to the corresponding epoch-intervals for the current (mjd-)day
        obs_loop: do i= 1, nr_data_lines
        
            ! only assign observations for *curr_mjd*, i.e. if day of observation is equal to *curr_mjd*
            if (floor(all_obs(i) % mjd) == curr_mjd) then
        
                ! assign data to epoch at 0h for observations in the interval [0h,6h[
                if (all_obs(i) % hour >= 0 .AND. all_obs(i) % hour < 6) then
                    
                    ! raise counter for epoch 0h
                    ! note: 1st dimension: day (mjd)= index_mjd, 2nd dimension: epoch0 = 1
                    count_epoch(index_mjd, 1)=count_epoch(index_mjd, 1)+1
                    
                    ! assign index of current observation to the assignment variable
                    ! note: 1st dimension: day (mjd)= index_mjd, 2nd dimension: epoch0 = 1, 3rd dimension: observation = i
                    asgm_epoch(index_mjd, 1, count_epoch(index_mjd, 1))= i
                        
                end if
    
                ! assign data to epoch at 6h for observations in the interval [0h,12h[
                if (all_obs(i) % hour >= 0 .AND. all_obs(i) % hour < 12) then
                        
                    ! raise counter for epoch 6h
                    ! note: 1st dimension: day (mjd)= index_mjd, 2nd dimension: epoch6 = 2
                    count_epoch(index_mjd, 2)=count_epoch(index_mjd, 2)+1
                    
                    ! assign index of current observation to the assignment variable
                    ! note: 1st dimension: day (mjd)= index_mjd, 2nd dimension: epoch6 = 2, 3rd dimension: observation = i
                    asgm_epoch(index_mjd, 2, count_epoch(index_mjd, 2))= i
                    
                end if
                
                ! assign data to epoch at 12h for observations in the interval [6h,18h[
                if (all_obs(i) % hour >= 6 .AND. all_obs(i) % hour < 18) then
                        
                    ! raise counter for epoch 12h
                    ! note: 1st dimension: day (mjd)= index_mjd, 2nd dimension: epoch12 = 3
                    count_epoch(index_mjd, 3)=count_epoch(index_mjd, 3)+1
                    
                    ! assign index of current observation to the assignment variable
                    ! note: 1st dimension: day (mjd)= index_mjd, 2nd dimension: epoch12 = 3, 3rd dimension: observation = i
                    asgm_epoch(index_mjd, 3, count_epoch(index_mjd, 3))= i
                
                end if  
                
                ! assign data to epoch at 18h for observations in the interval [12h,24h[
                if (all_obs(i) % hour >= 12 .AND. all_obs(i) % hour < 24) then
                        
                    ! raise counter for epoch 18h
                    ! note: 1st dimension: day (mjd)= index_mjd, 2nd dimension: epoch18 = 4
                    count_epoch(index_mjd, 4)=count_epoch(index_mjd, 4)+1
                    
                    ! assign index of current observation to the assignment variable
                    ! note: 1st dimension: day (mjd)= index_mjd, 2nd dimension: epoch18 = 4, 3rd dimension: observation = i
                    asgm_epoch(index_mjd, 4, count_epoch(index_mjd, 4))= i
                    
                end if
                
                ! assign data to epoch at 0h for observations in the interval [18h,0h[
                if (all_obs(i) % hour >= 18 .AND. all_obs(i) % hour < 24) then
                    
                    ! Attention: value of index_mjd must be raised by 1 at all times it is used in
                    ! order to combine all observations from
                    ! the intervals [18h,24h[ and [0h,6h[ splitted because of the datum change.
                    ! In this way all observations become gathered for epoch 0h.
                    
                    ! raise counter for epoch 0h
                    ! note: 1st dimension: day (mjd)= index_mjd, 2nd dimension: epoch0 = 1
                    count_epoch(index_mjd+1, 1)=count_epoch(index_mjd+1, 1)+1
                    
                    ! assign index of current observation to the assignment variable
                    ! note: 1st dimension: day (mjd)= index_mjd, 2nd dimension: epoch0 = 1, 3rd dimension: observation = i
                    asgm_epoch(index_mjd+1, 1, count_epoch(index_mjd+1, 1))= i
                    
                end if
            end if
        end do obs_loop
    end do curr_mjd_loop
    
    
    !----------------------------------------------------------------------------
    
    ! sort all epologs for epoch 0h alphabetically again using the station name
    ! --> only for epoch 0h necessary, because of the date switch between the first and the second half
    ! of the interval the assigned observations to the epoch 0h are not correctly sorted after the station-names
    ! (first for-loop receives station-wise sorted observations for [21h,24h[ next loop for next day appends again
    ! station-wise observations for [0h,3h[ --> mixed as two halfs)
    
    
    ! sort all epologs for 0h-epoch
        
    ! loop over all days
    do i= 1, size(count_epoch, 1) ! note: size of dimension 1 is the number of days and hence the number of epochs at 0h
        ! note: as set per definition epoch 0h is always the first epoch stored
        ! check if the epoch 0h for the current day has more than 1 entries (as 0 and 1 entry don't need to be sorted)
        if (count_epoch(i,1) > 1) then ! i... day, 1... epoch 0h
            
            ! allocate the structures used for sorting (see module_msort for help)
            allocate(sort_struct(count_epoch(i,1)))
            allocate(sort_struct_temp((count_epoch(i,1)+1)/2))
            
            ! assign station names according to the assigned observation number in the original data array using the asgm_epoch indices
            ! and create indices for sorting
            do j=1, count_epoch(i,1) ! i... day, 1... epoch 0h
                sort_struct(j) % value= all_obs(asgm_epoch(i, 1, j)) % station ! i... day, 1... epoch 0h, j... assigned original observation number
                sort_struct(j) % order= j
            end do
            
            ! sort the list with the station names
            ! use mergesort for sorting (note: original order is preserved in case of equal entries; mergesort is a stable sorting algorithm)
            ! note: chronological order per stations is preserved with msort (if AZEL-file was sorted chronological)
            call msort(sort_struct, count_epoch(i,1), sort_struct_temp)
            
            ! assigning of the sorted station names is not needed as only the sorting indices are necessary
            ! for sorting the assignment indices
            
            ! assign values to temp array for sorting the original values using the sort indices
            ! allocate the temp vector
            allocate(asgm_epoch_temp(count_epoch(i,1)))
            asgm_epoch_temp=asgm_epoch(i, 1, :)
            
            ! sort the assignment indices to the epoch using the received sorting indices
            ! note: due to fortran restrictions a direct assignment to the structure variable is not possible because of allocatable option
            !       error text: see ISO Fortran 2003 6.1.2. (C614)
            !       Therefore a loop for the assignments must be used, but in this case a temporal variable for the assignments is needed
            !       otherwise the sorting order would be destroyed by the loop.
            !       As the temporal variable is not a multidimensional array direct assignment without a loop is then possible!
            !       Alternative method: use a forall statement, with assignment a= a(sort_ind). Then the expression will use a constant copy
            !       (internal temp copy of a) of a to execute the expression and the correct results will be given just as if a= a_temp(sort_ind)
            !       would have been done in a normal do loop. No real speed gain when using forall, but it would make reading the program code quite tricky!
            do j=1, count_epoch(i,1) ! i... day, 1... epoch 0h
                asgm_epoch(i, 1, j)= asgm_epoch_temp(sort_struct(j) % order)
            end do
            
            ! deallocate temporal variable
            deallocate(asgm_epoch_temp)
            
            ! deacllocate the structure used for sorting
            deallocate(sort_struct)
            deallocate(sort_struct_temp)
            
        end if
        
    end do
    
    
    !----------------------------------------------------------------------------
    
    ! check if the assignments of the observations to the epochs were sucessful
        
    ! check if the number of assignments is not twice the number of observations since each
    ! observation should be assigned to two epochs
    ! note: sum() without dimension specification delivers sum of all elements
    if ( sum(count_epoch) /= 2 * nr_data_lines ) then
            
        ! report error message
        print '(a /)', 'Error: Not all observations are assigned correctly to the grib-file epochs! Program stopped!'
        ! stop the program
        stop
            
    end if
    
    
    !----------------------------------------------------------------------------
    
    ! write epologs according to assignments to the epochs
    
    ! determine which epolog has more than 0 entries and should therefore be output as an epolog-textfile
    ! and determine the filename for the session index file
    
    ! determine the number of epologs (epochs) that were created above by assigning the observations (more than 0 entries)
    nr_epologs= count(count_epoch > 0)
    
    ! allocate the variable for storing the epolog names
    allocate(epolog_names(nr_epologs))
    
    ! initialize index of epolog names
    ind_epolog= 0
    
    ! loop over all days
    do i= 1, size(count_epoch,1)
        ! loop over all epochs
        do j= 1, size(count_epoch,2)
            ! in case there are any entries for the specific epoch 
            if (count_epoch(i, j) > 0) then
                
                ! raise index of epolog names
                ind_epolog= ind_epolog + 1
                
                ! assign epolog name to the variable
                epolog_names(ind_epolog)= name_epoch(i, j)
                
                ! get grib-epoch name
                ! note: index(string, substring, back) where setting true for back delivers index of rightmost start of occurence of the substring
                grib_epoch= name_epoch(i, j)(index(name_epoch(i, j), '_', .TRUE.)+1:)
                
                ! create a new epolog file
                open(unit= file_unit_epolog, file= save_path_epolog // name_epoch(i, j) // '.txt', action= 'write', status= 'replace', iostat= open_status)
                
                ! Check if textfile can't be opened (iostat /= 0 --> error)
                if (open_status /= 0 ) then
                    ! report error message
                    print '(a, a, a /)', 'Error: Problem with creating epolog-textfile: ', name_epoch(i, j), '! Program stopped!'
                    ! stop the program
                    stop
                end if
                
                ! write header with information about session-name, epoch-name, total number of observations
                ! for the epoch and data description
                write(unit= file_unit_epolog, fmt= '(a)') '%******************************************************************************************************************'
                ! notes on output format:
                ! 1x outputs one blank
                ! first character in the line containing the session name is set to "S" as to faciliate the finding of the session name later when the file is read in
                write(unit= file_unit_epolog, fmt= '(a)') '% epolog2 for session:'
                write(unit= file_unit_epolog, fmt= '(a1, 1x, a)') 'S', session_name
                ! notes on output format:
                ! 1x outputs one blank
                ! first character in the line containing the grib-file epoch is set to "E" as to faciliate the finding of the epoch name later when the file is read in
                write(unit= file_unit_epolog, fmt= '(a)')'% epolog2 for epoch:'
                write(unit= file_unit_epolog, fmt= '(a1, 1x, a)')'E', grib_epoch
                ! notes on output format:
                ! 1x outputs one blank
                ! i0.1 outputs the minimal field width with at least one digit printed
                ! ss supresses plus sign output for all following formats
                ! first character in the line containing the number of observations set to "N" as to faciliate the finding of the number later when the file is read in
                write(unit= file_unit_epolog, fmt= '(a)') '% total number of observations in this epoch:'
                write(unit= file_unit_epolog, fmt= '(a1, 1x, ss, i0.1)') 'N', count_epoch(i, j)
                write(unit= file_unit_epolog, fmt= '(a)') '%'
                write(unit= file_unit_epolog, fmt= '(a)') '% Columns:'
                ! notes on output format:
                ! tr5 tabs 5 positions to the right
                write(unit= file_unit_epolog, fmt= '(a, tr4, a)') '%', '0  .... data specifier sign (S, E, N, D)'
                write(unit= file_unit_epolog, fmt= '(a, tr4, a)') '%', '1  .... scannumber'
                write(unit= file_unit_epolog, fmt= '(a, tr4, a)') '%', '2  .... mjd'
                write(unit= file_unit_epolog, fmt= '(a, tr4, a)') '%', '3  .... year'
                write(unit= file_unit_epolog, fmt= '(a, tr4, a)') '%', '4  .... day of year'
                write(unit= file_unit_epolog, fmt= '(a, tr4, a)') '%', '5  .... hour'
                write(unit= file_unit_epolog, fmt= '(a, tr4, a)') '%', '6  .... min'
                write(unit= file_unit_epolog, fmt= '(a, tr4, a)') '%', '7  .... sec'
                write(unit= file_unit_epolog, fmt= '(a, tr4, a)') '%', '8  .... station' 
                write(unit= file_unit_epolog, fmt= '(a, tr4, a)') '%', '9  .... azimuth in [rad]'
                write(unit= file_unit_epolog, fmt= '(a, tr4, a)') '%', '10 .... elevation in [rad]'
                write(unit= file_unit_epolog, fmt= '(a, tr4, a)') '%', '11 .... source'
                write(unit= file_unit_epolog, fmt= '(a, tr4, a)') '%', '12 .... temperature in [°C]'
                write(unit= file_unit_epolog, fmt= '(a, tr4, a)') '%', '13 .... pressure in [hPa]'
                write(unit= file_unit_epolog, fmt= '(a, tr4, a)') '%', '14 .... water vapour pressure in [hPa]'
                write(unit= file_unit_epolog, fmt= '(a)') '%******************************************************************************************************************'
                write(unit= file_unit_epolog, fmt= '(a)') '%'
                
                ! write the observations assigned to the current epoch to the epolog-file
                ! note: Implied loop for writing each line is used as it should be faster than a normal do loop.
                ! notes on output format:
                ! 1x outputs one blank
                ! first character in the line containing the observation data is set to "D" as to faciliate the finding of the data later when the file is read in
                write(unit= file_unit_epolog, fmt= data_fmt) ( 'D', &
                                                               all_obs(asgm_epoch(i, j, k)) % scannr, &
                                                               all_obs(asgm_epoch(i, j, k)) % mjd, &
                                                               all_obs(asgm_epoch(i, j, k)) % year, &
                                                               all_obs(asgm_epoch(i, j, k)) % doy, &
                                                               all_obs(asgm_epoch(i, j, k)) % hour, &
                                                               all_obs(asgm_epoch(i, j, k)) % min, &
                                                               all_obs(asgm_epoch(i, j, k)) % sec, &
                                                               all_obs(asgm_epoch(i, j, k)) % station, &
                                                               all_obs(asgm_epoch(i, j, k)) % az, &
                                                               all_obs(asgm_epoch(i, j, k)) % elev, &
                                                               all_obs(asgm_epoch(i, j, k)) % source, &
                                                               all_obs(asgm_epoch(i, j, k)) % temp, &
                                                               all_obs(asgm_epoch(i, j, k)) % pres, &
                                                               all_obs(asgm_epoch(i, j, k)) % wvpr, k= 1, count_epoch(i, j) )
                
                ! close the created textfile
                close(unit= file_unit_epolog)
                
            end if
        end do
    end do
    
    
    !----------------------------------------------------------------------------
    
    ! Create session index file
    ! The session index displays all epologs necessary for processing ray traced delays for the
    ! input session (azel-file).
    
    ! define filename for session index file
    save_filename_session_index= 'session_index2_' // session_name // '.txt'
    
    ! create new textfile for storing session index
    open(unit= file_unit_si, file= save_path_session_index // save_filename_session_index, action= 'write', status= 'replace', iostat= open_status)
    
    ! Check if textfile can't be opened (iostat /= 0 --> error)
    if (open_status /= 0) then
        ! report error message
        print '(a, a, a /)', 'Error: Problem with creating session index-textfile: ', save_filename_session_index, '! Program stopped!'
        ! stop the program
        stop
    end if
    
    ! write header with information about session-name, epoch-name, total number of observations
    ! for the epoch and data description
    write(unit= file_unit_si, fmt= '(a)') '%*****************************************************************************************************************'
    ! notes on output format:
    ! first character in the line containing the session name is set to "S" as to faciliate the finding of the session name later when the file is read in
    ! 1x outputs one blank
    write(unit= file_unit_si, fmt= '(a)') '% session index2 for session:'
    write(unit= file_unit_si, fmt= '(a, 1x, a)') 'S', session_name
    ! notes on output format:
    ! 1x outputs one blank
    ! i0.1 outputs the minimal field width with at least one digit printed
    ! ss supresses plus sign output for all following formats
    ! first character in the line containing the number of observations set to "N" as to faciliate the finding of the number later when the file is read in
    write(unit= file_unit_si, fmt= '(a)') '% total number of epologs2 necessary for processing this session:'
    write(unit= file_unit_si, fmt= '(a, 1x, ss, i0.1)') 'N', nr_epologs
    write(unit= file_unit_si, fmt= '(a)') '%*****************************************************************************************************************'
    write(unit= file_unit_si, fmt= '(a, a)') '%'
    
    ! write the epolog (epoch) names with specified extension (".txt") to the session index file
    ! note: Implied loop for writing each line is used as it should be faster than a normal do loop.
    write(unit= file_unit_si, fmt= '(a1, 1x, a)') ( 'D', epolog_names(i) // '.txt', i= 1, nr_epologs )
    
    ! close the created textfile
    close(unit= file_unit_si)
    
    
end subroutine create_session_index_and_epologs_dual