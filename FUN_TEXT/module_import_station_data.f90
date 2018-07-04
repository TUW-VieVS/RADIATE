! module_import_station_data.f90

!****************************************************************************
!
! PURPOSE:  Module containing subroutine to import station dependent data
!            
!           "import_station_data" loads specific data for all available stations from a textfile
!           specified in the input.
!
!           Furthermore a vector with station names of all stations is created: "name_list",
!           vectors containing latitude (ellipsoidal), longitude (ellipsoidal) and height (ellipsoidal)
!           are created.
!
!           Attention: In principle the interval for the latitude values is assumed to be [-90°, 90°] and
!                      for the longitude values [0°, 360°[.
!                      In case of encountered negative longitude values they are assumed to be in the interval [-180°, 180°]
!                      and to these negative values 360° is added to convert them to the interval [0°, 360°[.
!
!           All output variables are sorted alphabetically using the station names!
! 
!
! INPUT:
!         load_filename...   specifies the path and name of the textfile ("vlbi.ell" or files with equal layout)
! 
! 
! OUTPUT:
!         name_list......   vector containing only the station names sorted alphabetically
!         lat_ell........   vector containing ellipsoidal latitude values in [°] of all stations in the
!                           same order as in the variable "name_list"
!         lon_ell........   vector containing ellipsoidal longitude values in [°] of all stations in the
!                           same order as in the variable "name_list"
!         h_ell..........   vector containing ellipsoidal height values in [m] of all stations in the
!                           same order as in the variable "name_list"
!
!----------------------------------------------------------------------------
! 
! Fortran-file created by Armin Hofmeister
!   based on Matlab scripts created by Armin Hofmeister
! 
!----------------------------------------------------------------------------
! History:
! 
! 04.11.2014: create the Fortran-file based on the Matlab-file "import_stations_textscan.m"
! 05.11.2014: programming
! 06.11.2014: programming
! 10.11.2014: comments
! 11.11.2014: comments
! 12.11.2014: programming
! 13.11.2014: comments
! 27.11.2014: comments
! 29.01.2015: change to module as to save explicit interface block
! 12.02.2015: programming
!             use constant to define length of station names
! 20.04.2015: correct comments and messages
! 21.04.2015: correct comments
! 07.05.2015: enhance error message on opening station data file
! 11.05.2015: add comments
! 03.06.2015: add comments
! 09.06.2015: add check of longitude values if the fit to interval [0°, 360°[
!             and adapt if necessary
! 09.09.2015: add check for comment lines in the station file
! 14.09.2015: solve stack overflow problem at applying sorting order through temporal variable usage
!
!****************************************************************************

module module_import_station_data
    
contains
    
    subroutine import_station_data( load_path, &
                                    load_filename, &
                                    name_list, &
                                    lat_ell, &
                                    lon_ell, &
                                    h_ell )

        ! Define modules to be used
        ! note: usage of "module_sort_type_definitions" not necessary as use declared in module for sorting and types are public
        !use module_sort_type_definitions ! module for type definitions used in sorting subroutines
        use module_qsort ! module for sorting character strings or double values ascending (alphabetical)
        
        use module_constants, only: len_statname
    
        !----------------------------------------------------------------------------
    
        ! Variable definitions
        implicit none
    
    
        ! dummy arguments
        !----------------
    
        ! INPUT
        
        ! define variable for the path to the file
        character(len=*), intent(in) :: load_path
        
        ! define variable for the filename (possibly including the path to the file)
        character(len=*), intent(in) :: load_filename
        
        ! OUTPUT
        
        ! define variable for storing the station names
        ! note: see "module_constants" for length of station names
        character(len= len_statname), dimension(:), allocatable, intent(out) :: name_list ! attention: set length of each string to 8 characters, dimension will be specified through counting the number of lines in the file
    
        ! define variable for storing the ellipsoidal latitude values in [°]
        double precision, dimension(:), allocatable, intent(out) :: lat_ell ! attention: dimension will be specified through counting the number of lines in the file
    
        ! define variable for storing the ellipsoidal longitude values in [°]
        double precision, dimension(:), allocatable, intent(out) :: lon_ell ! attention: dimension will be specified through counting the number of lines in the file
    
        ! define variable for storing the ellipsoidal height values in [m]
        double precision, dimension(:), allocatable, intent(out) :: h_ell ! attention: dimension will be specified through counting the number of lines in the file
    
    
        ! local variables
        !----------------
    
        ! variable with the unit of the file
        integer, parameter :: file_unit=1
    
        ! variable for storing the opening status of the file
        integer :: open_status
    
        ! variable for storing the reading status of the file
        integer :: read_status
    
        ! variable for storing the number of lines present in the file
        integer :: nr_lines
        
        ! variable for storing the number of comment lines present in the file
        integer :: nr_comment_lines
    
        ! variable for storing the number of lines containing data present in the file
        integer :: nr_data_lines
    
        ! variable for storing the first character of a line of the file
        character :: first_char
    
        ! data index
        integer :: ind_data
        
        ! loop index
        integer :: i
    
        ! define structure that is used for sorting (see module_sort_type_definitions for definition)
        type(sort_char_type), dimension(:), allocatable :: sort_struct
        
        ! define temporal variable for storing the ellipsoidal latitude values in [°]
        ! note: The temporal variable is needed to avoid a stack overflow during assigning the new sorting order in case of a large station catalogue
        double precision, dimension(:), allocatable :: temp_lat_ell ! attention: dimension will be specified through counting the number of lines in the file
    
        ! define temporal variable for storing the ellipsoidal longitude values in [°]
        ! note: The temporal variable is needed to avoid a stack overflow during assigning the new sorting order in case of a large station catalogue
        double precision, dimension(:), allocatable :: temp_lon_ell ! attention: dimension will be specified through counting the number of lines in the file
    
        ! define temporal variable for storing the ellipsoidal height values in [m]
        ! note: The temporal variable is needed to avoid a stack overflow during assigning the new sorting order in case of a large station catalogue
        double precision, dimension(:), allocatable :: temp_h_ell ! attention: dimension will be specified through counting the number of lines in the file
    
        
        ! CONSTANTS
        
        ! see "module_constants" for "len_statname"
        
    
        !----------------------------------------------------------------------------
    
        !============================================================================
        ! Body of the subroutine
        !============================================================================
    
        ! Load data for VLBI stations from file
        ! 
        ! Open and read the textfile "load_filename" which contains names and ellipsoidal coordinates for
        ! the stations.
    
        ! Open the file
        open(unit= file_unit, file= load_path // load_filename, action= 'read', status= 'old', iostat= open_status)
        ! Check if textfile can be opened (iostat == 0 --> no error)
        if (open_status == 0) then
        
            ! Read in the file once in order to determine the number of lines = number of station entries
            ! This is necessary for allocating the variables for storing the data.
        
            ! Initialize the counter of number of lines in the file
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
            
            ! allocate the variables that will be read in using the determined number of data lines
            allocate(name_list(nr_data_lines), lat_ell(nr_data_lines), lon_ell(nr_data_lines), h_ell(nr_data_lines))
        
            ! rewind the file
            rewind(unit= file_unit)
            
            ! initialize index for storage of the data
            ind_data= 0
        
            ! loop over all lines in the file
            do i= 1, nr_lines
                
                ! read next line
                read(unit= file_unit, fmt= '(a)', iostat= read_status) first_char
                
                ! test if the first character of the line is not a comment sign ("%" or "!")
                if (first_char /= '%' .and. first_char /= '!') then
                    
                    ! change the file position backwards by one record (line) in order to read the line again
                    backspace(unit= file_unit)
                    
                    ! raise the counter of the index for the data
                    ind_data= ind_data + 1
                    
                    ! read the line again to import the station data
                    ! use automatic format as recommended to reduce risk of errors
                    ! read in the current line using automatic format as recommended to reduce risk of errors
                    read(unit= file_unit, fmt=*) name_list(ind_data), lat_ell(ind_data), lon_ell(ind_data), h_ell(ind_data)
                    
                end if
                
            end do
        
        else
            ! report error message
            print '(a, a, a, a /)', 'Error: Problem with opening textfile "', load_path, load_filename, &
                  '" for station information! Program stopped!'
            ! stop the program
            stop
        
        end if
    
        ! close the file
        close(unit= file_unit)
    
        !----------------------------------------------------------------------------
        
        ! Check if the longitude of the station position is given in the correct interval
        ! This section checks if the input longitudes of the station positions are in the interval [0°, 360°[.
        ! If a negative longitude is found, it is assumed that the interval [-180°, 180°] is used for this station
        ! and the value is corrected to the needed interval [0°, 360°[ by simply adding 360°. All positive longitudes
        ! should be correct no matter which interval, so all positive longitudes always apply to the interval [0°, 360°[.
        
        ! check all longitude values for negative values assuming then they are in the interval [-180°, 180°]
        where (lon_ell < 0)
            ! add 360° to get the longitude to interval [0°, 360°[
            lon_ell= lon_ell + 360
        end where
    
        !----------------------------------------------------------------------------
    
        ! Sort the station infos alphabetically using the station names
    
        ! allocate the structure used for sorting
        allocate(sort_struct(nr_data_lines))
    
        ! assign station names and indices to the structure for sorting
        do i=1, nr_data_lines
            sort_struct(i) % value= name_list(i)
            sort_struct(i) % order= i
        end do
    
        ! sort the list with the station names "name_list"
        ! use quicksort for sorting (note: original order is not preserved in case of equal entries; quicksort is not a stable sorting algorithm)
        call qsort(sort_struct, nr_data_lines)
    
        ! assign the sorted station names (note: loop is necessary because of allocatable type of structure)
        do i=1, nr_data_lines
            name_list(i)=sort_struct(i) % value
        end do
    
        ! allocate the temporal variables for assigning new sorting orders
        allocate( temp_lat_ell(nr_data_lines), &
                  temp_lon_ell(nr_data_lines), &
                  temp_h_ell(nr_data_lines) )
        
        ! assign the data with the new order to the temporal variables
        temp_lat_ell= lat_ell(sort_struct % order)
        temp_lon_ell= lon_ell(sort_struct % order)
        temp_h_ell= h_ell(sort_struct % order)
        
        ! assign the alphabetically sorted data to original variables
        lat_ell= temp_lat_ell
        lon_ell= temp_lon_ell
        h_ell= temp_h_ell
        
        ! deallocate the temporal variables
        deallocate( temp_lat_ell, &
                    temp_lon_ell, &
                    temp_h_ell )
    
        ! deacllocate the structure used for sorting
        deallocate(sort_struct)


    end subroutine import_station_data
    
end module module_import_station_data