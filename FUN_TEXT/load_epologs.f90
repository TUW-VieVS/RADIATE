! load_epologs.f90

!****************************************************************************
!
! PURPOSE:  Subroutine to import the epolog data created with "create_epologs".
!            
!           The output "log" is a structure containing all data contained in each epoch's epolog.
! 
!
! INPUT:
!         load_path.......   specifies path to the loading directory
!         load_filename...   filename of the epolog
!         delete_file..... variable for determining if file should be deleted at the closing:
!                          .TRUE. .... delete the file
!                          .FALSE. ... keep the file
! 
! 
! OUTPUT:
!         epolog............  structure array containing the epolog-data
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
! 25.11.2014: create the Fortran-file based on the Matlab-file "load_epologs.m"
! 26.11.2014: programming, changes due to the new dimensional architecture of the observations structure type
! 08.01.2015: correct comments, error message
! 12.01.2015: correct comment
! 12.02.2015: use constant "len_sessname" to specify the length of a session name
!             use constant "len_gribepochname" to specify the length of a grib epoch name
! 20.04.2015: correct comments and messages
! 07.05.2015: add possibility to delete the loaded epolog-file after the loading when it is closed
! 11.05.2015: add comments
! 17.12.2015: changes due to azel file extension by water vapour pressure
!
!****************************************************************************
    
    
subroutine load_epologs( load_path, &
                         load_filename, &
                         delete_file, &
                         epolog )

    ! Define modules to be used
    
    use module_type_definitions, only : epolog_type
    
    use module_constants, only: len_sessname, len_gribepochname
    
    !----------------------------------------------------------------------------
    
    ! Variable definitions
    implicit none
    
    
    ! dummy arguments
    !----------------
    
    ! INPUT
    
    ! define variable for the path to the file
    character(len=*), intent(in) :: load_path
    
    ! define variable for the filename
    character(len=*), intent(in) :: load_filename
    
    ! define variable for determining if file should be deleted at the closing
    logical, intent(in) :: delete_file
    
    
    ! OUTPUT
    
    ! define structure for storing the loaded epolog data (observations) that will be the output
    type(epolog_type), intent(out) :: epolog
    
    
    ! local variables
    !----------------
    
    ! variable with the unit of the file
    integer, parameter :: file_unit=1
    
    ! variable for storing the opening status of the file
    integer :: open_status
    
    ! variable for storing the reading status of the file
    integer :: read_status
    
    ! define variable for storing the number of data entries (observations) in the file
    integer :: nr_data
    
    ! variable for storing the first character of a line of the file
    character :: first_char
    
    ! define variables for checking if the "S", "E" and "N" entry has already been found in the file
    logical :: S_found, E_found, N_found
    
    ! define variable for storing the session name
    ! note: see "module_constants" for length of session name
    character(len= len_sessname) :: session_name
    
    ! define variable for storing the session name
    ! note: see "module_constants" for length of grib epoch name
    character(len= len_gribepochname) :: epoch_name
    
    ! index for data entries
    integer :: i
    
    
    ! CONSTANTS
        
    ! see "module_constants" for "len_sessname" and "len_gribepochname"
    
    
    !----------------------------------------------------------------------------
    
    !============================================================================
    ! Body of the subroutine
    !============================================================================
    
    
    ! Load data from epolog-file
    
    ! Columns in the epolog file:
    !    0     .... data specifier sign: S= session name, E= epoch of grib-file, N= number of data entries (observations), D= data entry (observation)
    !	 1     .... scannumber
    !	 2     .... mjd
    !	 3     .... year
    !	 4     .... day of year
    !    5     .... hour
    !    6     .... minute
    !    7     .... sec
    !	 8     .... station
    !	 9     .... azimuth [rad]
    !	 10    .... elevation [rad]
    !	 11    .... source
    !	 12    .... temperature [deg. C]
    !	 13    .... pressure [hPa]
    !    14    .... water vapour pressure [hPa]
    
    
    ! Open and read the epolog-textfile "filename" which contains observations of a specific session assigned to the specific grib-file epoch
    open(unit= file_unit, file= load_path // load_filename, action= 'read', status= 'old', iostat= open_status)
    
    ! Check if textfile can be opened (iostat == 0 --> no error)
    if (open_status == 0) then
        
        ! Read in the file once in order to determine the number of data (observation) entries
        ! This is necessary for allocating the variable for storing the filenames.
        
        ! initialize the variables for checking if the "S", "E" and "N" entry has already been found in the file
        S_found= .FALSE.
        E_found= .FALSE.
        N_found= .FALSE.
        
        ! initialize the variable of number of data entries (observations) in the file
        nr_data= 0
        
        ! loop do-if-exit over the file
        do
            
            ! check if all needed data specifications have been found (except for the real data entries signalled by "D")
            ! then exit the loop as task of finding the data specifications has been acomplished
            ! note: the check must be placed here as cyle-statements are use in the following checks
            if ( S_found .AND. E_found .AND. N_found ) exit
            
            
            ! read next line
            read(unit= file_unit, fmt= '(a)', iostat= read_status) first_char
            
            ! test if end of file is reached and exit loop in this case
            if (is_iostat_end(read_status)) exit
            
            ! check if current line is a comment line (with "%" or "!" as first sign)
            ! in this case cycle one loop step
            if (first_char == '%' .OR. first_char == '!') cycle
            
            ! check if first character is a "S" signalling the session name
            ! note: extra check for found flag of specifier may not be necessary, but if done so only the first appearance of the specifier will be used for reading in
            if ( (first_char == 'S') .AND. (.NOT. S_found)) then
                ! change the file position backwards by one record (line) in order to read the line again
                backspace(unit= file_unit)
                
                ! read in the line containing the number of data entries using using automatic format as recommended to reduce risk of errors
                read(unit= file_unit, fmt= *, iostat= read_status) first_char, session_name
                
                ! set logical variable reporting the finding of the data specifier to .true.
                S_found= .TRUE.
                
                ! cycle to next loop step as value has been found and no further statements in this loop step are needed
                cycle
                
            end if
            
            ! check if first character is a "E" signalling the epoch name (grib-epoch)
            ! note: extra check for found flag of specifier may not be necessary, but if done so only the first appearance of the specifier will be used for reading in
            if ((first_char == 'E') .AND. (.NOT. E_found)) then
                ! change the file position backwards by one record (line) in order to read the line again
                backspace(unit= file_unit)
                
                ! read in the line containing the number of data entries using using automatic format as recommended to reduce risk of errors
                read(unit= file_unit, fmt= *, iostat= read_status) first_char, epoch_name
                
                ! set logical variable reporting the finding of the data specifier to .true.
                E_found= .TRUE.
                
                ! cycle to next loop step as value has been found and no further statements in this loop step are needed
                cycle
                
            end if
            
            ! check if first character is a "N" signalling the number of data entries (observations)
            ! note: extra check for found flag of specifier may not be necessary, but if done so only the first appearance of the specifier will be used for reading in
            if ( (first_char == 'N') .AND. (.NOT. N_found)) then
                ! change the file position backwards by one record (line) in order to read the line again
                backspace(unit= file_unit)
                
                ! read in the line containing the number of data entries using using automatic format as recommended to reduce risk of errors
                read(unit= file_unit, fmt= *, iostat= read_status) first_char, nr_data
                
                ! set logical variable reporting the finding of the data specifier to .true.
                N_found= .TRUE.
                
                ! cycle to next loop step as value has been found and no further statements in this loop step are needed
                cycle
                
            end if
            
        end do
        
        ! check if nr_data is not higher than 0 --> this would be an error
        if (nr_data <= 0) then
            ! report error message
            print '(a, a, a, /)', 'Error: No data entries found in the epolog: ', load_filename, '! Program stopped!'
            ! stop the program
            stop
        end if
        
        
        ! allocate the structure observations in the structure epolog that will store the read in data using the determined number of lines
        allocate(epolog % observations(nr_data))
        
        
        ! rewind the file
        rewind(unit= file_unit)
        
        ! initialize index for data entries
        i= 0
        
        ! loop over the file until all expected data entries have been found
        do while (i < nr_data)
            ! read in the current line using automatic format as recommended to reduce risk of errors
            read(unit= file_unit, fmt= '(a)') first_char
            
            ! check if first character is a "D" signalling a data entry
            if (first_char == 'D') then
                ! change the file position backwards by one record (line) in order to read the line again
                backspace(unit= file_unit)
                
                ! raise index for data entries
                i= i + 1
                    
                ! read in the line containing the data entry using using automatic format as recommended to reduce risk of errors
                ! note: direct read in to character array
                read(unit= file_unit, fmt= *) first_char, &
                                              epolog % observations(i) % scannr, &
                                              epolog % observations(i) % mjd, &
                                              epolog % observations(i) % year, &
                                              epolog % observations(i) % doy, &
                                              epolog % observations(i) % hour, &
                                              epolog % observations(i) % min, &
                                              epolog % observations(i) % sec, &
                                              epolog % observations(i) % station, &
                                              epolog % observations(i) % az, &
                                              epolog % observations(i) % elev, &
                                              epolog % observations(i) % source, &
                                              epolog % observations(i) % temp, &
                                              epolog % observations(i) % pres, &
                                              epolog % observations(i) % wvpr
                
            end if
            
        end do
        
    else
        ! report error message
        print '(a, a, a, /)', 'Error: Problem with opening epolog-file: ', load_path // load_filename, ' Program stopped!'
        ! stop the program
        stop
        
    end if
    
    
    ! close the file
    
    ! check if the file should be deleted and the closing
    if ( delete_file ) then
            
        ! file will be closed and deleted
        ! note: "delete" only works if the file has not been opened as open(readonly)
        close(unit= file_unit, status= 'delete')
        
        ! report error message
        print '(tr4, a)', 'Epolog-file is deleted.'
            
    else
            
        ! file will be closed, but is kept existing
        close(unit= file_unit)
        
        ! report error message
        print '(tr4, a)', 'Epolog-file is kept stored.'
        
    end if
    
    
    ! assign the session name to the epolog structure
    epolog % session_name= session_name
    
    ! assign the epoch name to the epolog structure
    epolog % epoch_name= epoch_name
    
    ! assign the total number of observations in the epolog to the epolog structure
    epolog % nr_obs= nr_data
    
end subroutine load_epologs