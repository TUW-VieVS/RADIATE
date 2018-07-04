!  create_azel_name_list.f90 
!
!  This file reads the azel-file creation specifications from the file 
! azel_spec.txt and creates a text-file containing a list with all azel-files
! during the specified period. As a next step, this list can be read by a 
! batch-file in order to create ray-tracing for all these files.
!

!----------------------------------------------------------------------------
! 
! Fortran-file created by Daniel Landskron
!   based on Matlab scripts created by Daniel Landskron
!
!----------------------------------------------------------------------------
! History:
! 
! 2017/11/20: create the Fortran-file 
! 
!****************************************************************************

    
program create_azel_name_list

    
    !------------------------------------------------------------------------
    ! Define modules to be used
    
    ! Module for type definitions
    use module_date_type_definition
    
    ! module for the subroutine 'date2mjd'
    use module_date2mjd
    
    ! module for the subroutine 'mjd2date'
    use module_mjd2date
    
    !------------------------------------------------------------------------
    
    implicit none
    
    !------------------------------------------------------------------------
    ! Define variables
    
    ! define variable for storing the file unit of the azel_spec file (unit is always the same)
    integer, parameter :: file_unit_azel_spec = 1
    
    ! define variable for the azel_spec filename
    character(len=:), allocatable :: load_path_filename_azel_spec
    
    ! variable for storing the opening status of the file
    integer :: open_status

    ! variable for storing the first character of a line of the file
    character :: first_char
    
    ! variable for storing the reading status of the file
    integer :: read_status
    
    ! define variables for the read in dates
    character(len=1000) :: start_date_str
    character(len=1000) :: end_date_str
    
    ! variable for the number of epochs
    integer :: num_epochs
    
    ! variable for storing the epoch civil date for calculating the epoch in mjd
    type(date_type) :: date_of_epoch_start, date_of_epoch_end, date_of_epoch
    
    ! variable for the mjd's
    double precision :: mjd_start, mjd_end
    
    ! variables for the mjd's
    integer :: i_mjd, n_mjd
    
    ! variable for mjd
    double precision, dimension(:), allocatable :: mjd
    
    ! variables for the date_str specification
    character(len=4) :: year_str
    character(len=2) :: month_str
    character(len=2) :: day_str
    character(len=2) :: hour_str
    
    ! define variable for storing the file unit of the azel_list file (unit is always the same)
    integer, parameter :: file_unit_azel_list = 2
    
    ! define variable for the saving path of the azel_list.txt
    character(len=:), allocatable :: save_path_filename_azel_list
    
    
    !============================================================================
    ! Body of the program
    !============================================================================
    
    
    
    ! define the input path and output path + files
    load_path_filename_azel_spec = '../DATA/INPUT/AZEL_spec.txt'
    save_path_filename_azel_list = '../DATA/INPUT/AZEL_list.txt'
    
    ! Open the file azel_spec.txt
    open(unit= file_unit_azel_spec, file= load_path_filename_azel_spec, action= 'read', status= 'old', iostat= open_status)
    
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
        
        ! 1.) this argument corresponds to the start date
        read(unit= file_unit_azel_spec, fmt='(A)') start_date_str
         
        ! read first character of the next line, and if a comment, then jump over it to the next line
        do
            read(unit= file_unit_azel_spec, fmt= '(a)', iostat= read_status) first_char
            if (first_char /= '%' .and. first_char /= '!' .and. first_char /= '#') then
                backspace(unit= file_unit_azel_spec)
                exit
            end if
        end do
                
        ! 2.) this argument corresponds to the end date
        read(unit= file_unit_azel_spec, fmt='(A)') end_date_str
        
        ! read first character of the next line, and if a comment, then jump over it to the next line
        do
            read(unit= file_unit_azel_spec, fmt= '(a)', iostat= read_status) first_char
            if (first_char /= '%' .and. first_char /= '!' .and. first_char /= '#') then
                backspace(unit= file_unit_azel_spec)
                exit
            end if
        end do
        
        ! 3.) this argument corresponds to the number of epochs
        read(unit= file_unit_azel_spec, fmt=*) num_epochs

        if (num_epochs /= 4) then
            print '(a, a, a /)', 'The script determine_trop.f90 is hardcoded for num_epochs=4 (through the 0.25d0). ', &
                  'If this number changes, then determine_trop.f90 must first be adapted, most likely through hardcoding for ', &
                  'num_epochs=8. If this is done, then adjust this error message here! Program stopped!'
            stop
        end if
            
        ! close the file
        close(unit= file_unit_azel_spec)
            
    else
        
        ! report error message
        print '(a, a, a /)', 'Error: Problem with opening the azel specifications file "', load_path_filename_azel_spec, &
                             '" for observation information! Program stopped!'
        ! stop the program
        stop
        
    end if
    
    ! convert the date_strings to integer
    read(start_date_str(1:4),*)  date_of_epoch_start % year
    read(end_date_str(1:4),*)    date_of_epoch_end   % year
    read(start_date_str(6:7),*)  date_of_epoch_start % month
    read(end_date_str(6:7),*)    date_of_epoch_end   % month
    read(start_date_str(9:10),*) date_of_epoch_start % day
    read(end_date_str(9:10),*)   date_of_epoch_end   % day
    date_of_epoch_start % hour = 0d0
    date_of_epoch_end   % hour = 0d0
    date_of_epoch_start % min = 0d0
    date_of_epoch_end   % min = 0d0
    date_of_epoch_start % sec = 0.d0
    date_of_epoch_end   % sec = 0.d0
    
    
    ! get mjd from date
    call date2mjd(date_of_epoch_start, mjd_start)
    call date2mjd(date_of_epoch_end, mjd_end)
    
    ! if mjd_end < mjd_start, then report an error message
    if (mjd_end < mjd_start) then
        print '(a /)', 'Error: mjd_end must not be smaller than mjd_start! Program stopped!'
    end if
        
    ! determine the number of mjd's, which is consequently also the number of azel-files
    n_mjd = (mjd_end-mjd_start+1)*num_epochs
    allocate(mjd(n_mjd))
    
    
    ! open the azel_list.txt file
    open(unit= file_unit_azel_list, file= save_path_filename_azel_list , &
         action= 'write', status= 'replace', iostat= open_status)
    ! Check if textfile can't be opened (iostat /= 0 --> error)
    if (open_status /= 0 ) then
        ! report error message
        print '(a, a, a /)', 'Error: Problem with creating the file azel_list.txt. Program stopped!'
        ! stop the program
        stop
    end if
   
    
    ! determine the list of mjd's and consequently the azel-filenames and write them to the output file
    do i_mjd= 1,n_mjd
        
        ! determine the respective mjd
        if (i_mjd==1) then
            mjd(i_mjd) = mjd_start
        else
            mjd(i_mjd) = mjd(i_mjd-1) + 1.0/num_epochs
        end if
        
        ! convert the mjd to a date vector
        call mjd2date(mjd(i_mjd), date_of_epoch)
        
        ! convert the date vector elements to strings
        write(unit=year_str  , fmt='(ss, i4.4)') date_of_epoch % year
        write(unit=month_str , fmt='(ss, i2.2)') date_of_epoch % month
        write(unit=day_str   , fmt='(ss, i2.2)') date_of_epoch % day
        write(unit=hour_str  , fmt='(ss, i2.2)') date_of_epoch % hour
        
        ! write the azel-filename to the output text file
        write(unit= file_unit_azel_list, fmt= '(a)') 'azel_' // year_str // month_str // day_str // hour_str //'_UNI.txt'
        
    end do
    
    ! close the file
    close(unit= file_unit_azel_list)
    
    
end program create_azel_name_list

