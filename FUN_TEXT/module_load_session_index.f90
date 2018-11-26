! module_load_session_index.f90

!****************************************************************************
!
! PURPOSE:  Module containing subroutine to import the filenames of all epolog-files for the specific session created with "create_session_index_and_epologs"
!            
!           "load_session_index" loads the filenames of all epolog textfiles that are necessary for the
!           specific session (AZEL-file)
!
!           The output "epolog_filenames" is a character string array containing all epolog filenames for the specific session.
!           Note: As the creation of the epologs is done chronologically also the epolog filenames in the textfile should have chronological order!
! 
!
! INPUT:
!         load_path....... specifies path to the loading directory
!         load_filename... filename of the session index
!         delete_file..... variable for determining if file should be deleted at the closing:
!                          .TRUE. .... delete the file
!                          .FALSE. ... keep the file
! 
! 
! OUTPUT:
!         epolog_filenames... character string array containing all epolog filenames for one specific
!                             session defined through the load_filename
!         nr_files........... number of epologs that make up the specific session
!
!----------------------------------------------------------------------------
! 
! Fortran-file created by Armin Hofmeister
!   based on Matlab scripts created by Armin Hofmeister
! 
!----------------------------------------------------------------------------
! History:
! 
! 20.11.2014: create the Fortran-file based on the Matlab-file "load_session_index.m"
! 24.11.2014: programming
! 25.11.2014: comments
! 29.01.2015: change to module as to save explicit interface block
! 20.04.2015: correct comments and messages
! 07.05.2015: add possibility to delete the loaded session index-file after the loading when it is closed
! 31.08.2015: correct comments
! 10.09.2015: change the string length of the output variable "epolog_filenames" from dynamic allocation to
!             to avoid errors with gfortran fixed length
!             use constants for string length determination of epolog filenames
!
!****************************************************************************

module module_load_session_index

contains
    
    subroutine load_session_index( load_path, &
                                   load_filename, &
                                   delete_file, &
                                   epolog_filenames, &
                                   nr_files )

        ! Define modules to be used
        
        use module_constants, only: len_epolog_filename
    
    
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
        
        ! define variable for storing the station names
        ! note: see "module_constants" for "len_epolog_filename"
        character(len=len_epolog_filename), dimension(:), allocatable, intent(out) :: epolog_filenames
    
        ! define variable for storing the number of epologs that make up the specific session
        integer, intent(out) :: nr_files
    
        
        ! local variables
        !----------------
    
        ! variable with the unit of the file
        integer, parameter :: file_unit=1
    
        ! variable for storing the opening status of the file
        integer :: open_status
    
        ! variable for storing the reading status of the file
        integer :: read_status
    
        ! variable for storing the first character of a line of the file
        character :: first_char
    
        ! index for data entries
        integer :: i
        
        
        ! CONSTANTS
        
        ! see "module_constants" for "len_epolog_filename"
    
    
        !----------------------------------------------------------------------------
    
        !============================================================================
        ! Body of the subroutine
        !============================================================================
    
        ! Columns in the session index file:
        !    0     .... data specifier sign: S= session name, N= number of data entries (filenames), D= data entry (filename)
        !	 1     .... filename of epolog
    
    
        ! Open and read the session index textfile
        open(unit= file_unit, file= load_path // load_filename, action= 'read', status= 'old', iostat= open_status)
    
        ! Check if textfile can be opened (iostat == 0 --> no error)
        if (open_status == 0) then
        
            ! Read in the file once in order to determine the number of filename entries
            ! This is necessary for allocating the variable for storing the filenames.
        
            ! Initialize the variable of number of data entries (filenames) in the file
            nr_files= 0
        
            ! loop do-if-exit over the file
            do
                ! read next line
                read(unit= file_unit, fmt= '(a)', iostat= read_status) first_char
            
                ! test if end of file is reached and exit loop in this case
                if (is_iostat_end(read_status)) exit
            
                ! check if first character is a "N" signalling the number of data entries (filenames)
                if (first_char == 'N') then
                    ! change the file position backwards by one record (line) in order to read the line again
                    backspace(unit= file_unit)
                
                    ! read in the line containing the number of data entries using using automatic format as recommended to reduce risk of errors
                    read(unit= file_unit, fmt= *, iostat= read_status) first_char, nr_files
                
                    ! exit the loop as task of finding the number of data entries has been acomplished
                    exit
                
                end if
            
            end do
        
            ! check if nr_files is not higher than 0 --> this would be an error
            if (nr_files <= 0) then
                ! report error message
                print '(a, a, a, /)', 'Error: No filenames found in the session index! ', load_filename, '! Program stopped!'
                ! stop the program
                stop
            end if
        
            ! rewind the file
            rewind(unit= file_unit)
        
            
            ! allocate the variable for storing the epolog filenames
            allocate(epolog_filenames(nr_files))
            
            ! initialize index for data entries
            i= 0
        
            ! loop over the file until all expected data entries have been found
            do while (i < nr_files)
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
                    read(unit= file_unit, fmt= *) first_char, epolog_filenames(i)
                                
                end if
            
            end do
        
        else
            ! report error message
            print '(a, a, a, /)', 'Error: Problem with opening session index-file: ', load_path // load_filename, ' Program stopped!'
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
            print '(tr4, a)', 'Session index-file is deleted.'
            
        else
            
            ! file will be closed, but is kept existing
            close(unit= file_unit)
            
            ! report error message
            print '(tr4, a)', 'Session index-file is kept stored.'
        
        end if

    end subroutine load_session_index
    
end module module_load_session_index