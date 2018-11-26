! get_gribdata_txt.f90

!****************************************************************************
!
! PURPOSE:  Subroutine to import the grib-file data that have been stored in a textfile (by MATLAB function).
! 
!
! INPUT:
!         load_path.......   specifies path to the loading directory
!         load_filename...   filename of the epolog
! 
! 
! OUTPUT:
!         gribdata........  structure array containing the grib-data
!
!
!----------------------------------------------------------------------------
! 
! Fortran-file created by Armin Hofmeister
! 
!----------------------------------------------------------------------------
! History:
! 
! 08.01.2015: create the Fortran-file
! 12.01.2015: programming
! 19.01.2015: correct right alignment of rows, delete usage of read_status as not needed
! 20.04.2015: correct comments and messages
! 25.08.2015: correct comments
! 26.08.2015: correct comments
! 21.12.2015: add check of correct representation of the year of the grib epoch by 4 digits
!             and add overall check of correct loaded grib epoch through comparison of
!             grib-file filename and epoch time from grib data
! 13.01.2016: correct comments
!
!****************************************************************************
    
    
subroutine get_gribdata_txt( load_path, &
                             load_filename, &
                             gribdata )

    ! Define modules to be used
    use module_type_definitions, only: gribdata_type
    use module_constants, only: len_year, len_month, len_day, len_hour
    
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
    
    
    ! OUTPUT
    
    ! define structure for storing the grib-data retrieved from the textfile
    type(gribdata_type), intent(out) :: gribdata
    
    
    ! local variables
    !----------------
    
    ! variable with the unit of the file
    integer, parameter :: file_unit=1
    
    ! variable for storing the opening status of the file
    integer :: open_status
    
    ! loop index for pressure level
    integer :: ind_lev
    
    ! define variables for storing the epoch time extracted from the grib-epoch name from the filename
    integer :: year_filename, month_filename, day_filename, hour_filename
    
    
    !----------------------------------------------------------------------------
    
    !============================================================================
    ! Body of the subroutine
    !============================================================================
    
    
    ! Load grib-data from textfile
    
    ! Open and read the grib-textfile "filename" which contains the data necessary for ray-tracing that have been extracted from a grib-file
    open(unit= file_unit, file= load_path // load_filename, action= 'read', status= 'old', iostat= open_status)
    
    ! Check if textfile can be opened (iostat == 0 --> no error)
    if (open_status == 0) then
        
        ! read lines
        ! additional text after the values will be ignored as reading is set to default advancing (note: list direct input with asterisk does not support non-advancing option)
            
        ! first line: nlat, nlon
        read(unit= file_unit, fmt= *) gribdata % nlat, gribdata % nlon
            
        ! second line: lat_first, lat_last
        read(unit= file_unit, fmt= *) gribdata % lat_first, gribdata % lat_last
            
        ! third line: lon_first, lon_last
        read(unit= file_unit, fmt= *) gribdata % lon_first, gribdata % lon_last
            
        ! fourth line: dint_lat, dint_lon
        read(unit= file_unit, fmt= *) gribdata % dint_lat, gribdata % dint_lon
            
        ! fifth line: nr_pres_lev
        read(unit= file_unit, fmt= *) gribdata % nr_pres_lev
            
            
        ! allocate pressure variable
        allocate(gribdata % p(gribdata % nr_pres_lev))
            
            
        ! sixth line: pressure
        read(unit= file_unit, fmt= *) gribdata % p
            
        ! seventh line: nr_param
        read(unit= file_unit, fmt= *) gribdata % nr_param
            
            
        ! check if number of parameters is 3 (Z, Q, T records required)
        ! note: this subroutine only works if exactly the 3 different parameters necessary for ray-tracing are contained in the textfile
        !       errors would occur in case more or less than the Z, Q and T records are contained
        if (gribdata % nr_param /= 3) then
            ! report error message
            print '(a, /)', 'Error: Number of parameters in grib-textfile is not 3 as required! Program stopped!'
            ! stop the program
            stop 
        end if
        
            
        ! eighth line: time of grib epoch
        read(unit= file_unit, fmt= *) gribdata % epoch % year, gribdata % epoch % month, gribdata % epoch % day, gribdata % epoch % hour, gribdata % epoch % min, gribdata % epoch % sec
        
        ! correct the year of the epoch if necessary
        ! Note: Years before 2001 are displayed without the century, i.e.
        !       1999 = 99 and 2000 = 100 !!!
        if (gribdata % epoch % year <= 100) then
            gribdata % epoch % year = gribdata % epoch % year + 1900
        end if
        
        
        ! check if the contained epoch in the file is the same as desired by comparing the time to the grib-file filename
        
        ! extract the epoch time from the grib-file filename
        ! note: The grib epoch is the filename without the extension.
        !       Thus, "len_year", "len_month", "len_day" and "len_hour" can be used to extract the year, month, day and epoch hour
        !       from the filename by correct indexing.
        read(unit= load_filename(1 : len_year), fmt=*) year_filename
        read(unit= load_filename(len_year + 1 : len_year + len_month), fmt=*) month_filename
        read(unit= load_filename(len_year + len_month + 1 : len_year + len_month + len_day), fmt=*) day_filename
        read(unit= load_filename(len_year + len_month + len_day + 1 : len_year + len_month + len_day + len_hour), fmt=*) hour_filename
        
        ! check if the desired grib epoch from the filename is not the same as contained in the file
        if  ( (year_filename /= gribdata % epoch % year) .OR. (month_filename /= gribdata % epoch % month) .OR. (day_filename /= gribdata % epoch % day) .OR. (hour_filename /= gribdata % epoch % hour) ) then
                    
            ! report error message
            print '(a, a, a, i0.2, a, i0.2, a, i0.4, a, i0.2, a /)', 'Error: Grib-file "', load_filename, '" contains data for a wrong epoch (', gribdata % epoch % day, '.', gribdata % epoch % month, '.', gribdata % epoch % year, ' ', gribdata % epoch % hour, 'h). Program stopped!'
            ! stop the program
            stop
         
        end if
        
        
        ! read in the grib-records for the specific pressure levels and parameters
        
        ! allocate the variables used for storing the parameter records
        ! note: the size is dependent on the number of pressure levels and the total numbero grid points = nlat * nlon
        allocate(gribdata % Z(gribdata % nr_pres_lev, gribdata % nlat * gribdata % nlon), &
                 gribdata % Q(gribdata % nr_pres_lev, gribdata % nlat * gribdata % nlon), &
                 gribdata % T(gribdata % nr_pres_lev, gribdata % nlat * gribdata % nlon))
        
        
        ! loop over all pressure levels
        do ind_lev= 1, gribdata % nr_pres_lev
                
            ! read the Z record for the current pressure level
            read(unit= file_unit, fmt= *) gribdata % Z(ind_lev, :)
                
            ! read the Q record for the current pressure level
            read(unit= file_unit, fmt= *) gribdata % Q(ind_lev, :)
                
            ! read the T record for the current pressure level
            read(unit= file_unit, fmt= *) gribdata % T(ind_lev, :)
                
        end do
        
        
    else
        ! report error message
        print '(a, a, a, /)', 'Error: Problem with opening grib-data textfile: ', load_path // load_filename, ' Program stopped!'
        ! stop the program
        stop
        
    end if
    
    
    ! close the file
    close(unit= file_unit)
    
    end subroutine get_gribdata_txt