! get_global_undulation.f90

!****************************************************************************
!
! PURPOSE:  Subroutine to to load the global geoid undulation values with a specific
!           resolution in latitude and longitude.
!
!           The outputs of this function are the geoid undulation values for a global
!           grid with latitude and longitude resolution specified through the input variables.
! 
!
! INPUT:
!         input_path_undulation... path to the directory with the global geoid undulation files with
!                                  different grid resolutions
!         dint_lat................ grid interval (resolution) in the grib-file for latitude in [°]
!         dint_lon................ grid interval (resolution) in the grib-file for longitude in [°]
!         nlat.................... number of points along a latitude circle in the grib-file
!         nlon.................... number of points along a (longitude) meridian in the grib-file
! 
! 
! OUTPUT:
!         undulation_global_grid...... global geoid undulations in [m] in grid format dimension(nlat, nlon)
!
!----------------------------------------------------------------------------
! 
! Fortran-file created by Armin Hofmeister
!   based on Matlab scripts created by Armin Hofmeister
! 
!----------------------------------------------------------------------------
! History:
! 
! 19.01.2015: create the Fortran-file based on the Matlab-file "get_global_undulation.m"
! 20.01.2015: refine error message
! 21.01.2015: comments
! 20.04.2015: correct comments and messages
! 22.04.2015: correct comments
! 03.06.2015: change to implied loop for reading in the data lines
!
!****************************************************************************
    
    
subroutine get_global_undulation( input_path_undulation, &
                                  dint_lat, &
                                  dint_lon, &
                                  nlat, &
                                  nlon, &
                                  undulation_global_grid )

    ! Define modules to be used
    
    !----------------------------------------------------------------------------
    
    ! Variable definitions
    implicit none
    
    
    ! dummy arguments
    !----------------
    
    ! INPUT
    
    ! define variable for the path to the directory with the global geoid undulation files
    ! with different grid resolutions
    character(len=*), intent(in) :: input_path_undulation
    
    ! define variable for the grid interval (resolution) in the grib-file for latitude in [°]
    double precision, intent(in) :: dint_lat
    
    ! define variable for the grid interval (resolution) in the grib-file for longitude in [°]
    double precision, intent(in) :: dint_lon
    
    ! define variable for the number of points along a latitude circle in the grib-file
    integer, intent(in) :: nlat
    
    ! define variable for the number of points along a (longitude) meridian in the grib-file
    integer, intent(in) :: nlon
    
    
    ! OUTPUT
    
    ! define variable for storing the gridded undulation data
    ! note: the dimension is determined by (nlat, nlon)
    double precision, dimension(nlat, nlon), intent(out) :: undulation_global_grid
    
    
    ! local variables
    !----------------
    
    ! variable for storing the latitude resolution as a string
    character(len=5) :: dint_lat_str
    
    ! variable for storing the longitude resolution as a string
    character(len=5) :: dint_lon_str
    
    ! define variable for storing the filename of the undulation-textfile
    ! note: length is set to 50, this combines text and fixed resolution output, remember this when changing desired filename
    character(len=50) :: filename_undulation
    
    ! variable with the unit of the file
    integer, parameter :: file_unit=1
    
    ! variable for storing the opening status of the file
    integer :: open_status
    
    ! variable for storing the reading status of the file
    integer :: read_status
    
    ! variable for storing the number of lines in the file
    integer :: nr_lines
    
    ! loop index
    integer :: i
    
    ! temporal variables for latitude and longitude
    double precision :: temp_lat, temp_lon
    
    ! variable for storing the undulation data as vector when read in
    ! note: the size is determined by nlat * nlon
    double precision, dimension(nlat * nlon) :: undulation_global_vec
    
    ! variable for temporal storage of the gridded undulations after reshaping the vector, but prior to the tranposition
    ! note: allocatable attribute is not neccessary as size can be determined at declaration, but then the variable can not be deallocated after use
    double precision, dimension(:, :), allocatable :: temp_undulation_global_grid
    
    
    !----------------------------------------------------------------------------
    
    !============================================================================
    ! Body of the subroutine
    !============================================================================
    
    
    ! Load global geoid undulations with the desired resolution in latitude and longitude
    
    ! Define the correct filename of the file where the geoid undulations are stored in the desired
    ! latitude/longitude resolution
    
    ! convert input resolutions in latitude and longitude from numbers to to character strings
    ! convert dint_lat to a string of 5 characters in the form x.yyy (fmt: ss --> no plus sign, f5.3 --> 1 leading zero in case of shorter value)
    write (unit= dint_lat_str, fmt='(ss, f5.3)') dint_lat
    ! convert dint_lon to a string of 5 characters in the form x.yyy (fmt: ss --> no plus sign, f5.3 --> 1 leading zero in case of shorter value)
    write (unit= dint_lon_str, fmt='(ss, f5.3)') dint_lon
    
    ! create filename
    ! note: length is set to 50, this combines text and fixed resolution output, remember this when changing desired filename
    filename_undulation = 'global_undulations_dint_lat' // dint_lat_str // '_dint_lon' // dint_lon_str // '.txt'
    
    
    ! open and read the textfile "filename_undulation", which contains the global geoid undulations in
    ! the desired latitude/longitude resolution
    open(unit= file_unit, file= input_path_undulation // filename_undulation, action= 'read', status= 'old', iostat= open_status)
    
    ! check if textfile can be opened (iostat == 0 --> no error)
    if (open_status == 0) then
                    
        ! import the data using textscan (automatic detection of different entries) with the correct data format
        ! columns:
        !         1... latitude in [°]
        !         2... longitude in [°]
        !         3... geoid undulation in [m]
        
        
        ! read the file
        
        ! Read in the file once in order to determine the number of lines = number of station entries
        ! This is necessary for allocating the variables for storing the data.
        
        ! Initialize the counter of number of lines in the file
        nr_lines= 0
        
        ! loop do-if-exit
        do
            ! read next line
            read(unit= file_unit, fmt= '(a)', iostat= read_status)
            
            ! test if end of file is reached and exit loop in this case
            if (is_iostat_end(read_status)) exit
            
            ! raise counter for number of lines
            nr_lines= nr_lines + 1
        end do
        
        ! check if the nr_lines are equal to the necessary number of data = nlat * nlon that need to be present in the file
        ! if true, proceed with reading the data
        if (nr_lines == (nlat * nlon)) then
            ! rewind the file
            rewind(unit= file_unit)
                
            ! import the data
            ! note: Implied loop for reading the consecutive data lines is used as it should be faster than a normal do loop.
            ! note: skip latitude and longitude columns by just using temporal variables that will be overwritten in each loop
            !       to save RAM when working with high resolutions
            read(unit= file_unit, fmt= *) ( temp_lat, temp_lon, undulation_global_vec(i), i= 1, nr_lines ) ! undulations as vector in [m]
                
        else
            ! report error message
            print '(a /)', 'Error: Number of lines in the undulation-file is not equal to the expected number of geoid undulations! Program stopped!'
            ! stop the program
            stop           
        end if
        
    else
        ! report error message
        print '(a)', 'Error: Problem with opening textfile for global geoid undulations!'
        print '(tr4, 5a)', 'Make sure that there really exists a textfile for the global geoid undulations with the needed resolutions: res_lat= ', dint_lat_str, ', res_lon= ', dint_lon_str, '!'
        print '(tr4, a /)', 'Program stopped!'
        ! stop the program
        stop
        
    end if
    
    ! close the textfile
    close(unit= file_unit)
    
    
    ! create grid out of the vectorized grid points for undulation
    ! note: reshape fills array columnwise
    ! note: temporal variable is necessary to split step of reshape and transpose to avoid stack overflow
    
    ! allocate temporal variable
    ! note: prior to transpose the dimensions of the grid are (nlon, nlat)
    allocate(temp_undulation_global_grid(nlon, nlat))
    
    ! reshape the undulations from a vector to a rectangular grid
    temp_undulation_global_grid= reshape(undulation_global_vec, [nlon, nlat])
    
    ! transpose temporal variable and assign values to receive output in dimension(nlat, nlon)
    undulation_global_grid= transpose(temp_undulation_global_grid) ! transpose is necessary to get grid format according to grib-file grid!!!
    
    ! deallocate temporal variable
    deallocate(temp_undulation_global_grid)
    
    end subroutine get_global_undulation