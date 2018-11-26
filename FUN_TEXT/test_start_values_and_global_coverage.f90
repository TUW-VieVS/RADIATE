! test_start_values_and_global_coverage.f90

!****************************************************************************
!
! PURPOSE:  Subroutine to test input grid on global coverage
!            
!           "test_start_values_and_global_coverage" is a subroutine to test if an input grid delivers data
!           starting at a specific latitude and longitude value and it tests if the input data has a global
!           coverage in latitude and longitude.
!
!
! INPUT:
!         nlat.......... number of points along a latitude circle in the grid
!         nlon.......... number of points along a longitude meridian in the grid
!         dint_lat...... grid intervall in the grib-file for latitude in [°]
!         dint_lon...... grid intervall in the grib-file for longitude in [°]
!         lat_first..... latitude of first grid point in the grib-file
!         lon_first..... longitude of first grid point in the grib-file
!         reference_start_lat... value against which the first latitude value of the grid should be
!                                tested
!         reference_start_lon... value against which the first longitude value of the grid should be
!                                tested
! 
! 
! OUTPUT:
!         global_lat_check......... logical value telling if input grid is global in latitude (1... global, 0... not global)
!         global_lon_check......... logical value telling if input grid is global in longitude (1... global, 0... not global)
!         start_lat_check.......... logical value telling if input grid has the desired starting value in latitude (1... correct starting value, 0... other starting value)
!         start_lon_check.......... logical value telling if input grid has the desired starting value in longitude (1... correct starting value, 0... other starting value)
!         start_and_global_check... Logical value telling if all pre-test have been true. If so, then this variable delivers true (=1) otherwise false (=0)
!
!----------------------------------------------------------------------------
! 
! Fortran-file created by Armin Hofmeister
!   based on Matlab scripts created by Armin Hofmeister
! 
!----------------------------------------------------------------------------
! History:
! 
! 12.01.2015: create the Fortran-file based on the Matlab-file "test_start_values_and_global_coverage.m"
! 28.01.2015: comments
! 02.06.2015: correct comments
!
!****************************************************************************
    
    
subroutine test_start_values_and_global_coverage( nlat, &
                                                  nlon, &
                                                  dint_lat, &
                                                  dint_lon, &
                                                  lat_first, &
                                                  lon_first, &
                                                  reference_start_lat, &
                                                  reference_start_lon, &
                                                  global_lat_check, &
                                                  global_lon_check, &
                                                  start_lat_check, &
                                                  start_lon_check, &
                                                  start_and_global_check )

    ! Define modules to be used
    
    !----------------------------------------------------------------------------
    
    ! Variable definitions
    implicit none
    
    
    ! dummy arguments
    !----------------
    
    ! INPUT
    
    ! variables for storing the number of grid points in latitude and longitude
    integer, intent(in) :: nlat, nlon
    
    ! variables for storing the grid resolution in latitude and longitude
    double precision, intent(in) :: dint_lat, dint_lon
    
    ! variable for storing the latitude of the first grid point
    double precision, intent(in) :: lat_first
    
    ! variable for storing the longitude of the first grid point
    double precision, intent(in) :: lon_first
    
    ! define variables for storing the reference values of the desired starting position of the grid in latitude and longitude
    double precision, intent(in) :: reference_start_lat, reference_start_lon
    
    
    ! OUTPUT
    
    ! define control variable for check if grid coverage in latitude is global
    logical, intent(out) :: global_lat_check
    
    ! define control variable for check if grid coverage in longitude is global
    logical, intent(out) :: global_lon_check
    
    ! define control variable for check if latitude of first grid point is equal to requirement
    logical, intent(out) :: start_lat_check
    
    ! define control variable for check if longitude of first grid point is equal to requirement
    logical, intent(out) :: start_lon_check
    
    ! define control variable for check if starting values of first grid point (latitude and longitude) apply to requirements
    ! as well as the grid coverage overall (latitude and longitude) is global
    logical, intent(out) :: start_and_global_check
    
    
    !----------------------------------------------------------------------------
    
    !============================================================================
    ! Body of the subroutine
    !============================================================================
    
    
    ! check if it is a global grid in latitude
    ! use the number of latitude grid points and the latitude grid resolution to check for global grid
    ! in latitude direction
    ! a global latitude grid has one entry more in nlat than needed to cover 180°, because of the equator --> (nlat-1)*dint_lat should be 180
    global_lat_check= ((nlat-1)*dint_lat==180) ! global latitude grid: .TRUE. ... global, .FALSE. ... not global
    
    ! check if it is a global grid in longitude
    ! a global grid doesn't report lon_last as 360° (=0°) --> lon_last would then be 360°- [grid interval in longitude dint_lon]
    ! use the number of longitude grid points and the longitude grid resolution to check for global grid
    ! in longitude direction
    global_lon_check= (nlon*dint_lon==360); ! global longitude grid: .TRUE. ... global, .FALSE. ... not global
    
    ! check if grid starts with latitude value as specified by reference_start_lat (usually 90°)
    ! this check is needed to control the correct assignment of the undulations in a later step
    start_lat_check= (lat_first==reference_start_lat) ! start_lat_check: .TRUE. ... correct starting value, .FALSE. ... other starting value
    
    ! check if grid starts with longitude value as specified by reference_start_lon (usually 0°)
    ! this check is needed to control the correct assignment of the undulations in a later step
    start_lon_check= (lon_first==reference_start_lon) ! start_lon_check: .TRUE. ... correct starting value, .FALSE. ... other starting value
                                    
    
    ! check if the grid is global in latitude and longitude and starting with correct latitude and
    ! longitude
    start_and_global_check= (global_lat_check .AND. global_lon_check .AND. start_lat_check .AND. start_lon_check) ! note: result of .AND.-expression is only .TRUE. in case all connected items are .TRUE.

    
end subroutine test_start_values_and_global_coverage