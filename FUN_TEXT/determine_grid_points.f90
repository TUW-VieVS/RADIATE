! module_determine_grid_points.f90

!****************************************************************************
!
! PURPOSE:  Module containing subroutine to 
!           
!           Note: Created as module to avoid explicit inferface for the subroutine as needed because of assumed-shape arrays in subroutine.
!           
!           "determine_grid_points" is a subroutine to determine the four surrounding grid points around a
!           specific point of interest (POI).
!
!           The coordinates of the POI need to fullfill the following conditions (adapt this prior to calling
!           this function):
!                  lat=[90°;-90°]
!                  lon=[0°;360°[
! 
!           The function supports the search of grid points around a POI that may lie across the
!           orginal grid borders of the input grid (0° in longitude and +-90° in latitude) in case that the
!           input is a global grid with starting value 0° in longitude and 90° in latitude.
!           Attention: Grids that are global only in one dimension are not treated as global in the following!
! 
!           The output of this function are the values (this is optional and only in case the respective parameter is present at the call)
!           and the indices of the four points in the grid surrounding the POI.
! 
! 
! INPUT:
!         POI_lat............. latitude of POI in [°]
!         POI_lon............. longitude of POI in [°]
!         dint_lat............ grid interval for latitude in [°]
!         dint_lon............ grid interval for longitude in [°]
!         grid_lat............ grid containing the latitude nodes around the station in [°]
!         grid_lon............ grid containing the longitude nodes around the station in [°]
!         grid_size........... size of the grid
!         start_and_global_check... logical value telling if input grid has the desired starting
!                                   values in latitude and longitude and if it is a global grid in
!                                   latitude and longitude.
!                                   .TRUE. ... all checks are true
!                                   .FALSE. ... at least one check is false
! 
! 
! OUTPUT:
!         optional output: in case output parameters are present at the call:
!
!         optional: lat1lon1............ [lat,lon] of point lat1lon1 in the grid in [°]
!         optional: lat1lon2............ [lat,lon] of point lat1lon2 in the grid in [°]
!         optional: lat2lon2............ [lat,lon] of point lat2lon2 in the grid in [°]
!         optional: lat2lon1............ [lat,lon] of point lat2lon1 in the grid in [°]

!         ind_lat1lon1........ index [row,column] of point lat1lon1 in the grid
!         ind_lat1lon2........ index [row,column] of point lat1lon2 in the grid
!         ind_lat2lon2........ index [row,column] of point lat2lon2 in the grid
!         ind_lat2lon1........ index [row,column] of point lat2lon1 in the grid
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
! 04.02.2015: create the Fortran-file based on the Matlab-file "determine_grid_points.m"
! 13.05.2015: comments
! 02.06.2015: correct comments
! 08.09.2015: correct wrong index calculation for lon1 in case of 0° longitude transition
!             if the 0° longitude is in between the grid and not at the start
!
!****************************************************************************

module module_determine_grid_points

contains
    
    subroutine determine_grid_points( POI_lat, &
                                      POI_lon, &
                                      dint_lat, &
                                      dint_lon, &
                                      grid_lat, &
                                      grid_lon, &
                                      grid_size, &
                                      start_and_global_check, &
                                      ind_lat1lon1, &
                                      ind_lat1lon2, &
                                      ind_lat2lon2, &
                                      ind_lat2lon1, &
                                      lat1lon1_out, & ! optional
                                      lat1lon2_out, & ! optional
                                      lat2lon2_out, & ! optional
                                      lat2lon1_out ) ! optional

        ! Define modules to be used
        
        
        !----------------------------------------------------------------------------
    
        ! Variable definitions
        implicit none
    
    
        ! dummy arguments
        !----------------
    
        ! INPUT
        
        ! input variable for ellips. latitude of a single POI
        double precision, intent(in) :: POI_lat
    
        ! input variable for ellips. longitude of a single POI
        double precision, intent(in) :: POI_lon
        
        ! input variable for storing the grid resolution in latitude
        double precision, intent(in) :: dint_lat
        
        ! input variable for storing the grid resolution in longitude
        double precision, intent(in) :: dint_lon
        
        ! input variable for storing grid of latitude values
        ! note: dimension = nlat, nlon
        double precision, dimension(:, :), intent(in) :: grid_lat
        
        ! input variable for storing grid of longitude values
        ! note: dimension = nlat, nlon
        double precision, dimension(:, :), intent(in) :: grid_lon
        
        ! input variable for storing the grid size in latitude and longitude
        integer, dimension(:), intent(in) :: grid_size
        
        ! input variable for check if starting values of first grid point (latitude and longitude) apply to requirements
        ! as well as the grid coverage overall (latitude and longitude) is global
        logical, intent(in) ::  start_and_global_check
        
    
        ! OUTPUT
        
        ! output variable for index [row,column] of point lat1lon1 in the grid
        integer, dimension(2), intent(out) :: ind_lat1lon1
        
        ! output variable for index [row,column] of point lat1lon2 in the grid
        integer, dimension(2), intent(out) :: ind_lat1lon2
        
        ! output variable for index [row,column] of point lat2lon2 in the grid
        integer, dimension(2), intent(out) :: ind_lat2lon2
        
        ! output variable for index [row,column] of point lat2lon1 in the grid
        integer, dimension(2), intent(out) :: ind_lat2lon1
        
        ! OPTIONAL OUTPUT
        
        ! optional output variable for [lat,lon] of point lat1lon1 in the grid in [°]
        double precision, dimension(2), intent(out), optional :: lat1lon1_out
        
        ! optional output variable for [lat,lon] of point lat1lon1 in the grid in [°]
        double precision, dimension(2), intent(out), optional :: lat1lon2_out
        
        ! optional output variable for [lat,lon] of point lat1lon1 in the grid in [°]
        double precision, dimension(2), intent(out), optional :: lat2lon2_out
        
        ! optional output variable for [lat,lon] of point lat1lon1 in the grid in [°]
        double precision, dimension(2), intent(out), optional :: lat2lon1_out
        
    
        ! local variables
        !----------------
    
        ! local variable for storing [lat,lon] of point lat1lon1 in the grid in [°]
        double precision, dimension(2) :: lat1lon1
        
        ! local variable for storing [lat,lon] of point lat1lon2 in the grid in [°]
        double precision, dimension(2) :: lat1lon2
        
        ! local variable for storing [lat,lon] of point lat2lon2 in the grid in [°]
        double precision, dimension(2) :: lat2lon2
        
        ! local variable for storing [lat,lon] of point lat2lon1 in the grid in [°]
        double precision, dimension(2) :: lat2lon1
        
        ! variable for storing the first latitude value in the grid (starting latitude) in [°]
        double precision :: lat_first
        
        ! variable for storing the first longitude value in the grid (starting longitude) in [°]
        double precision :: lon_first
        
        ! variable for storing the last longitude value in the grid (end longitude) in [°]
        double precision :: lon_last
        
        ! variable for storing the number of grid nodes in latitude
        integer :: nlat
        
        ! variable for storing the number of grid nodes in longitude
        integer :: nlon
        
        ! variable for storing the value of lat1 in [°]
        double precision :: lat1
        
        ! variable for storing the value of lon1 in [°]
        double precision :: lon1
        
        ! variable for storing the index of lat1
        integer :: ind_lat1
        
        ! variable for storing the index of lon1
        integer :: ind_lon1
        
        ! auxiliary variable for index of lon2 for ind_lon2_lat1
        integer :: ind_lon2_for_lat1
        
        ! variable for storing the index of lat2
        integer :: ind_lat2
        
        ! auxiliary variable for index of lon1 for ind_lon1_lat2
        integer :: ind_lon1_for_lat2
        
        ! auxiliary variable for index of lon2 for ind_lon2_lat2
        integer :: ind_lon2_for_lat2
        
        ! variable for storing the index of lon2
        integer :: ind_lon2
    
        !----------------------------------------------------------------------------
    
        !============================================================================
        ! Body of the subroutine
        !============================================================================
    
        ! get the first latitude and longitude entries in the input grid
        
        ! determine the first values in the (sub-)grid
        lat_first= grid_lat(1, 1)
        lon_first= grid_lon(1, 1)
        
        
        ! get the number of grid nodes in latitude and longitude
        
        ! determine the number of grid points in latitude and longitude
        nlat= grid_size(1)
        nlon= grid_size(2)
        
        ! get the last longitude entry in the input grid

        ! determine the first values in the (sub-)grid
        lon_last= grid_lon(1,nlon)
        
        
        ! get the values lat1 and lon1
        ! The point lat1_lon1 fullfills the following conditions compared to the POI:
        ! lat1 < POI_lat
        ! lon1 < POI_lon
        
        ! get previous latitude value in the grid (latitude below the POI)
        ! note: floor(POI/dint) determines the full number of times dint is contained in the POI value
        lat1= floor( POI_lat / dint_lat ) * dint_lat ! in [°], use floor to calculate the last (towards negative infinite) integer solution of the division --> multiplied with the resolution in latitude delivers the grid latitude below the POI
        
        ! get previous longitude value in the grid (lower longitude than the POI)
        ! note: floor(POI/dint) determines the full number of times dint is contained in the POI value
        lon1= floor( POI_lon / dint_lon ) * dint_lon ! in [°], use floor to calculate the last (towards negative infinite) integer solution of the division --> multiplied with the resolution in longitude delivers the previous grid longitude before the POI
        
        ! find the indices for lat1 in the subgrid for latitude
        ! note: nint() returns nearest integer and returns an integer value, but input must be a kind of real
        ! usage of nint() determines which grid point is the nearest to the POI
        ! note: round() is just used to be sure that an integer value as strictly expected is received (as using lat1, which is an exact grid point) --> just needed to avoid numerical problems
        !       round() is not part of the formula to determine the index, like if the nearest point in grid would be searched!
        ind_lat1= nint( (lat_first - lat1) / dint_lat ) + 1 ! +1 because otherwise only the number of intervals is determined
        
        ! find the indices for lon1 in the subgrid for longitude
        ! note: nint() returns nearest integer and returns an integer value, but input must be a kind of real
        ! usage of nint() determines which grid point is the nearest to the POI
        ! note: round() is just used to be sure that an integer value as strictly expected is received (as using lon1, which is an exact grid point) --> just needed to avoid numerical problems
        !       round() is not part of the formula to determine the index, like if the nearest point in grid would be searched!
        ind_lon1= nint( (lon1 - lon_first) / dint_lon ) + 1 ! +1 because otherwise only the number of intervals is determined
        
        ! Check if the grid contains the transition of 0° longitude somewhere in between the grid and not at
        ! the beginning.
        ! Note: Usual the global grids contain start with 0° longitude, but subgrids may have the transition
        !       across 0° longitude.
        ! note: This leads to a wrong index calculation of the longitude indices for both the subgrid and
        !       the global grid mode.
        !       The wrong index would be treated like an out of subgrid point and would be reset to the grid
        !       boundary!
        ! 
        ! Test 1:
        ! In case the first longitude is greater or equal to the last longitude then the 0° meridian is in
        ! between the grid.
        ! Note: A check for greater or equal is not needed as usually start and end value should not be
        !       the same as otherwise the data is contained twice.
        if (lon_first > lon_last) then
                
            ! Test 2:
            ! With positive test 1, the next test is to check if the searched lat1 is part of the grid on
            ! the right side meaning at a longitude after 0°, but lower or equal to the right grid boundary.
            ! In this case the previous calculated ind_lon1 must be corrected as a negative index has been
            ! calculated suggesting an out of grid index, but this is only due to the 0° longitude
            ! transition. 
            if (lon1 <= lon_last) then
                ! add the number of indices of 360/dint_lon to "ind_lon1" in order to treat the occurring 0° longitude transition 
                ind_lon1= ind_lon1 + (360 / dint_lon)
        
            ! In case the test 2 says that lon1 is larger than the right grid boundary it may have two
            ! resons:
            !   a) lon1 is really outside the grid on the rightside as lon1>lon_last_sg, but also outside
            !      the leftside as lon1<lon_first_sg
            !   b) lon1 is part of the left side of the grid before 0° longitude as lon1>lon_first_sg
            ! Therefore an additional check is needed to set a treatment in case of a), but only for grids
            ! that are not global as this case is not possible in a global grid, where lon1 must be inside
            ! the grid on one side.
            ! Case b) can be left without any treatment as the previous calculated ind_lon1 is correct for
            ! this case.
            !
            ! The next test needs to be made to distinguish to which grid boundary the lon1 value is nearer,
            ! although it is outside both sides.
            ! This test is necessary to decide to which boundary index ind_lon1 should be set to deliver a
            ! more correct positional location
            ! Use the differences between (lon_first_sg-lon1) and (lon1-lon_last_sg) to decide to which
            ! boundary lon1 is nearer
            ! In case (lon1-lon_last_sg) is smaller than (lon_first_sg-lon1) then lon1 is nearer to the
            ! right boundary, so ind_lon1 must be set to the right boundary-1. For the other case no action
            ! is necessry as ind_lon1 will be set to the left boundary automatically later due to the
            ! negative original calculated ind_lon1.
            ! 
            ! Attention: The following check only applies to grids that are not global as in global grids
            !            lon1 can never be out of bound on both sides, so if check (lon1<=lon_last_sg) fails
            !            lon1 must be >= lon_first_sg leading to directly correct results for ind_lon1.
            else if ( ( .NOT. start_and_global_check) .AND. ( (lon_first - lon1) > (lon1 - lon_last) ) ) then
                ! set ind_lon1 to the right boundary
                ! note: nlon-1 to be able to set ind_lon2= nlon
                ind_lon1= nlon - 1
            end if
        end if
        
        
        ! correct lat1 and lon1 if necessary and get the values lat2 and lon2
        
        ! in case that it is a global grid starting at 0° longitude and 90° latitude
        if (start_and_global_check) then
            
            ! check if indices for lat1 and lon1 are correct
                
            ! in case that input grid is global (grids that are global only in one dimension are not treated
            ! as global)
            ! --> nothing needs to be done because the indices of lat1 and lon1 will always be found correctly
            !     at least due to the previous correction step!
            
            ! combine indices for row and columns
            ind_lat1lon1=[ ind_lat1, ind_lon1 ]
            ! get the values for lat1_lon1
            lat1lon1= [ grid_lat(ind_lat1, 1), grid_lon(1, ind_lon1) ]
            
            
            ! get correct indices for lat1lon2, lat2lon2 and lat2lon1 separately
            ! Attention: Because the value of lon2 would be different for the points lat1lon2 and lat2lon2
            ! in case that the pole is crossed, it is necessary to determine the values for each of the
            ! remaing 3 points separately!!!
            
            
            ! get the value for lat1lon2
            ! check if the index of ind_lon1 is the last index in the grid
            ! --> crossing the zero-meridian
            if ( ind_lon1 == nlon ) then
                ! get correct ind_lon2_for_lat1
                ind_lon2_for_lat1= 1 ! --> set ind_lon2 to the first entry
            else ! no zero-meridian crossing
                ! get correct ind_lon2_for_lat1
                ind_lon2_for_lat1= ind_lon1 + 1 ! raise ind_lon1 by 1 to get ind_lon2
            end if
            
            ! ind_lat1lon2
            ! combine indices for row and columns
            ind_lat1lon2= [ ind_lat1, ind_lon2_for_lat1 ]
            ! get the values for lat1_lon2
            lat1lon2= [ grid_lat(ind_lat1, 1), grid_lon(1, ind_lon2_for_lat1) ]
            
            
            ! case of pole crossing:
            ! get the values for lat2lon1 and lat2lon2
            ! check if the index of ind_lat1 is the first index in the grid
            ! note: in this case the pole is crossed!
            if (ind_lat1 == 1) then
                ! get correct ind_lat2
                ind_lat2= 2 ! --> set ind_lat2 to the second entry as the "higher" latitude across the pole is actually the next lower one (lat1=90° --> lat2=89° if dint_lat=1°)
                
                ! get correct ind_lon1_for_lat2
                ! change longitude value of lon1 by 180° due to crossing the pole
                ! correct longitude: -180°, if longitude is equal or above 180°
                if ( lat1lon1(2) >= 180 ) then ! check value of lon1
                    ind_lon1_for_lat2= ind_lon1 - 180 / dint_lon
                ! correct longitude: +180°, if longitude is below 180°
                else if ( lat1lon1(2) < 180 ) then ! check value of lon1
                    ind_lon1_for_lat2= ind_lon1 + 180 / dint_lon
                end if
                
                ! get correct ind_lon2_for_lat2
                ! change longitude value of lon2 by 180° due to crossing the pole
                ! correct longitude: -180°, if longitude is equal or above 180°
                if ( lat1lon2(2) >= 180 ) then ! check value of lon2
                    ind_lon2_for_lat2= ind_lon2_for_lat1 - 180 / dint_lon
                ! correct longitude: +180°, if longitude is below 180°
                else if ( lat1lon2(2) < 180 ) then ! check value of lon2
                    ind_lon2_for_lat2= ind_lon2_for_lat1 + 180 / dint_lon
                end if
                
            else ! no pole crossing
                
                ! get correct ind_lat2
                ind_lat2= ind_lat1 - 1 ! -1 because the node lies "above" in the grid therefore the index is lower
                
                ! get correct ind_lon1_for_lat2
                ind_lon1_for_lat2= ind_lon1 ! index for lon1 stays the same if pole is not crossed
                
                ! get correct ind_lon2_for_lat2
                ind_lon2_for_lat2= ind_lon2_for_lat1 ! index for lon2 stays the same if pole is not crossed
                
            end if
            
            ! ind_lat2lon1
            ! combine indices for row and columns
            ind_lat2lon1= [ ind_lat2, ind_lon1_for_lat2 ]
            ! get the values for lat1_lon2
            lat2lon1= [ grid_lat(ind_lat2, 1), grid_lon(1, ind_lon1_for_lat2) ]
            
            ! ind_lat2lon2
            ! combine indices for row and columns
            ind_lat2lon2= [ ind_lat2, ind_lon2_for_lat2]
            ! get the values for lat1_lon2
            lat2lon2= [ grid_lat(ind_lat2, 1), grid_lon(1, ind_lon2_for_lat2) ]
            
            
        else ! in any other case = no global grid in latitude and longitude or wrong starting values in latitude or longitude
            
            ! check if indices for lat1 and lon1 are correct
            
            ! check if found indices for lat1 and lon1 are out of bounds and reset them to the edges of the grid
            ! --> check if indices of lat1 and lon1 are set correctly and do not cross borders of input grid
            
            ! if the index of latitude is smaller than the second entry in the grid (index 1 is necessary for lat2)
            if (ind_lat1 < 2) then
                ind_lat1= 2
            end if
            
            ! if the index of latitude is higher than the last index in the grid
            if (ind_lat1 > nlat) then
                ind_lat1= nlat
            end if
            
            ! if the index of longitude is smaller than the first index in the grid
            if (ind_lon1 < 1) then
                ind_lon1= 1
            end if
            
            ! if the index of longitude is higher than the second last index in the grid (last index is used for ind_lon2)
            if (ind_lon1 > nlon-1) then
                ind_lon1= nlon - 1
            end if
            
            
            ! get correct indices for lat2 and lon2 without the possibility of crossing grid borders
            
            ! get the grid-index for the lat2 using the above determined index for lat1
            ind_lat2= ind_lat1 - 1 ! -1 because the node lies "above" in the grid therefore the index is lower
            
            ! get the grid-index for lon2 using the above determined index for lon1
            ind_lon2= ind_lon1 + 1 ! raise ind_lon1 by 1 to get ind_lon2
            
            
            ! assign values (values and indices each with combined latitude and longitude domain)
            
            ! determine values for lat1lon1, lat1lon2, lat2lon2 and lat2lon1 using the found indices
            lat1lon1= [ grid_lat(ind_lat1, 1), grid_lon(1, ind_lon1) ]
            lat1lon2= [ grid_lat(ind_lat1, 1), grid_lon(1, ind_lon2) ]
            lat2lon2= [ grid_lat(ind_lat2, 1), grid_lon(1, ind_lon2) ]
            lat2lon1= [ grid_lat(ind_lat2, 1), grid_lon(1, ind_lon1) ]
            
            ! combine indices for row and columns
            ind_lat1lon1= [ ind_lat1, ind_lon1 ]
            ind_lat1lon2= [ ind_lat1, ind_lon2 ]
            ind_lat2lon2= [ ind_lat2, ind_lon2 ]
            ind_lat2lon1= [ ind_lat2, ind_lon1 ]
            
        end if
        
        
        ! determine which optional output arguments are needed
        if ( present(lat1lon1_out) ) then
            lat1lon1_out= lat1lon1
        end if
        
        if ( present(lat1lon2_out) ) then
            lat1lon2_out= lat1lon2
        end if
        
        if ( present(lat2lon2_out) ) then
            lat2lon2_out= lat2lon2
        end if
        
        if ( present(lat2lon1_out) ) then
            lat2lon1_out= lat2lon1
        end if
        
        
    end subroutine determine_grid_points        
    
end module module_determine_grid_points