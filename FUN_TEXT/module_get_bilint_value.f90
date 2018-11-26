! module_get_bilint_value.f90

!****************************************************************************
!
! PURPOSE:  Module containing subroutine to bilinear interpolate values at a specific point of interest (POI).
!           
!           Note: Created as module to avoid explicit inferface for the subroutine as needed because of assumed-shape arrays in subroutine.
!           
!           "get_bilint_value" is a subroutine to determine bilinear interpolated values at a specific point of
!           interest (POI) for three different variables using the four surrounding points of the POI from a grid.
!
!           Note: There are some cases of optional input and output parameters:
!                 1) already known indices of surrounding points: If input of "ind_lat1lon1_in", "ind_lat1lon2_in", "ind_lat2lon2_in" and "ind_lat2lon1_in" is provided then the input of
!                           "grid_lat", "grid_lon", "grid_size" and "start_and_global_check" is not needed as the determination of the indices of the
!                           surrounding grid points is skipped.
!                 2) unknown indices of surrounding points: If no input of "ind_lat1lon1_in", "ind_lat1lon2_in", "ind_lat2lon2_in" and "ind_lat2lon1_in" is provided then the input of
!                           "grid_lat", "grid_lon", "grid_size" and "start_and_global_check" is required as the determination of the indices of the
!                           surrounding grid points needs to be done using the subroutine "determine_grid_points".
!                 3) output of indices of surrounding points: the output of the indices of the surrounding grid points around the POI is optional. This means that for (only for) each
!                           provided output parameter the value will be reported.
! 
!      
! INPUT:
!         v1.................. grid "v1": gridded values with the dimensions: (level,lat,lon) with which bilinear interpolation should be done
!         v2.................. grid "v2": gridded values with the dimensions: (level,lat,lon) with which bilinear interpolation should be done
!         v3.................. grid "v3": gridded values with the dimensions: (level,lat,lon) with which bilinear interpolation should be done
!         level............... index of (height) level for interpolation
!         POI_lat............. latitude of POI in [°]
!         POI_lon............. longitude of POI in [°]
!         dint_lat............ grid interval for latitude in [°]
!         dint_lon............ grid interval for longitude in [°]
!
!         optional input: not needed, if indices of surrounding grid points are provided as input:
!
!         optional: grid_lat............ grid containing the latitude nodes around the station in [°]
!         optional: grid_lon............ grid containing the longitude nodes around the station in [°]
!         optional: grid_size........... size of the grid
!         optional: start_and_global_check... logical value telling if input grid has the desired starting
!                                             values in latitude and longitude and if it is a global grid in
!                                             latitude and longitude.
!                                             .TRUE. ... all checks are true
!                                             .FALSE. ... at least one check is false
!
!         optional input: in this case no determination of indices is made
!
!         optional: ind_lat1lon1_in........ index [row,column] of point lat1lon1 in the grid
!         optional: ind_lat1lon2_in........ index [row,column] of point lat1lon2 in the grid
!         optional: ind_lat2lon2_in........ index [row,column] of point lat2lon2 in the grid
!         optional: ind_lat2lon1_in........ index [row,column] of point lat2lon1 in the grid
! 
! 
! OUTPUT:
!         v1_bilint............ bilinear interpolated value from "v1" at the POI
!         v2_bilint............ bilinear interpolated value from "v2" at the POI
!         v3_bilint............ bilinear interpolated value from "v3" at the POI
!
!         optional output: these parameters will only be part of the output if their variable is present in the call of the subroutine:
!
!         optional: ind_lat1lon1_out........ index [row,column] of point lat1lon1 in the grid
!         optional: ind_lat1lon2_out........ index [row,column] of point lat1lon2 in the grid
!         optional: ind_lat2lon2_out........ index [row,column] of point lat2lon2 in the grid
!         optional: ind_lat2lon1_out........ index [row,column] of point lat2lon1 in the grid
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
! 03.02.2015: create the Fortran-file based on the Matlab-file "get_bilint_value.m"
! 04.02.2015: programming
! 20.04.2015: correct comments and messages
! 13.05.2015: comments
! 02.06.2015: correct comments
!
!****************************************************************************

module module_get_bilint_value

contains
    
    subroutine get_bilint_value( v1, &
                                 v2, &
                                 v3, &
                                 level, &
                                 POI_lat, &
                                 POI_lon, &
                                 dint_lat, &
                                 dint_lon, &
                                 grid_lat, & ! optional: see cases where required
                                 grid_lon, & ! optional: see cases where required
                                 grid_size, & ! optional: see cases where required
                                 start_and_global_check, & ! optional: see cases where required
                                 ind_lat1lon1_in, & ! optional: see cases
                                 ind_lat1lon2_in, & ! optional: see cases
                                 ind_lat2lon2_in, & ! optional: see cases
                                 ind_lat2lon1_in, & ! optional: see cases
                                 v1_bilint, &
                                 v2_bilint, &
                                 v3_bilint, &
                                 ind_lat1lon1_out, & ! optional: see cases
                                 ind_lat1lon2_out, & ! optional: see cases
                                 ind_lat2lon2_out, & ! optional: see cases
                                 ind_lat2lon1_out ) ! optional: see cases
    
        ! Define modules to be used
        use module_determine_grid_points
        
        
        !----------------------------------------------------------------------------
    
        ! Variable definitions
        implicit none
    
    
        ! dummy arguments
        !----------------
    
        ! INPUT
    
        ! input variable for grid "v1": gridded values with the dimensions: (level,lat,lon) with which bilinear interpolation should be done
        double precision, dimension(:, :, :), intent(in) :: v1
        
        ! input variable for grid "v2": gridded values with the dimensions: (level,lat,lon) with which bilinear interpolation should be done
        double precision, dimension(:, :, :), intent(in) :: v2
        
        ! input variable for grid "v3": gridded values with the dimensions: (level,lat,lon) with which bilinear interpolation should be done
        double precision, dimension(:, :, :), intent(in) :: v3
        
        ! input variable for index of (height) level for interpolation
        integer, intent(in) :: level
        
        ! input variable for ellips. latitude of a single POI
        double precision, intent(in) :: POI_lat
    
        ! input variable for ellips. longitude of a single POI
        double precision, intent(in) :: POI_lon
        
        ! input variable for storing the grid resolution in latitude
        double precision, intent(in) :: dint_lat
        
        ! input variable for storing the grid resolution in longitude
        double precision, intent(in) :: dint_lon
    
        ! OPTIONAL INPUT
        ! note: optional, but definitely required in case no input support of "ind_lat1lon1_in", "ind_lat1lon2_in", "ind_lat2lon2_in" and "ind_lat2lon1_in"
        
        ! optional input variable for storing grid of latitude values
        ! note: dimension = nlat, nlon
        double precision, dimension(:, :), intent(in), optional :: grid_lat
        
        ! optional input variable for storing grid of longitude values
        ! note: dimension = nlat, nlon
        double precision, dimension(:, :), intent(in), optional :: grid_lon
        
        ! optional input variable for storing the grid size in latitude and longitude
        integer, dimension(:), intent(in), optional :: grid_size
        
        ! optional input variable for check if starting values of first grid point (latitude and longitude) apply to requirements
        ! as well as the grid coverage overall (latitude and longitude) is global
        logical, intent(in), optional ::  start_and_global_check
        
        
        ! OPTIONAL INPUT
        ! note: optional, but definitely required in case no input support of "grid_lat", "grid_lon" "grid_size" and "start_and_global_check"
        
        ! optional input variable for index [row,column] of point lat1lon1 in the grid
        integer, dimension(2), intent(in), optional :: ind_lat1lon1_in
        
        ! optional input variable for index [row,column] of point lat1lon2 in the grid
        integer, dimension(2), intent(in), optional :: ind_lat1lon2_in
        
        ! optional input variable for index [row,column] of point lat2lon2 in the grid
        integer, dimension(2), intent(in), optional :: ind_lat2lon2_in
        
        ! optional input variable for index [row,column] of point lat2lon1 in the grid
        integer, dimension(2), intent(in), optional :: ind_lat2lon1_in
        
    
        ! OUTPUT
    
        ! output variable for bilinear interpolated value from "v1" at the POI
        double precision, intent(out) :: v1_bilint
        
        ! output variable for bilinear interpolated value from "v2" at the POI
        double precision, intent(out) :: v2_bilint
        
        ! output variable for bilinear interpolated value from "v3" at the POI
        double precision, intent(out) :: v3_bilint
        
        ! OPTIONAL OUTPUT
        
        ! optional output variable for index [row,column] of point lat1lon1 in the grid
        integer, dimension(2), intent(out), optional :: ind_lat1lon1_out
        
        ! optional output variable for index [row,column] of point lat1lon2 in the grid
        integer, dimension(2), intent(out), optional :: ind_lat1lon2_out
        
        ! optional output variable for index [row,column] of point lat2lon2 in the grid
        integer, dimension(2), intent(out), optional :: ind_lat2lon2_out
        
        ! optional output variable for index [row,column] of point lat2lon1 in the grid
        integer, dimension(2), intent(out), optional :: ind_lat2lon1_out
        
    
        ! local variables
        !----------------
    
        ! local variable for index [row,column] of point lat1lon1 in the grid
        integer, dimension(2) :: ind_lat1lon1
        
        ! local variable for index [row,column] of point lat1lon2 in the grid
        integer, dimension(2) :: ind_lat1lon2
        
        ! local variable for index [row,column] of point lat2lon2 in the grid
        integer, dimension(2) :: ind_lat2lon2
        
        ! local variable for index [row,column] of point lat2lon1 in the grid
        integer, dimension(2) :: ind_lat2lon1
        
        
        ! local variable for value of grid point lat1lon1 in grid v1
        double precision :: v1_lat1lon1
        
        ! local variable for value of grid point lat1lon2 in grid v1
        double precision :: v1_lat1lon2
        
        ! local variable for value of grid point lat2lon2 in grid v1
        double precision :: v1_lat2lon2
        
        ! local variable for value of grid point lat2lon1 in grid v1
        double precision :: v1_lat2lon1
        
        
        ! local variable for value of grid point lat1lon1 in grid v2
        double precision :: v2_lat1lon1
        
        ! local variable for value of grid point lat1lon2 in grid v2
        double precision :: v2_lat1lon2
        
        ! local variable for value of grid point lat2lon2 in grid v2
        double precision :: v2_lat2lon2
        
        ! local variable for value of grid point lat2lon1 in grid v2
        double precision :: v2_lat2lon1
        
        
        ! local variable for value of grid point lat1lon1 in grid v3
        double precision :: v3_lat1lon1
        
        ! local variable for value of grid point lat1lon2 in grid v3
        double precision :: v3_lat1lon2
        
        ! local variable for value of grid point lat2lon2 in grid v3
        double precision :: v3_lat2lon2
        
        ! local variable for value of grid point lat2lon1 in grid v3
        double precision :: v3_lat2lon1
    
        !----------------------------------------------------------------------------
    
        !============================================================================
        ! Body of the subroutine
        !============================================================================
        
        ! Search for correct indices
        
        ! in case all indices (ind_lat1lon1_in, ind_lat1lon2_in, ind_lat2lon2_in, ind_lat2lon1_in) are provided as input estimation of the indices is not necessary
        if ( present(ind_lat1lon1_in) .AND. present(ind_lat1lon2_in) .AND. present(ind_lat2lon2_in) .AND. present(ind_lat2lon1_in) ) then
            ! just assign the input indices to the local variables
            ind_lat1lon1= ind_lat1lon1_in
            ind_lat1lon2= ind_lat1lon2_in
            ind_lat2lon2= ind_lat2lon2_in
            ind_lat2lon1= ind_lat2lon1_in
        
        ! in case that not all indices are provided as input, they have to be determined
        else
            
            ! check if all here needed optional inputs "grid_lat", "grid_lon" "grid_size" and "start_and_global_check" are present
            if ( present(grid_lat) .AND. present(grid_lon) .AND. present(grid_size) .AND. present(start_and_global_check) ) then
                
                ! get indices of the four grid points around the POI in the original grid
                ! note: see subroutine for optional arguments (output of latlon values: here not needed)
                call determine_grid_points( POI_lat, &
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
                                            ind_lat2lon1 )
                
            ! if at least one needed optional input is missing
            else
                ! report error message
                write(unit= *, fmt= '(a)') 'Error in "get_bilint_value": Missing of optional, but in this way of calling needed input parameter(s)! Program stopped!'
                ! stop the program
                stop
                
            end if
            
        end if
        
        
        ! Bilinear Interpolation for values of the grids "v1", "v2" and "v3"
        
        ! grid "v1"
        
        ! get the values from "v1" for the determined nodes used for the bilinear interpolation at the POI
        v1_lat1lon1= v1(level, ind_lat1lon1(1), ind_lat1lon1(2))
        v1_lat1lon2= v1(level, ind_lat1lon2(1), ind_lat1lon2(2))
        v1_lat2lon2= v1(level, ind_lat2lon2(1), ind_lat2lon2(2))
        v1_lat2lon1= v1(level, ind_lat2lon1(1), ind_lat2lon1(2))
        
        ! bilinear interpolation for deriving value at the POI in grid "v1"
        call bilinear_interpolation( POI_lat, &
                                     POI_lon, &
                                     dint_lat, &
                                     dint_lon, &
                                     v1_lat1lon1, &
                                     v1_lat1lon2, &
                                     v1_lat2lon2, &
                                     v1_lat2lon1, &
                                     v1_bilint )
        
        
        ! grid "v2"
        
        ! get the values from "v2" for the determined nodes used for the bilinear interpolation at the POI
        v2_lat1lon1= v2(level, ind_lat1lon1(1), ind_lat1lon1(2))
        v2_lat1lon2= v2(level, ind_lat1lon2(1), ind_lat1lon2(2))
        v2_lat2lon2= v2(level, ind_lat2lon2(1), ind_lat2lon2(2))
        v2_lat2lon1= v2(level, ind_lat2lon1(1), ind_lat2lon1(2))
        
        ! bilinear interpolation for deriving value at the POI in grid "v2"
        call bilinear_interpolation( POI_lat, &
                                     POI_lon, &
                                     dint_lat, &
                                     dint_lon, &
                                     v2_lat1lon1, &
                                     v2_lat1lon2, &
                                     v2_lat2lon2, &
                                     v2_lat2lon1, &
                                     v2_bilint )
        
        
        ! grid "v3"
        
        ! get the values from "v3" for the determined nodes used for the bilinear interpolation at the POI
        v3_lat1lon1= v3(level, ind_lat1lon1(1), ind_lat1lon1(2))
        v3_lat1lon2= v3(level, ind_lat1lon2(1), ind_lat1lon2(2))
        v3_lat2lon2= v3(level, ind_lat2lon2(1), ind_lat2lon2(2))
        v3_lat2lon1= v3(level, ind_lat2lon1(1), ind_lat2lon1(2))
        
        ! bilinear interpolation for deriving value at the POI in grid "v3"
        call bilinear_interpolation( POI_lat, &
                                     POI_lon, &
                                     dint_lat, &
                                     dint_lon, &
                                     v3_lat1lon1, &
                                     v3_lat1lon2, &
                                     v3_lat2lon2, &
                                     v3_lat2lon1, &
                                     v3_bilint )
    
        
        ! determine if indices are subject of output
        ! if present then assign to output variable
        if (present(ind_lat1lon1_out)) then
            ind_lat1lon1_out= ind_lat1lon1
        end if
        
        if (present(ind_lat1lon2_out)) then
            ind_lat1lon2_out= ind_lat1lon2
        end if
        
        if (present(ind_lat2lon2_out)) then
            ind_lat2lon2_out= ind_lat2lon2
        end if
        
        if (present(ind_lat2lon1_out)) then
            ind_lat2lon1_out= ind_lat2lon1
        end if
        
    end subroutine get_bilint_value
    
end module module_get_bilint_value