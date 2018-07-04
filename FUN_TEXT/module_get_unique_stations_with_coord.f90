! module_get_unique_stations_with_coord.f90

!****************************************************************************
!
! PURPOSE:  Module containing subroutine to determine the unique station names and their coordinates and additional data
!           of stations carrying out observations.
!            
!           "get_unique_stations_with_coord" is a subroutine to determine from a bunch of station data the unique stations according to the station names
!           and their according data values like coordinates and additional data.
!
!           Attention: Within this subroutine the input station names and data will be sorted alphabetically after the station name.
!                      Otherwise the unique station determination method would fail. So after sorting the same station names are grouped
!                      (eg. BADARY, BADARY, ..., BADARY, FORTLEZA, ..., FORTLEZA).
!
!           Note: Because of the alphabetically sorting of the input stations, the output of the unique stations and their according data
!                 is also alphabetically sorted after the station names.
! 
!           Note: Created as module to avoid explicit inferface for the subroutine as needed because of assumed-shape dummy argument in subroutine.
!
! 
! INPUT:
!         dupl_stat........ structure of type "observing_stations_type" containing the observing station data including possible duplicates
!                           The structure contains the following variables:
!
!           % station.............. variable for storing the station name
!           % nr_obs_per_obsstat... variable for storing the number of observations per station, e.g. in a specific epolog
!           % suspend_raytr........ variable for checking if ray-tracing should (needs to) be suspended for the observing station
!           % lat_ell.............. variable for storing the ellipsoidal latitude coordinate
!           % lon_ell.............. variable for storing the ellipsoidal longitude coordinate
!           % h_ell................ variable for storing the ellipsoidal height coordinate
! 
! 
! OUTPUT:
!         unique_stat...... structure of type "observing_stations_type" containing the unique observing station data
!                           The structure contains the following variables:
!
!           % station.............. variable for storing the station name
!           % nr_obs_per_obsstat... variable for storing the number of observations per station, but here used to store the number of times this unique station
!                                   is observing station in the epologs
!           % suspend_raytr........ variable for checking if ray-tracing should (needs to) be suspended for the observing station,
!                                   here: value is .TRUE. if it was .TRUE. for any entry of this station in an epolog
!           % lat_ell.............. variable for storing the ellipsoidal latitude coordinate
!           % lon_ell.............. variable for storing the ellipsoidal longitude coordinate
!           % h_ell................ variable for storing the ellipsoidal height coordinate
!
!
!----------------------------------------------------------------------------
! 
! Fortran-file created by Armin Hofmeister
!   based on Fortran scripts created by Armin Hofmeister
! 
!----------------------------------------------------------------------------
! History:
! 
! 12.05.2015: create the Fortran-file based on the Fortran-file "module_get_observing_stations.m"
! 13.05.2015: comments
!
!****************************************************************************

module module_get_unique_stations_with_coord

contains
    
    subroutine get_unique_stations_with_coord( dupl_stat, &
                                               unique_stat )

        ! Define modules to be used
    
        use module_type_definitions, only: observing_stations_type
        use module_msort ! module for sorting character strings or double or integer values ascending (alphabetical)
        
        ! note: usage of "module_sort_type_definitions" not necessary as use declared in module for sorting and types are public
        !use module_sort_type_definitions ! module for type definitions used for sorting with module_qsort and module_msort
    
        !----------------------------------------------------------------------------
    
        ! Variable definitions
        implicit none
    
    
        ! dummy arguments
        !----------------
    
        ! INPUT
        
        ! define structure for the input of station data containing possible duplicates
        type(observing_stations_type), dimension(:), intent(in) :: dupl_stat
    
        
        ! OUTPUT
        
        ! define structure for the output of the unique stations and their data
        ! note: see "module_constants" for length of station names
        type(observing_stations_type), dimension(:), allocatable, intent(out) :: unique_stat
    
    
        ! local variables
        !----------------
    
        ! define variable for the number of entries in "dupl_stat" (number of stations including possible duplicates)
        integer :: nr_total
    
        ! define structure that is used for sorting after station names (see module_sort_type_definitions for definition)
        type(sort_char_type), dimension(:), allocatable :: sort_struct_char
    
        ! define structure that is used for sorting after station names (see module_sort_type_definitions for definition)
        type(sort_char_type), dimension(:), allocatable :: sort_struct_char_temp
        
        ! define structure for the alphabetically sorted input data
        ! note: size is the same size as "dupl_stat"
        ! note: This variable is needed as the original input variable can not be newly defined as it has the intent(in) attribute.
        type(observing_stations_type), dimension(size(dupl_stat)) :: dupl_stat_sorted
        
        ! define variable for loop index
        integer :: i
        
        ! define variable for index of current observation that is looked at
        integer :: ind_entry
    
        ! define logical vector for storing comparison results of station names
        ! note: size is the number of stations in "dupl_stat" including possible duplicates
        logical, dimension(size(dupl_stat)) :: check_for_curr_stat
    
        ! define variable for storing the index of the observing (unique) station name
        integer :: ind_unique
        
        ! define logical vector for the values of "suspend_raytr" at all appearances of the current compared station
        logical, dimension(:), allocatable :: curr_suspend_raytr
        
        
        !----------------------------------------------------------------------------
    
        !============================================================================
        ! Body of the subroutine
        !============================================================================
    
        ! Sort the input station data with possibly contained duplicates alphabetically after the station name
        
        ! determine the number of entries in "dupl_stat" (number of stations possibly containing duplicates)
        nr_total= size(dupl_stat)
        
        ! allocate the structures used for sorting (see module_msort for help)
        allocate(sort_struct_char(nr_total))
        allocate(sort_struct_char_temp((nr_total+1)/2))
    
        ! assign station names and indices to the structure for sorting
        do i= 1, nr_total
            sort_struct_char(i) % value= dupl_stat(i) % station
            sort_struct_char(i) % order= i
        end do
    
        ! sort the list with the station names
        ! use mergesort for sorting (note: original order is preserved in case of equal entries; mergesort is a stable sorting algorithm)
        call msort(sort_struct_char, nr_total, sort_struct_char_temp)
    
        ! sort all variables in the structure according to the new alphabetical order
        ! note: also the "station"-field needs to be sorted as it has not been done yet
        
        ! assign data in sorted order to structure "dupl_stat_sorted"
        ! note: Variable "dupl_stat_sorted" is needed as input variable "dupl_stat" can not be newly defined because of the inten(in) attribute.
        ! note: Furthermore an assignment of the type a= a(order) would lead to a stack overflow in case of a large dataset, so a temporal variable
        !       would be needed anyway.
        dupl_stat_sorted = dupl_stat(sort_struct_char % order)
        
        ! deacllocate the structure used for sorting after station names
        deallocate(sort_struct_char)
        deallocate(sort_struct_char_temp)
        
        
        !-----------------------------------------------------------------------------------------------
        
        ! Get all unique stations (unique station names and their coordinates)
    
        ! set index for the current entry to 1
        ind_entry= 1
    
        ! set index for unique stations (start at 0!!!)
        ind_unique= 0
    
    
        ! find names of all unique stations and determine the number of appearances per station
        ! loop over the entries in "dupl_stat_sorted"
        ! note: At this stage entries in "dupl_stat_sorted" must be sorted alphabetically after the station names.
        do while (ind_entry <= nr_total) ! do while the input list has not reached the end
            ! compare first station name entry of "dupl_stat_sorted" to all entries in "dupl_stat_sorted" and receive
            ! logical array with .TRUE. at equal names
            check_for_curr_stat= (dupl_stat_sorted(ind_entry) % station == dupl_stat_sorted % station)
        
            ! raise index for unique station
            ind_unique=ind_unique + 1
        
            ! assign name of current compared station as new found unique station
            ! note: Assigning the current found unique station name and following the according station data to the structure 
            !       currently used for the unique station determination is possible as this can overwrite only an entry that has
            !       already been treated and whose data will not be needed or searched for any more.
            !       Also the logical check "check_for_curr_stat" for the number of appearances of the current station is not affected
            !       as only that station data is overwritten and newly written that will not be caught by this check any more in further loop cycles.
            !       All this holds true only in case the data has been sorted alphabetically after the station name before. Otherwise the algorithm will
            !       not work anyway and overwriting would destroy needed data.
            dupl_stat_sorted(ind_unique) % station= dupl_stat_sorted(ind_entry) % station
            
            ! assign according coordinates of current compared station to the new found unique station
            dupl_stat_sorted(ind_unique) % lat_ell= dupl_stat_sorted(ind_entry) % lat_ell
            dupl_stat_sorted(ind_unique) % lon_ell= dupl_stat_sorted(ind_entry) % lon_ell
            dupl_stat_sorted(ind_unique) % h_ell= dupl_stat_sorted(ind_entry) % h_ell
            
            ! determine number of total appearances for this station
            ! note: count determines the number of true elements
            dupl_stat_sorted(ind_unique) % nr_obs_per_obsstat= count(check_for_curr_stat)
            
            ! determine if ray-tracing has been suspended for any of the entries of the current unique station
            ! note: set value for "suspend_raytr" to .TRUE. if any of the entries for the current unique station is .TRUE.
            
            ! check if the current station has been found more than once
            ! note: this is necessary as the function any() needs the input of an array to work
            if ( dupl_stat_sorted(ind_unique) % nr_obs_per_obsstat > 1 ) then
                
                ! note: direct use of any() with the logical variable "" as indices is not possible
                
                ! allocate the variable for the values of "suspend_raytr" at all appearances of the current compared station
                ! note: size is determined through the number of appearances of the current compared station
                allocate(curr_suspend_raytr(dupl_stat_sorted(ind_unique) % nr_obs_per_obsstat) )
                
                ! get the values for "suspend_raytr" at all appearances
                ! note: Use the "check_for_curr_stat"-variable for determining for which entries the value should be taken.
                !       This mask leads to assigning only the .TRUE. elements to the result.
                curr_suspend_raytr= pack( dupl_stat_sorted % suspend_raytr, check_for_curr_stat )
                
                ! use any() and assign to structure
                dupl_stat_sorted(ind_unique) % suspend_raytr= any( curr_suspend_raytr )
                
                ! deallocate "curr_suspend_raytr" in order to be allocatable again for the next loop cycle
                deallocate(curr_suspend_raytr)
                
            ! in case there is only one appearance for the current unique station
            else
                ! assign the value of the index "ind_entry" as this is the only appearance of the current compared station
                dupl_stat_sorted(ind_unique) % suspend_raytr= dupl_stat_sorted(ind_entry) % suspend_raytr
            end if
            
            
            ! determine index of next station in duplicate list that is not equal to the current station
            ! note: in order to do this correct the input list of station names from the observations has to be sorted alphabetically as stated before
            !       the next entry is the one after the last .TRUE. entry (additive list of last ind_entry value + new .TRUE. entries)
            ind_entry= ind_entry + dupl_stat_sorted(ind_unique) % nr_obs_per_obsstat
        
        end do
    
        ! get the number of different observing stations
        ! note: use current value of "ind_unique" for this as it reports the total number of unique stations
        !nr_unique= ind_unique
    
        ! allocate the output structure according to the found number of different observing stations
        allocate(unique_stat(ind_unique))
    
        ! assign the unique station data to the output structure
        ! note: use the number of found unique stations for getting the upper boundary
        unique_stat= dupl_stat_sorted(1:ind_unique)
    
    end subroutine get_unique_stations_with_coord
    
end module module_get_unique_stations_with_coord