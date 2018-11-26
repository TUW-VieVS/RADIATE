! module_combine_and_sort_rd.f90
    

!****************************************************************************
!
! PURPOSE:  Module containing subroutine to combine and sort all observations and the calculated
!           delays from all epologs that set up one complete session.
!           The sorting is done in a way to get the observations in the order of:
!                                            1.) chronological mjd
!                                            2.) if same mjd then alphabetical order of station name
!                                            3.) if same mjd and same station name then alphabetical order of source names
!
!           Note: Created as module to avoid explicit inferface for the subroutine as needed because of assumed-shape arrays in subroutine.
!
!
! INPUT:
!        epolog......................... structure of single epoch specific ray-tracing data and calculations (comprehension of substructures)
!
!
! OUTPUT:
!        rdlog_epoch... structure containing the following variables and substructures:
!
!           % session_name... name of the session to which the data belongs
!
!           % observations... structure containing the observation data in the following variables:
!               % scannr... scannumber
!               % mjd...... modified julian date
!               % year..... year
!               % doy...... day of year
!               % hour..... hour
!               % min...... minute
!               % sec...... second
!               % station.. station name
!               % az....... azimuth in [rad]
!               % elev..... elevation in [rad]
!               % source... source name
!               % temp..... temperature in [°C]
!               % pres..... pressure in [hPa]
!               % wvpr..... water vapour pressure in [hPa]
!
!           % delay.......... structure containing the delay data in the following variables:
!               % dz_total........ zenith total delay in [m]
!               % dz_h............ zenith hydrostatic delay in [m]
!               % dz_w............ zenith wet delay in [m]
!               % ds_total_geom... slant total delay including geometric bending effect in [m]
!               % ds_total........ slant total delay in [m]
!               % ds_h_geom....... slant hydrostatic delay including geometric bending effect in [m]
!               % ds_h............ slant hydrostatic delay in [m]
!               % ds_w............ slant wet delay in [m]
!               % e_stat.......... elevation angle at the station in [rad]
!               % e_outgoing_rt... (iteratively) ray-traced outgoing elevation angle in [rad]
!               % diff_e.......... difference: outgoing (theoretical) - outgoing (ray-traced) elevation angle in [rad]
!               % dgeo............ geometric bending effect in [m]
!               % mf_total_geom... value for total mapping factor (includes treatment of geometric bending effect)
!               % mf_h_geom....... value for hydrostatic mapping factor (includes treatment of geometric bending effect)
!               % mf_w............ value for wet mapping factor
!               % break_elev...... logical signaling if a break in the while loop for calculating the
!                                  outgoing elevation angle has occured (.FLASE.= no break, .TRUE.= break)
!               % break_layer..... logical signaling if a break in the while loop for calculating the
!                                  next intersection point has occured (at least for one intersection point) (.FLASE.= no break, .TRUE.= break)
!
!           % meteo_stat_out... structure containing the meteorological parameters at the station position from the NWM:
!               % p...... value for the pressure value at the station position in [hPa]
!               % T...... value for the temperature value at the station position in [K]
!               % wvpr... value for the water vapour pressure value at the station position in [hPa]
!
!           % grib_epoch_mjd... vector containing the mjd of the grib-file epoch which has been used to calculate the delay of each observation
!
!----------------------------------------------------------------------------
! 
! Fortran-file created by Armin Hofmeister
!   based on Matlab scripts created by Armin Hofmeister
! 
!----------------------------------------------------------------------------
! History:
! 
! 21.04.2015: create the Fortran-file based on the Matlab-file "combine_and_sort_calc_sessionwise.m"
! 22.04.2015: programming
! 23.04.2015: comments
! 07.05.2015: correct comments
! 12.05.2015: comments
! 13.05.2015: correct comments
! 25.06.2015: add sorting of combined observation data alphabetically after source names as first step
!             of sorting to get the duplicates consecutively
! 17.12.2015: correct comments due to adding of the total mapping factor calculation to the "delay" substructure
!             changes due due to adding of the meteorological data at the station position from the NWM
! 11.01.2016: correct comments due to adding of the water vapour pressure from the azel-file to the "observations" substructure
! 12.01.2016: add comments
!
!****************************************************************************

module module_combine_and_sort_rd

contains
    
    subroutine combine_and_sort_rd( epolog, &
                                    rdlog_epoch )
    
        ! Define modules to be used
        use module_type_definitions, only : epolog_type, rdlog_epoch_type
        use module_msort ! module for sorting character strings or double or integer values ascending (alphabetical)
        
        ! note: usage of "module_sort_type_definitions" not necessary as use declared in module for sorting and types are public
        !use module_sort_type_definitions ! module for type definitions used for sorting with module_qsort and module_msort
    
    
        !----------------------------------------------------------------------------
    
        ! Variable definitions
        implicit none
    
    
        ! dummy arguments
        !----------------
    
        ! INPUT
    
        ! variable for input of structure of single epoch specific ray-tracing data and calculations (comprehension of substructures)
        type(epolog_type), dimension(:), intent(in) :: epolog
    
    
        ! OUTPUT
    
        ! variable for storing the combination of all observations and their delays (including duplicate observations with delays calculated at other epochs)
        type(rdlog_epoch_type), intent(out) :: rdlog_epoch
    
    
        ! local variables
        !----------------
        
        ! define variable for storing the number of epologs for the session
        integer :: nr_epologs
        
        ! define loop variable for storing the epolog index
        integer :: curr_epolog_ind
        
        ! define variable for storing the total number of observations (including possible duplicates) and delays in all epologs
        integer :: nr_obs_total
        
        ! define variable for storing the index of the beginning position where to add epolog data to the combined structure
        integer :: ind_begin
        
        ! define variable for storing the index of the end position where to add epolog data to the combined structure
        integer :: ind_end
        
        ! variable for temporal storing the combination of all observations and their delays
        ! note: this temporal variable is needed as implicit reordering the entries using the found sorting order leads to stack overflow error for large datasets
        type(rdlog_epoch_type) :: temp_rdlog_epoch
        
        ! define variable for loop index
        integer :: i
        
        ! define structure that is used for sorting after station names (see module_sort_type_definitions for definition)
        type(sort_char_type), dimension(:), allocatable :: sort_struct_char
    
        ! define structure that is used for sorting after station names (see module_sort_type_definitions for definition)
        type(sort_char_type), dimension(:), allocatable :: sort_struct_char_temp
        
        ! define structure that is used for sorting after mjd (see module_sort_type_definitions for definition)
        type(sort_double_type), dimension(:), allocatable :: sort_struct_double
    
        ! define structure that is used for sorting after mjd (see module_sort_type_definitions for definition)
        type(sort_double_type), dimension(:), allocatable :: sort_struct_double_temp
        
        
        !----------------------------------------------------------------------------
        
        !============================================================================
        ! Body of the subroutine
        !============================================================================
        
        ! combine all epologs in the new structure
        
        ! get number of epologs for the session
        ! note: size() without dimension specification returns the total number total elements in the array
        nr_epologs= size(epolog)
        
        ! determine the sum of observations (including possible duplicates) and delays in all epologs
        ! note: "nr_obs" is the total number of observations in one epolog. This number must be the same for the number of delays,
        !       although "nr_obs" has not been used in the ray-tracing part as the number of observations per station has been used,
        !       which has been indirectly derived from "nr_obs".
        !       sum() returns the sum of total elements in an array if no dimension is specified
        nr_obs_total= sum(epolog(1:nr_epologs) % nr_obs)
        
        ! allocate the substructures in structures "rdlog_epoch" and "temp_rdlog_epoch"
        ! note: a temporal variable is used as due to possible stackoverflow at the later sorting stage for
        !       large datasets, the implicit re-assignment to the variable itself when the sorting order is established is not possible
        allocate( rdlog_epoch % observations(nr_obs_total), &
                  rdlog_epoch % delay(nr_obs_total), &
                  rdlog_epoch % meteo_stat_out(nr_obs_total), &
                  rdlog_epoch % grib_epoch_mjd(nr_obs_total) )
        allocate( temp_rdlog_epoch % observations(nr_obs_total), &
                  temp_rdlog_epoch % delay(nr_obs_total), &
                  temp_rdlog_epoch % meteo_stat_out(nr_obs_total), &
                  temp_rdlog_epoch % grib_epoch_mjd(nr_obs_total) )
        
        
        ! combine all epologs
        
        ! initialize variable for index of end for added data
        ! note: setting "ind_end" to 0 at begin of loop is necessary as it is used when determining "ind_begin"
        ind_end= 0
        
        ! loop over all epologs
        ! combine all observations and the respective delays
        do curr_epolog_ind= 1, nr_epologs
            
            ! determine beginning index for adding data
            ind_begin= ind_end + 1
            
            ! determine index of end of added data
            ! note: use the number of observations in the current epolog to define the end index
            ind_end= ind_end + epolog(curr_epolog_ind) % nr_obs
            
            ! add the data of the current epolog to the combined structure "rdlog_epoch"
            temp_rdlog_epoch % observations(ind_begin:ind_end)= epolog(curr_epolog_ind) % observations
            temp_rdlog_epoch % delay(ind_begin:ind_end)= epolog(curr_epolog_ind) % delay
            temp_rdlog_epoch % meteo_stat_out(ind_begin:ind_end)= epolog(curr_epolog_ind) % meteo_stat_out
            temp_rdlog_epoch % grib_epoch_mjd(ind_begin:ind_end)= epolog(curr_epolog_ind) % epoch % mjd
            
        end do
        
        
        !------------------------------------------------------------------------------------------
        
        ! sort the combined data
        ! Sorting of the data first after the source name, then the station name and then the mjd using a
        ! stable sorting algorithm like sort() leads to an mjd sorted output where duplicates of the same
        ! observation are consecutively ordered.
        ! Note: Sorting the data in case of epolog mode "one_epoch_per_obs" is in principle not necessary,
        !       but output of results would then not be ordered!
        
        
        ! 1.)
        ! sort all data alphabetically using the source name
        
        ! allocate the structures used for sorting (see module_msort for help)
        allocate(sort_struct_char(nr_obs_total))
        allocate(sort_struct_char_temp((nr_obs_total+1)/2))
    
        ! assign source names and indices to the structure for sorting
        do i= 1, nr_obs_total
            sort_struct_char(i) % value= temp_rdlog_epoch % observations(i) % source
            sort_struct_char(i) % order= i
        end do
    
        ! sort the list with the source names
        ! use mergesort for sorting (note: original order is preserved in case of equal entries; mergesort is a stable sorting algorithm)
        call msort(sort_struct_char, nr_obs_total, sort_struct_char_temp)
    
        ! sort all variables in the structure according to the new alphabetical order
        ! note: also the "source"-field needs to be sorted as it has not been done yet
        
        ! note: The temporal variable must be used for the assignments in the new order here anyway, but in case of large datasets
        !       an assignment of a= a(order) would lead to a stack overflow.
        rdlog_epoch % observations = temp_rdlog_epoch % observations(sort_struct_char % order)
        rdlog_epoch % delay = temp_rdlog_epoch % delay(sort_struct_char % order)
        rdlog_epoch % meteo_stat_out = temp_rdlog_epoch % meteo_stat_out(sort_struct_char % order)
        rdlog_epoch % grib_epoch_mjd = temp_rdlog_epoch % grib_epoch_mjd(sort_struct_char % order)
        
        ! deacllocate the structure used for sorting after station names
        deallocate(sort_struct_char)
        deallocate(sort_struct_char_temp)
        
        
        ! 2.)
        ! sort all data alphabetically using the station name
        ! note: Data becomes sorted alphabetically after the station names, previous sorting
        !       order is preserved if possible as mergesort algorithm is used for the sorting and this is a stable sorting
        !       algorithm.
        
        ! allocate the structures used for sorting (see module_msort for help)
        allocate(sort_struct_char(nr_obs_total))
        allocate(sort_struct_char_temp((nr_obs_total+1)/2))
    
        ! assign station names and indices to the structure for sorting
        do i= 1, nr_obs_total
            sort_struct_char(i) % value= temp_rdlog_epoch % observations(i) % station
            sort_struct_char(i) % order= i
        end do
    
        ! sort the list with the station names
        ! use mergesort for sorting (note: original order is preserved in case of equal entries; mergesort is a stable sorting algorithm)
        call msort(sort_struct_char, nr_obs_total, sort_struct_char_temp)
    
        ! sort all variables in the structure according to the new alphabetical order
        ! note: also the "station"-field needs to be sorted as it has not been done yet
        
        ! note: The temporal variable must be used for the assignments in the new order here anyway, but in case of large datasets
        !       an assignment of a= a(order) would lead to a stack overflow.
        rdlog_epoch % observations = temp_rdlog_epoch % observations(sort_struct_char % order)
        rdlog_epoch % delay = temp_rdlog_epoch % delay(sort_struct_char % order)
        rdlog_epoch % meteo_stat_out = temp_rdlog_epoch % meteo_stat_out(sort_struct_char % order)
        rdlog_epoch % grib_epoch_mjd = temp_rdlog_epoch % grib_epoch_mjd(sort_struct_char % order)
        
        ! deacllocate the structure used for sorting after station names
        deallocate(sort_struct_char)
        deallocate(sort_struct_char_temp)
    
        
        ! 3.)
        ! sort all data ascending after mjd
        ! note: Data becomes sorted chronologically after the mjd in ascending order, previous sorting
        !       order is preserved if possible as mergesort algorithm is used for the sorting and this is a stable sorting
        !       algorithm.
        
        ! assign the currently just alphabetically ordered data to the temporal variable
        ! note: this is necessary as the old content of the temporal variable is not yet sorted alphabetically
        temp_rdlog_epoch= rdlog_epoch
        
        
        ! allocate the structures used for sorting (see module_msort for help)
        allocate(sort_struct_double(nr_obs_total))
        allocate(sort_struct_double_temp((nr_obs_total+1)/2))
    
        ! assign station names and indices to the structure for sorting
        do i= 1, nr_obs_total
            sort_struct_double(i) % value= temp_rdlog_epoch % observations(i) % mjd
            sort_struct_double(i) % order= i
        end do
    
        ! sort the list with the station names
        ! use mergesort for sorting (note: original order is preserved in case of equal entries; mergesort is a stable sorting algorithm)
        call msort(sort_struct_double, nr_obs_total, sort_struct_double_temp)
    
        ! sort all variables in the structure according to the new alphabetical order
        ! note: also the "mjd"-field needs to be sorted as it has not been done yet
        
        ! note: The temporal variable must be used for the assignments in the new order here anyway, but in case of large datasets
        !       an assignment of a= a(order) would lead to a stack overflow.
        rdlog_epoch % observations = temp_rdlog_epoch % observations(sort_struct_double % order)
        rdlog_epoch % delay = temp_rdlog_epoch % delay(sort_struct_double % order)
        rdlog_epoch % meteo_stat_out = temp_rdlog_epoch % meteo_stat_out(sort_struct_double % order)
        rdlog_epoch % grib_epoch_mjd = temp_rdlog_epoch % grib_epoch_mjd(sort_struct_double % order)
    
        ! deacllocate the structure used for sorting after mjd
        deallocate(sort_struct_double)
        deallocate(sort_struct_double_temp)
        
        
        ! add the session name to the output structure
        ! note: use the first epolog to get the session name as it is equal for all epologs
        rdlog_epoch % session_name= epolog(1) % session_name
        
    
    end subroutine combine_and_sort_rd
    
end module module_combine_and_sort_rd    