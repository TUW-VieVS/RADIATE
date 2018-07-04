! module_get_observing_stations.f90

!****************************************************************************
!
! PURPOSE:  Module containing subroutine to determine the unique station names of stations carrying out observations.
!            
!           "get_observing_stations" is a subroutine to determine from all observations, which stations have
!           been observing. Additionally the number of observations per station is determined.
!
!           Attention: Input of station names from the observations must be sorted alphabetically.
!                      This means that the observations in the epolog have to be sorted by station names,
!                      so that the same station names are grouped (eg. BADARY, BADARY, ..., BADARY, FORTLEZA, ..., FORTLEZA).
!                      Otherwise the determination of the unique station names will be corrupted and also
!                      the following ray-tracing will be corrupted as wrong station information may be used for
!                      ray-tracing as the assignment of station values and observation values may not
!                      correspond for the same actually observing station.
!
!           Note: Alphabetically sorted input list leads to alphabetically sorted output of the observing
!                      stations (unique station names).
!
!           Note: Created as module to avoid explicit inferface for the subroutine as needed because of assumed-shape dummy argument in subroutine.
! 
!
! 
! INPUT:
!         nr_obs................... number of observations
!         stations_from_obs........ vector containing names of stations for all observations,
!                                   contains duplicates according to observations
!                                   attention: has to be already sorted alphabetically!
! 
! 
! OUTPUT:
!         obs_stations........... vector containing names of all observing stations (no duplicates)
!         nr_obs_per_obsstat..... vector containing number of observations for each observing station
!         nr_observing_stations.. number of different observing stations
!
!----------------------------------------------------------------------------
! 
! Fortran-file created by Armin Hofmeister
!   based on Matlab scripts created by Armin Hofmeister
! 
!----------------------------------------------------------------------------
! History:
! 
! 26.11.2014: create the Fortran-file based on the Matlab-file "load_session_index.m"
! 27.11.2014: programming
! 01.12.2014: comments
! 26.01.2015: comments
! 29.01.2015: change to module as to save explicit interface block
! 12.02.2015: use constant to define length of station names
! 21.04.2015: comments
! 11.05.2015: add comments
! 12.05.2015: programming, add comments
! 13.05.2015: comments
! 14.09.2015: enhance program code by using nr_obs as input to determine dimensions of variables
!             in the declaration part
!
!****************************************************************************

module module_get_observing_stations

contains
    
    subroutine get_observing_stations( nr_obs, &
                                       stations_from_obs, &
                                       obs_stations, &
                                       nr_obs_per_obsstat, &
                                       nr_observing_stations )

        ! Define modules to be used
    
        use module_constants, only: len_statname
    
        !----------------------------------------------------------------------------
    
        ! Variable definitions
        implicit none
    
    
        ! dummy arguments
        !----------------
    
        ! INPUT
        
        ! define variable for the input of the number of entries in stations_from_obs (number of observations)
        integer, intent(in) :: nr_obs
        
        ! define variable for the input of station names from the observations (containing possible duplicates)
        ! note: see "module_constants" for length of station names
        character(len= len_statname), dimension(nr_obs), intent(in) :: stations_from_obs
    
        
        ! OUTPUT
        
        ! define variable for the unique station names
        ! note: see "module_constants" for length of station names
        character(len= len_statname), dimension(:), allocatable, intent(out) :: obs_stations
    
        ! define vector for the number of observations per (unique) station
        integer, dimension(:), allocatable, intent(out) :: nr_obs_per_obsstat
    
        ! define variable for storing the number of different observing stations
        integer, intent(out) :: nr_observing_stations
    
    
    
        ! local variables
        !----------------
    
        ! define variable for index of current observation that is looked at
        integer :: ind_obs
    
        ! define logical vector for storing comparison results of station names
        ! note: size is the number of station names from the observations including possible duplicates
        logical, dimension(nr_obs) :: obs_for_curr_stat
    
        ! define temporal vector for the unique station names
        ! note: see "module_constants" for length of station names
        ! note: size is the maximum possibe needed size which is the number of station names from the observations including possible duplicates
        ! note: This variable is necessary as the found unique stations can not be written to the input variable "stations_from_obs" because of
        !       the intent(in) attribute that prohibits new definition of input variables. Otherwise an overwriting would be possible as this
        !       would not disturb the searching algorithm for the unique stations.
        character(len= len_statname), dimension(nr_obs) :: temp_obs_stations
    
        ! define temporal vector for the number of observations per unique station name
        ! note: size is the maximum possibe needed size which is the number of station names from the observations including possible duplicates
        integer, dimension(nr_obs) :: temp_nr_obs_per_obsstat
    
        ! define variable for storing the index of the observing (unique) station name
        integer :: ind_obsstat
    
        
        ! CONSTANTS
        
        ! see "module_constants" for "len_statname"
        
    
        !----------------------------------------------------------------------------
    
        !============================================================================
        ! Body of the subroutine
        !============================================================================
    
    
        ! Get all observing stations (unique station names)
    
        ! set index for the current observation entry to 1
        ind_obs= 1
    
        ! set index for observing stations (start at 0!!!)
        ind_obsstat= 0
    
    
        ! find names of all observing stations and determine the number of observations per station
        ! loop over the observations (input list of station names with duplicates needs to be sorted alphabetically)
        do while (ind_obs <= nr_obs) ! do while the input list has not reached the end
            ! compare first entry of stations_from_obs to all entries in stations_from_obs and receive
            ! logical array with .TRUE. at equal names
            obs_for_curr_stat= (stations_from_obs(ind_obs) == stations_from_obs)
        
            ! raise index for observing station
            ind_obsstat=ind_obsstat + 1
        
            ! assign name of current compared station as new found observing station
            temp_obs_stations(ind_obsstat)= stations_from_obs(ind_obs)
            
            ! determine number of observations for this station
            ! note: count determines the number of true elements
            temp_nr_obs_per_obsstat(ind_obsstat)= count(obs_for_curr_stat)
        
            ! determine index of next station in duplicate list that is not equal to the current station
            ! note: in order to do this correct the input list of station names from the observations has to be sorted alphabetically as stated before
            !       the next entry is the one after the last .TRUE. entry (additive list of last ind_obs value + new .TRUE. entries)
            ind_obs= ind_obs + temp_nr_obs_per_obsstat(ind_obsstat)
        
        end do
    
        ! get the number of different observing stations
        nr_observing_stations= ind_obsstat
    
        ! allocate the output vectors according to the found number of different observing stations
        allocate(obs_stations(nr_observing_stations))
        allocate(nr_obs_per_obsstat(nr_observing_stations))
    
        ! assign the unique station names and their observation numbers
        ! note: use the number of found unique stations for getting the upper boundary
        obs_stations= temp_obs_stations(1:nr_observing_stations)
        nr_obs_per_obsstat= temp_nr_obs_per_obsstat(1:nr_observing_stations)
    
    
    end subroutine get_observing_stations
    
end module module_get_observing_stations