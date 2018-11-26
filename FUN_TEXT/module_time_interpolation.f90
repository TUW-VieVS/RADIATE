! module_time_interpolation.f90
    

!****************************************************************************
!
! PURPOSE:  Module containing subroutine to prepare input observations and delays for time domain interpolation of duplicate observations
!           with delays at different grib-epochs. If there are two or more same observations for the same mjd (of the observation) for
!           one specific station then a time interpolation of the delays at different grib-epochs will be done using an external function.
!           
!           Note: For VLBI: Due to erroneous NGS- and therefore also wrong Azel-files, the check for duplicate observations is enhanced by a check
!                 of the same observed source as there sometimes occur errors in the NGS-file, where a specific station observes two different
!                 sources at the same mjd, which is in principle not possible as long it is not a twin telescope which shares its
!                 name with the other one.
!                 Furthermore in case of GNSS observations this check is needed as there are always more than
!                 one "sources" = satellites observed at the same time by one station!
!                 Also a check of identical azimuth and elevation for duplicate observations is used for safety reasons in case (artificial) observations
!                 are present that suggest that a specific source is observed at the same mjd from the same station under different azimuth and/or
!                 elevation angles.
!
!           Note: Created as module to avoid explicit inferface for the subroutine as needed because of assumed-shape dummy arguments in subroutine.
!
!
! INPUT and OUTPUT:
!        rdlog_epoch... structure containing all observations including duplicates on the input,
!                       but is reduced by the duplicates on the output (time domain interpolation of delays already done)
!                           Structure has the following variables and substructures:
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
! 23.04.2015: create the Fortran-file based on the Matlab-file "time_interpolation.m"
! 27.04.2015: programming
! 30.04.2015: programming
! 07.05.2015: correct comments
! 13.05.2015: correct comments
! 25.06.2015: add check of same source name for duplicate finding, needed in VLBI in case
!             that (impossible but happening in NGS and Azel files ) a specific station
!             observes two different sources at the same time, but also generally needed
!             for GNSS observations;
!             add check of same azimuth and elevation for duplicate finding (in principle
!             only necessary for artificial observations as the same observed source at one
!             specific mjd should not have been observed under different azimuth and/or
!             elevation angles)
! 02.12.2015: change the direct time interpolation of the mapping factors to a new
!             calculation based on the time interpolated slant and zenith delays
! 17.12.2015: add the total mapping factor calculation
!             changes due due to adding of the meteorological data at the station position from the NWM
! 11.01.2016: correct comments due to adding of the water vapour pressure from the azel-file to the "observations" substructure
! 12.01.2016: add comments
!
!****************************************************************************

module module_time_interpolation

contains
    
    subroutine time_interpolation( rdlog_epoch )
    
        ! Define modules to be used
        use module_type_definitions, only : rdlog_epoch_type, duplicate_obs_data_type
        use module_lagrange_int
    
        !----------------------------------------------------------------------------
    
        ! Variable definitions
        implicit none
    
    
        ! dummy arguments
        !----------------
    
        ! INPUT and OUTPUT
    
        ! define input and output variable for storing observations and their delays
        ! note: on input the structure contains duplicate observations with delays at different grib-epochs
        !       in output the structure contains only unique observations with time interpolated delays
        type(rdlog_epoch_type), intent(in out) :: rdlog_epoch
    
    
        ! local variables
        !----------------
        
        ! define variable for storing the total number of observations (including possible duplicates) and delays in all epologs
        integer :: nr_obs_total
        
        ! define temporal variable for storing observations and their delays
        ! note: as the input is the output variable, but with reduced number of observations as duplicates have been interpolated
        !       a temporal variable is needed to do the processing steps during the subroutine
        type(rdlog_epoch_type) :: temp_rdlog_epoch
        
        ! define variable for storing the duplicate observation check results
        logical, dimension(:), allocatable :: duplicate_check
        
        ! define variable for storing the check if specific observations have already been treated by the duplicate check
        logical, dimension(:), allocatable :: already_treated_obs
        
        ! define variable for storing the number of untreated appearances of the current observation
        integer :: nr_obs_untreated_found
        
        ! define variable for loop index
        integer :: ind_obs
        
        ! define variable for index of unique observations
        integer :: ind_unique
        
        ! define structure for assigning the values needed for time domain interpolation of appearances of the same observation.
        type(duplicate_obs_data_type) :: duplicate_obs
        
        ! define variable for storing the Lagrange basis polynomials
        double precision, dimension(:), allocatable :: L
        
        !----------------------------------------------------------------------------
        
        !============================================================================
        ! Body of the subroutine
        !============================================================================

        ! check for multiple delays of the same observation
        
        ! compare mjd, station name, source name, azimuth and elevation of observations
        ! note: Due to the sorting after the function "combine_and_sort_rd", where all observations are sorted alphabetically after the source name,
        !       then the station name and then after the mjd, duplicates of the same observation with ray-traced delays at different grib-epochs should
        !       be sorted consecutively.
        !       This means that all duplicates should follow each other consecutively and should not be separated
        !       by other observations that are not duplicates of the currently looked at observation.
        !       Nevertheless the search of duplicates always considers all observations when searching, so
        !       there will be no error in case of not consecutively sorted duplicates!
        !       An exceptional case where the input observations are not already sorted in a way that all duplicates of
        !       one observation are following consecutively may occur if there are two or more observations that have the
        !       same mjd, station and source name, but azimuth and/or elevation are different though it is the same source!
        
        
        ! get number of observations (including the duplicates)
        ! note: size() without dimension specification returns the total number total elements in the array
        nr_obs_total= size(rdlog_epoch % observations)
        
        ! allocate the temporal variable substructures according to the number of observations in the input structure
        allocate( temp_rdlog_epoch % observations(nr_obs_total), &
                  temp_rdlog_epoch % delay(nr_obs_total), &
                  temp_rdlog_epoch % meteo_stat_out(nr_obs_total), &
                  temp_rdlog_epoch % grib_epoch_mjd(nr_obs_total) )
        
        ! assing the input structure content to the temporal structure
        temp_rdlog_epoch= rdlog_epoch
        
        ! deallocate the substructures of the input and output structure
        ! note: these will be allocated again in a later step when their size without duplicates is known
        deallocate( rdlog_epoch % observations, &
                    rdlog_epoch % delay, &
                    rdlog_epoch % meteo_stat_out, &
                    rdlog_epoch % grib_epoch_mjd )
        
        ! allocate variable for the duplicate observation check (logical vector of same size as nr_obs_total)
        allocate( duplicate_check(nr_obs_total) )
        
        ! allocate variable for check if specific observations have already been treated by the duplicate check (logical vector of same size as nr_obs_total)
        allocate( already_treated_obs(nr_obs_total) )
        ! initialize "already_treated_obs" with .FALSE.
        already_treated_obs= .FALSE.
        
        ! initialize index variable for unique observations
        ind_unique= 0
        
        ! loop over all observations to find duplicates of observations and calculate time interpolation for
        ! them
        do ind_obs= 1, nr_obs_total
            
            ! determine the duplicates of the current (ind_obs) observation
            ! check which observations have the same mjd (of the observation) and station name as the current (ind_obs) observation
            ! note: assignments from a previous loop cycle will be overwritten
            ! attention: The check for the same observed source is necessary in case that a specific station observes two different sources at the same time!
            !            In principle, this is not possible in VLBI as long there are not twin telescopes with the same station name, but neveretheless due to
            !            errors in the NGS-files this phenomenon may occur and be introduced to the Azel-files.
            !            For GNSS observations this check is always necessary as there are always more than
            !            one "sources" = satellites observed at the same time by one station!
            !            Not checking for the source name in the duplicate finding would lead then to found duplicates that consist of different observations for
            !            the same station at the same time and the time interpolation would create a mixture of totally different delays as the source has not been
            !            the same!
            !            The check for equality of azimuth and elevation is for safety reasons generally in case of artificial observations where the same named
            !            source is observed under different azimuth and/or elevation angles!
            duplicate_check= ( ( temp_rdlog_epoch % observations(ind_obs) % mjd == temp_rdlog_epoch % observations(:) % mjd ) .AND. &
                               ( temp_rdlog_epoch % observations(ind_obs) % station == temp_rdlog_epoch % observations(:) % station ) .AND. &
                               ( temp_rdlog_epoch % observations(ind_obs) % source == temp_rdlog_epoch % observations(:) % source ) .AND. &
                               ( temp_rdlog_epoch % observations(ind_obs) % az == temp_rdlog_epoch % observations(:) % az ) .AND. &
                               ( temp_rdlog_epoch % observations(ind_obs) % elev == temp_rdlog_epoch % observations(:) % elev ) )
            
            ! reduce the found duplicates by the already treated duplicates
            ! note: this leads only to .TRUE. if "duplicate_check" was .TRUE. "already_treated_obs" was .FALSE.
            duplicate_check= ( ( duplicate_check .AND. ( .NOT. already_treated_obs ) ) )
            
            ! update "already_treated_obs"
            already_treated_obs= ( already_treated_obs .OR. duplicate_check )
            
            ! determine the number of untreated appearances of the current observation
            nr_obs_untreated_found= count(duplicate_check)
            
            
            ! check the number of appearances of the current observation and if it is more than once do time domain interpolation
            
            ! in case that all appearances of the current observation have already been treated, so "nr_obs_untreated_found" is 0
            ! note: testing this case may speed up the code as otherwise two remaining cases would be tested
            if (nr_obs_untreated_found == 0) then
                
                ! cycle to next loop step as processing of all appearances for the current observation has been done
                cycle
            
            ! in case the current observation is appearing only once (there should be no other appearances that have already been treated)
            ! then directly assign the unique observation as no time domain interpolation is needed
            else if (nr_obs_untreated_found == 1) then
                
                ! raise the index for unique observations
                ind_unique= ind_unique + 1
                
                ! assign the only found observation and all related data to the temporal structure at the next index for unique observations
                ! note: This overwrites the old structure content of the temporal variable, but this is desired and not problematic
                !       as "ind_unique" can at the max be the same as "ind_obs"! So no values of up to now unchecked entries for duplicates will be overwritten!
                ! attention: Also the logical check should not be affected by this change of values as there is the variable "already_treated_obs" that would prevent
                !            undesired re-use of the same observation!
                ! note: The index of the single found observation must be "ind_obs".
                temp_rdlog_epoch % observations(ind_unique)= temp_rdlog_epoch % observations(ind_obs)
                temp_rdlog_epoch % delay(ind_unique)= temp_rdlog_epoch % delay(ind_obs)
                temp_rdlog_epoch % meteo_stat_out(ind_unique)= temp_rdlog_epoch % meteo_stat_out(ind_obs)
                temp_rdlog_epoch % grib_epoch_mjd(ind_unique)= temp_rdlog_epoch % grib_epoch_mjd(ind_obs)
                
            ! in case the current observation is appearing more than once
            ! then time domain interpolation is needed
            else if (nr_obs_untreated_found > 1) then
                
                ! raise the index for unique observations
                ind_unique= ind_unique + 1
                
                ! Allocate the substructures for assigning the values needed for time domain interpolation of
                ! appearances of the same observation. Size of substructures is according to the number of same observations found.
                ! note: duplicate_obs % mjd_obs must not be allocated as it is defined as a scalar
                allocate( duplicate_obs % delay(nr_obs_untreated_found), &
                          duplicate_obs % meteo_stat_out(nr_obs_untreated_found), &
                          duplicate_obs % grib_epoch_mjd(nr_obs_untreated_found) )
                               
                ! assign the duplicates (all appearances of the same observation) to the structure used for time domain interpolation
                ! note: the mask of duplicate_check leads to assigning only the .TRUE. elements to the result
                duplicate_obs % delay= pack(temp_rdlog_epoch % delay, duplicate_check)
                duplicate_obs % meteo_stat_out= pack(temp_rdlog_epoch % meteo_stat_out, duplicate_check)
                duplicate_obs % grib_epoch_mjd= pack(temp_rdlog_epoch % grib_epoch_mjd, duplicate_check)
                ! note: It is sufficient to assign the mjd of one appearance of the current checked observation as all duplicates
                !       should have the same mjd as it is one of the parameters in the check of same observations!
                !       Therefore the index "ind_obs" can be used for the assignment as this mjd was part of the check and must be true for all other appearances!
                duplicate_obs % mjd_obs= temp_rdlog_epoch % observations(ind_obs) % mjd
                
                ! allocate variable for Lagrange basis polynomials for interpolation at x_int= duplicate_obs % mjd_obs
                ! note: the size of L is determined through the number of duplicates = nr_obs_untreated_found
                allocate( L(nr_obs_untreated_found) )
                
                ! call subroutine to do time domain interpolation of found duplicates for one observation using the
                ! different delays from the different grib-file epochs
                ! note: The call of the lagrange interpolation is done for only one variable (dz_total) of the delay-structure and the interpolated value is
                !       assigned to the "temp_rdlog_epoch % delay"-structure at the next index for unique observations ("ind_unique").
                !       This overwrites the old structure content of the temporal variable, but this is desired and not problematic as "ind_unique" can at
                !       the max be the same as "ind_obs"! So no values of up to now unchecked entries for duplicates will be overwritten!
                !       The output of the Lagrange basis polynomials is then later used to do the interpolation for all remaining values in the delay-structure.
                ! attention: "duplicate_obs % mjd_obs" is a scalar
                call lagrange_int( x_int= duplicate_obs % mjd_obs, &
                                   x= duplicate_obs % grib_epoch_mjd, &
                                   y= duplicate_obs % delay % dz_total, &
                                   y_int= temp_rdlog_epoch % delay(ind_unique) % dz_total, &
                                   L_out= L ) ! note: output of optional Lagrange basis polynomials in order to calculate interpolation for other vales in "% delay" and in "meteo_stat_out"
                
                ! interpolate the remaining values in the delay-structure using the above determined Lagrange basis polynomials
                ! note: The interpolated values are assigned to the "temp_rdlog_epoch % delay"-structure at the next index for unique observations ("ind_unique").
                !       This overwrites the old structure content of the temporal variable, but this is desired and not problematic as "ind_unique" can at
                !       the max be the same as "ind_obs"! So no values of up to now unchecked entries for duplicates will be overwritten!
                ! note: dot_product() does the following: sum(vector_a * vector_b), where vector_a * vector_b is the element-wise multiplication
                temp_rdlog_epoch % delay(ind_unique) % dz_h= dot_product(L, duplicate_obs % delay % dz_h)
                temp_rdlog_epoch % delay(ind_unique) % dz_w= dot_product(L, duplicate_obs % delay % dz_w)
                temp_rdlog_epoch % delay(ind_unique) % ds_total_geom= dot_product(L, duplicate_obs % delay % ds_total_geom)
                temp_rdlog_epoch % delay(ind_unique) % ds_total= dot_product(L, duplicate_obs % delay % ds_total)
                temp_rdlog_epoch % delay(ind_unique) % ds_h_geom= dot_product(L, duplicate_obs % delay % ds_h_geom)
                temp_rdlog_epoch % delay(ind_unique) % ds_h= dot_product(L, duplicate_obs % delay % ds_h)
                temp_rdlog_epoch % delay(ind_unique) % ds_w= dot_product(L, duplicate_obs % delay % ds_w)
                temp_rdlog_epoch % delay(ind_unique) % e_stat= dot_product(L, duplicate_obs % delay % e_stat)
                temp_rdlog_epoch % delay(ind_unique) % e_outgoing_rt= dot_product(L, duplicate_obs % delay % e_outgoing_rt)
                temp_rdlog_epoch % delay(ind_unique) % dgeo= dot_product(L, duplicate_obs % delay % dgeo)
                temp_rdlog_epoch % delay(ind_unique) % diff_e= dot_product(L, duplicate_obs % delay % diff_e)
                
                ! calculate the "interpolated" values of the total, hydrostatic and wet mapping factor
                ! note: Due to numeric resons a direct linear interpolation of the mapping factors may lead to significant differences.
                !       Therefore the mapping factors are newly calculated based on the interpolated slant and zenith delays.
                temp_rdlog_epoch % delay(ind_unique) % mf_total_geom= temp_rdlog_epoch % delay(ind_unique) % ds_total_geom / temp_rdlog_epoch % delay(ind_unique) % dz_total
                temp_rdlog_epoch % delay(ind_unique) % mf_h_geom= temp_rdlog_epoch % delay(ind_unique) % ds_h_geom / temp_rdlog_epoch % delay(ind_unique) % dz_h
                temp_rdlog_epoch % delay(ind_unique) % mf_w= temp_rdlog_epoch % delay(ind_unique) % ds_w / temp_rdlog_epoch % delay(ind_unique) % dz_w
                
                
                ! interpolate the values in the meteo_stat_out-structure using the above determined Lagrange basis polynomials
                ! note: The interpolated values are assigned to the "temp_rdlog_epoch % meteo_stat_out"-structure at the next index for unique observations ("ind_unique").
                !       This overwrites the old structure content of the temporal variable, but this is desired and not problematic as "ind_unique" can at
                !       the max be the same as "ind_obs"! So no values of up to now unchecked entries for duplicates will be overwritten!
                ! note: dot_product() does the following: sum(vector_a * vector_b), where vector_a * vector_b is the element-wise multiplication
                temp_rdlog_epoch % meteo_stat_out(ind_unique) % p= dot_product(L, duplicate_obs % meteo_stat_out % p)
                temp_rdlog_epoch % meteo_stat_out(ind_unique) % T= dot_product(L, duplicate_obs % meteo_stat_out % T)
                temp_rdlog_epoch % meteo_stat_out(ind_unique) % wvpr= dot_product(L, duplicate_obs % meteo_stat_out % wvpr)
                
                
                ! determination of "interpolated" "grib_epoch_mjd"-value
                ! note: As time domain interpolation was carried out the new value for "grib_epoch_mjd" is determined as the mean of the duplicate entries.
                temp_rdlog_epoch % grib_epoch_mjd(ind_unique)= sum(duplicate_obs % grib_epoch_mjd) / nr_obs_untreated_found
                
                ! determination of "interpolated" logicals "break_elev" and "break_layer" in "delay"
                ! note: In order to retain possible breaks from ray-tracing function any() is used to determine if the new time domain interpolated delay values
                !       result from using any delay values where a break occured in the calculation of the delay.
                !       Therefore .TRUE. is received if any of the duplicates was .TRUE.
                temp_rdlog_epoch % delay(ind_unique) % break_elev= any( duplicate_obs % delay % break_elev )
                temp_rdlog_epoch % delay(ind_unique) % break_layer= any( duplicate_obs % delay % break_layer )
                
                ! assign the observational data of the current (ind_obs) observation to the temporal structure at the next index for unique observations
                ! note: This overwrites the old structure content of the temporal variable, but this is desired and not problematic
                !       as "ind_unique" can at the max be the same as "ind_obs"! So no values of up to now unchecked entries for duplicates will be overwritten!
                ! attention: Also the logical check should not be affected by this change of values as there is the variable "already_treated_obs" that would prevent
                !            undesired re-use of the same observation!
                ! note: The index of the first duplicate of the current observation must be "ind_obs".
                temp_rdlog_epoch % observations(ind_unique)= temp_rdlog_epoch % observations(ind_obs)
                
                ! deallocate the substructures of the duplicates
                ! note: "duplicate_obs % mjd_obs" can not be deallocated as it is a not allocatable scalar.
                !       Its value will be overwritten when variable is used the next time.
                deallocate( duplicate_obs % delay, &
                            duplicate_obs % meteo_stat_out, &
                            duplicate_obs % grib_epoch_mjd )
                
                ! deallocate variable for Lagrange basis polynomials for interpolation at x_int= duplicate_obs % mjd_obs
                deallocate( L )
                
            end if
            
        end do
        
        ! allocate the output substructures using "ind_unique", which determines the size of the substructures
        allocate( rdlog_epoch % observations(ind_unique), &
                  rdlog_epoch % delay(ind_unique), &
                  rdlog_epoch % meteo_stat_out(ind_unique), &
                  rdlog_epoch % grib_epoch_mjd(ind_unique) )
        
        ! assign the values from "temp_rdlog_epoch" to "rdlog_epoch" using only the really unique observations (1: ind_unique)
        rdlog_epoch % observations= temp_rdlog_epoch % observations(1: ind_unique)
        rdlog_epoch % delay= temp_rdlog_epoch % delay(1: ind_unique)
        rdlog_epoch % meteo_stat_out= temp_rdlog_epoch % meteo_stat_out(1: ind_unique)
        rdlog_epoch % grib_epoch_mjd= temp_rdlog_epoch % grib_epoch_mjd(1: ind_unique)
        
        ! note: deallocation of "temp_rdlog_epoch" is done automatically at the exit of this subroutine
        
    end subroutine time_interpolation
    
end module module_time_interpolation