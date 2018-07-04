! module_create_trp_global.f90
    

!****************************************************************************
!
! PURPOSE:  Module containing subroutine for the creation of the .trp-file.
!
!           This subroutine carries creates a .trp-file for the session containing all trp-format-specific data.
!           The output format is consistent with the format defined by the Goddard Space Flight Center (GSFC) except for
!           some additional comment lines.
!           The sorting order of the station data is alphabetically and the sorting order of the observations is by mjd and secondly alphabetically.
!
!           Note: Created as module to avoid explicit inferface for the subroutine as needed because of assumed-shape dummy argument in subroutine.
!
!           Attention: The station names and sources names are always written to the trp-file using a width of 8 characters. In case of station and sources names
!                      that are longer, the names will be truncated.
!                      The length of 8 characters is defined by the trp-format.
!
!
! INPUT:
!        observing_stations_session... structure of type "observing_stations_type" containing the unique observing station data for the whole session
!                                      The structure contains the following variables:
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
!        rdlog_epoch... structure containing all observations and delays
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
!           % delay......... structure containing the delay data in the following variables:
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
!           % grib_epoch_mjd... note: deallocated at this stage:
!                                     prior used as vector containing the mjd of the grib-file epoch which has been used to calculate
!                                     the delay of each observation
!
!        save_path.............. specifies path to the saving directory for the .trp-file
!        RADIATE_version........ version info of the RADIATE program
!        RADIATE_subversion..... sub-version info of the RADIATE program
!        epolog_mode............ variable reports mode of epolog creation
!        interpolation_method... variable reports method of vertical interpolation
!        raytr_method........... variable reports method of determining the ray path
!        filename_stat_info..... filename of file with station data
!        epoch_name............. vector containing the epolog epoch names for each epoch = epolog
!
!        gridres................ structure containing for each epoch = epolog:
!           % res_lat................ latitude resolution in [°] of the grib-file used to calculate
!                                     meteorological parameters
!           % res_lon................ longitude resolution in [°] of the grib-file used to calculate
!                                     meteorological parameters
!
!----------------------------------------------------------------------------
! 
! Fortran-file created by Armin Hofmeister
!   based on Matlab scripts created by Armin Hofmeister
! 
!----------------------------------------------------------------------------
! History:
! 
! 12.05.2015: create the Fortran-file based on the Matlab-file "create_trp_global.m"
! 13.05.2015: programming
! 02.06.2015: programming
! 03.06.2015: programming
! 09.06.2015: programming
! 11.06.2015: programming
! 17.12.2015: correct comments due to adding of the total mapping factor calculation to the "delay" substructure
! 11.01.2016: correct comments due to adding of the water vapour pressure from the azel-file to the "observations" substructure
! 12.01.2016: add comments
! 28.01.2016: correct the writing of the observation seconds to the file for correct representation with possible needed leading zeros
!
!****************************************************************************

module module_create_trp_global

contains
    
    subroutine create_trp_global( observing_stations_session, &
                                  rdlog_epoch, &
                                  save_path, &
                                  RADIATE_version, &
                                  RADIATE_subversion, &
                                  epolog_mode, &
                                  interpolation_method, &
                                  raytr_method, &
                                  filename_stat_info, &
                                  epoch_name, &
                                  gridres )
    
        ! Define modules to be used
        use module_type_definitions, only: observing_stations_type, rdlog_epoch_type, gridres_type
        use module_constants, only: c, rad2deg
        use module_date_type_definition
        use module_mjd2date
        use module_ell2xyz
        
    
        !----------------------------------------------------------------------------
    
        ! Variable definitions
        implicit none
    
    
        ! dummy arguments
        !----------------
    
        ! INPUT
    
		! define input variable for storing the observing station data for the whole session
		type(observing_stations_type), dimension(:), intent(in) :: observing_stations_session
	
        ! define input variable for storing observations and their delays
        type(rdlog_epoch_type), intent(in) :: rdlog_epoch
        
        ! define input variable for the saving path
        character(len=*), intent(in) :: save_path
        
        ! define input variable for the version info of the RADIATE program
        character(len=*), intent(in) :: RADIATE_version
        
        ! define input variable for the sub-version info of the RADIATE program
        character(len=*), intent(in) :: RADIATE_subversion
        
        ! define input variable for the epolog mode
        character(len=*), intent(in) :: epolog_mode
        
        ! define input variable for the vertical interpolation method
        character(len=*), intent(in) :: interpolation_method
        
        ! define input variable for the ray-tracing method
        character(len=*), intent(in) :: raytr_method
        
        ! define input variable for the filename of file with station data
        character(len=*), intent(in) :: filename_stat_info
        
        ! define input variable for the epolog epoch names
        ! note: size is determined through input of each "epoch_name" variables from all epologs
        character(len=*), dimension(:), intent(in) :: epoch_name
        
        ! define input variable for the grid resolution (grib-file horizontal resolution)
        ! note: size is determined through input of all "gridres"-structures from the epologs and must be equal to "epoch_name"
        type(gridres_type), dimension(size(epoch_name)), intent(in) :: gridres
        
        
        ! local variables
        !----------------
        
        ! define variable for storing the total number of observations
        integer :: nr_obs
        
        ! define variable for storing the slant total delay (including geometric bending effect) in [s]
        double precision, dimension(:), allocatable :: slant_total_delay_s
        
        ! define variable for storing the zenith hydrostatic delay in [s]
        double precision, dimension(:), allocatable :: zenith_h_delay_s
        
        ! define variable for storing the zenith wet delay in [s]
        double precision, dimension(:), allocatable :: zenith_w_delay_s
        
        ! define variable for storing the azimuth in [°]
        double precision, dimension(:), allocatable :: azimuth_deg
        
        ! define variable for storing the elevation angle in [°]
        double precision, dimension(:), allocatable :: elevation_deg
        
        ! define struture for storing date of the observations from conversion of mjd
        type(date_type), dimension(:), allocatable :: date_obs
        
        ! define variable for storing the total number of observating stations in the session
        integer :: nr_stat
        
        ! define variable for storing the geocentric x-coordinates of the observing stations
        double precision, dimension(:), allocatable :: x_stat
        
        ! define variable for storing the geocentric y-coordinates of the observing stations
        double precision, dimension(:), allocatable :: y_stat
        
        ! define variable for storing the geocentric z-coordinates of the observing stations
        double precision, dimension(:), allocatable :: z_stat
        
        ! store the names of the observing stations to a new character variable with specific length
        ! dependent on the format of the .trp-file
        ! note: Set the character string length to 8 characters as it is the fixed length used to write the data to the .trp-file.
        character(len= 8), dimension(:), allocatable :: names_fixl_observing_stat
            
        ! store the names of the stations in the observations to a new character variable with specific length
        ! dependent on the format of the .trp-file
        ! note: Set the character string length to 8 characters as it is the fixed length used to write the data to the .trp-file.
        character(len= 8), dimension(:),allocatable :: names_fixl_observations_stat
            
        ! store the names of the sources in the observations to a new character variable with specific length
        ! dependent on the format of the .trp-file
        ! note: Set the character string length to 8 characters as it is the fixed length used to write the data to the .trp-file.
        character(len= 8), dimension(:),allocatable :: names_fixl_observations_sou
        
        ! store the seconds of each observation to a new character variable with specific length
        ! dependent on the format of the .trp-file
        ! note: Set the character string length to 4 characters as it is the fixed length used to write the data to the .trp-file.
        !       The string representation is necessary for the workaround that enables the output of the seconds with possible
        !       leading zeros.
        character(len= 4), dimension(:),allocatable :: obs_sec_str
        
        ! define loop variable
        integer :: ind_obs
        
        ! define format for writing the station header line
        character(len=:), allocatable :: header_stat_fmt
        
        ! define format for writing the observation header line
        character(len=:), allocatable :: header_obs_fmt
        
        ! define format for writing the station data lines
        character(len=:), allocatable :: data_stat_fmt
        
        ! define format for writing the observation data lines
        character(len=:), allocatable :: data_obs_fmt
        
        ! define variable for storing the current date
        character(len=8) :: date_value
        character(len=10) :: date_formatted
    
        ! define variable for storing the current time
        character(len=10) :: time_value
        character(len=12) :: time_formatted
        
        ! variable for storing the unit of the file
        integer, parameter :: file_unit_trp= 1
    
        ! variable for storing the opening status of the file
        integer :: open_status
        
        ! variable for storing the record of the file header line
        character(len=:), allocatable :: fileheader
        
        ! define loop variable
        integer :: ind_epoch
        
        ! define loop variable
        integer :: ind_stat
        
        
        ! CONSTANTS
        
        ! see "module_constants" for "c" and "rad2deg"
        
        
        !----------------------------------------------------------------------------
        
        !============================================================================
        ! Body of the subroutine
        !============================================================================
        

        ! Define coefficients
        
        ! define constant for speed of light
		! see "module_constants" for "c" in [m/s]
        
        ! transformation constant from [rad] to [deg]
        ! see "module_constants" for "rad2deg"
        
        !-----------------------------------------------------------------------------------------------
        
        ! Create .trp-file
        ! The .trp-file contains information and results of the ray-tracing processing for all stations and observations of the session
        ! in the specific trp-format.
        
		! get the current date and time
		! get the current date and time as a time stamp for the created trp-file
        call DATE_AND_TIME(date= date_value, time= time_value)
		
        ! format the date and time
        date_formatted= date_value(1:4) // '-' // date_value(5:6) // '-' // date_value(7:8)
        time_formatted= time_value(1:2) // ':' // time_value(3:4) // ':' // time_value(5:)
		
		! create new .trp-ascii-file for storing the data
        ! note: an old file is replaced by the new one as status is set to 'replace'
        open(unit= file_unit_trp, file= save_path // rdlog_epoch % session_name // '.trp', action= 'write', status= 'replace', iostat= open_status)
        
        ! Check if textfile can't be opened (iostat /= 0 --> error)
        if (open_status /= 0) then
            
            ! report warning message
            write(unit= *, fmt= '(tr4, a, a, a, /)') 'Warning: Problem with creating .trp-file: "', save_path // rdlog_epoch % session_name // '.trp"', '! No file created!'    
        
        else
            
            ! Calculate output parameters for the trp-file
            
            ! get the total number of observations (the same number as the number of delays) in the
            ! specific session = size of the combined (and time interpolated) epologs = rdlog_epoch
            ! note: use size() of substructure "% observations"
            nr_obs= size(rdlog_epoch % observations)
            
            ! allocate variables for data for the output to the trp-file, which are converted to other domains
            ! note: Storing the output data to new variables is necessary as the new definition of a variable that has the intent(in)
            !       attribute is not possible.
            allocate( slant_total_delay_s(nr_obs), &
                      zenith_h_delay_s(nr_obs), &
                      zenith_w_delay_s(nr_obs), &
                      azimuth_deg(nr_obs), &
                      elevation_deg(nr_obs), &
                      date_obs(nr_obs) )
            
            ! convert slant total delay from [m] to [sec] using speed of light
            slant_total_delay_s= rdlog_epoch % delay % ds_total_geom / c
            
            ! convert zenith hydrostatic delay from [m] to [sec] using speed of light
            zenith_h_delay_s= rdlog_epoch % delay % dz_h / c
            
            ! convert zenith wet delay from [m] to [sec] using speed of light
            zenith_w_delay_s= rdlog_epoch % delay % dz_w / c
            
            ! convert azimuth from [rad] to [°]
            azimuth_deg= rdlog_epoch % observations % az * rad2deg
            
            ! convert elevation from [rad] to [°]
            elevation_deg= rdlog_epoch % observations % elev * rad2deg
            
            
            ! convert mjd to date using function *mjd2date*
            ! note: only month and day are necessary, because they are not provided by the original AZEL-file and therefore not in the "observations"-structure
            ! note: don't use the seconds later, because they are not equal to the seconds provided by the AZEL-file (rdlog_epoch % observations % sec)!
            call mjd2date( rdlog_epoch % observations % mjd, &
                           date_obs )
            
            
            ! determine the geocentric x-, y-, z-coordinates of the observing stations for the S-records
            
            ! get the number of observing stations in the specific session
            nr_stat= size(observing_stations_session)
            
            ! allocate the variables for storing the geocentric coordinates
            ! note: size is according to the number of observing stations
            allocate(x_stat(nr_stat), &
                     y_stat(nr_stat), &
                     z_stat(nr_stat) )
            
            ! transform the geodetic to geocentric cartesian coordinates.
            ! call subroutine to do the transformation
            call ell2xyz( observing_stations_session % lat_ell, &
                          observing_stations_session % lon_ell, &
                          observing_stations_session % h_ell, &
                          x_stat, &
                          y_stat, &
                          z_stat )
            
            !-----------------------------------------------------------------------------------------------
            
            ! Define variables used to write data to the .trp-file
            ! Note: This section is necessary as input variables for observing station names, name of the station for each observation
            !       and the source name for each observation may have a character length that is different from the format length
            !       in the .trp-file, which is fixed to 8 characters for each of the mentioned variables.
            !       In case an input variable has a length shorter than the output format, then leading spaces would be introduced when
            !       writing the variable to the file. E.g. variable has a length of only 3 and the output format is a8 then 5 leading
            !       spaces would be introduced destroying the required left-justification of the character variables in the file.
            !       Therefore the variables are written to new variables with the same length that is used in the output for
            !       writing the file. In this case names that are shorter than the new length receive spaces at the end and not in front leading
            !       to left-justification as desired. Nevertheless input variables that are longer than the output format will be cut off at the last
            !       field position of the variable and will not be represented fully.
            ! Note: Also the seconds of the observation time need to be preprocessed since Fortran does not support the output of leading zeros in
            !       case of floating point variables. Therefore a workaround is needed.
            ! Attention: The splitting of the writing of the seconds into an integer part and a fractional part is not possible since rounding effects for a value of
            !            e.g. 46.99 would lead to the output of 461.0 instead of 47.0!
            
            ! allocate the variable for storing the names of the observing stations with specific length
            ! dependent on the format of the .trp-file
            allocate( names_fixl_observing_stat(nr_stat) )
            
            ! allocate the variables for storing the names of the stations and sources in the observations with specific length
            ! dependent on the format of the .trp-file
            allocate( names_fixl_observations_stat(nr_obs), &
                      names_fixl_observations_sou(nr_obs) )
            
            
            ! store the names of the observing stations to a new character variable with specific length
            ! dependent on the format of the .trp-file
            names_fixl_observing_stat= observing_stations_session % station
            
            ! store the names of the stations in the observations to a new character variable with specific length
            ! dependent on the format of the .trp-file
            names_fixl_observations_stat= rdlog_epoch % observations % station
            
            ! store the names of the sources in the observations to a new character variable with specific length
            ! dependent on the format of the .trp-file
            names_fixl_observations_sou= rdlog_epoch % observations % source
            
            
            ! allocate the variable for storing the strings of the observation time seconds
            allocate( obs_sec_str(nr_obs) )
            
            ! write the observation time seconds to the variable using the floating point representation with the
            ! desired field width and precision as desired by the .trp-file format, but with possible leading blanks
            ! instead of possible needed leading zeros
            ! note: A possible negative sign would lead to the output of asterisks as the field width would be to short then.
            !       Therefore abs() is used before writing the string to avoid this as seconds should anyway never be less than 0.
            !       Additionally a positive sign output is suppressed by using the format ss at writing the value to the string.
            
            ! loop over all observations
            do ind_obs= 1, nr_obs
                
                write( unit= obs_sec_str(ind_obs), fmt= '(ss, f4.1)' ) abs(rdlog_epoch % observations(ind_obs) % sec)
                
                ! in case the written value is below 10
                ! check if the first character is a blank and replace it by a zero if necessary
                ! note: Due to the floating point format f4.1 with which the string has been written the first character
                !       in the string has to be checked if it is a blank. In case the value that has been written is
                !       below 10 a blank appears at the position of the first character.
                if ( obs_sec_str(ind_obs)(1:1) == ' ' ) then
                    
                    ! replace the blank at the first position by a 0
                    obs_sec_str(ind_obs)(1:1)= '0'
                    
                    ! in case the written value is below 1
                    ! check if the second character is a blank and replace it by a zero if necessary
                    ! note: Depending on the processor of the system a 0 or no 0 is placed in front of the comma in case the value
                    !       that has been written is below 0. Therefore also a check of the second character is needed.
                    !       This check is only needed if already the first character has been a blank, otherwise the check is
                    !       not needed.
                    if ( obs_sec_str(ind_obs)(2:2) == ' ' ) then
                        ! replace the blank at the second position by a 0
                        obs_sec_str(ind_obs)(2:2)= '0'
                    end if
                    
                end if
                
            end do
            
            
            !-----------------------------------------------------------------------------------------------
            
            ! Define formats for writing the headers and the data to the .trp-file
            
            ! define format for writing the station header line
            ! notes on output format:
            ! 1x outputs one blank
            header_stat_fmt= '(a1, a11, a15, a14, a14, a10, a9, a7)'
            
            ! define format for writing the observation header line
            ! notes on output format:
            ! 1x outputs one blank
            header_obs_fmt= '(a1, a8, a12, a26, a10, a11, a9, a8, a6, a17, a16, a16, a15)'
            
            ! define format for writing the station data lines
            ! notes on output format:
            ! 1x outputs one blank
            ! ss supresses plus sign output for all following formats
            data_stat_fmt= '(ss, a1, 2x, a8, 2x, f13.4, 1x, f13.4, 1x, f13.4, 2x, f8.4, 1x, f8.4, 1x, f7.2)'
            
            ! define format for writing the observation data lines
            ! notes on output format:
            ! 1x outputs one blank
            ! ss supresses plus sign output for all following formats
            ! Note: In order to write the seconds with leading zeros, a workaround is necessary as leading zero output
            !       is only possible for integer and not for floating point format in Fortran.
            !       Formats used for this: The seconds are written to a string that is then modified to have leading zeros in case it is necessary.
            !                              A possible negative sign would lead to the output of asterisks as the field width would be to short then.
            !                              Therefore abs() is used before writing the string to avoid this as seconds should anyway never be less than 0.
            !                              Additionally a positive sign output is suppressed by using the format ss at writing the value to the string.
            !                              The format a4 prints the string of the seconds with the field width as specified in the trp-file format.
            ! Attention: The splitting of the writing of the seconds into an integer part and a fractional part is not possible since rounding effects for a value of
            !            e.g. 46.99 would lead to the output of 461.0 instead of 47.0!
            data_obs_fmt= '(ss, a1, 2x, i5, 4x, a8, 5x, i4.4, a1, i2.2, a1, i2.2, a1, i2.2, a1, i2.2, a1, a4, 2x, a8, 2x, f9.5, 1x, f8.5, 2x, f6.1, 1x, f5.1, 2x, es15.7, 1x, es15.7, 1x, es15.7, 1x, es15.7)'
            
            !-----------------------------------------------------------------------------------------------
            
            ! write to .trp-file
			
            ! define file header
            fileheader= 'TROPO_PATH_DELAY  Exchange format  v 1.2_TUVienna  Format version of 2014.07.10'
            
            ! write file header
			write(unit= file_unit_trp, fmt= '(a)') fileheader
			write(unit= file_unit_trp, fmt= '(a)') '#'
			
			! write file description
			write(unit= file_unit_trp, fmt= '(a)') '# This file contains tropospheric slant and zenith path delays and wet mapping factor.'
			write(unit= file_unit_trp, fmt= '(a)') '#'
			write(unit= file_unit_trp, fmt= '(a)') '# Created by: Armin Hofmeister, Technische Universität Wien'
			write(unit= file_unit_trp, fmt= '(a)') '#'
            
            ! write information about the program version
			write(unit= file_unit_trp, fmt= '(a)') '# Created with RADIATE program developed by Armin Hofmeister'
			write(unit= file_unit_trp, fmt= '(a, tr4, a, a)') '#', 'Version: ', RADIATE_version
			write(unit= file_unit_trp, fmt= '(a, tr4, a, a)') '#', 'Subversion: ', RADIATE_subversion
			write(unit= file_unit_trp, fmt= '(a)') '#'
            
            ! write time information about creation of the file
			write(unit= file_unit_trp, fmt= '(a)') '# File created on:'
            write(unit= file_unit_trp, fmt= '(a, tr4, a, 1x, a1, 1x, a)') '#',  date_formatted, '|', time_formatted
            write(unit= file_unit_trp, fmt= '(a)') '#'
            
            ! write information about the data content of the file
            write(unit= file_unit_trp, fmt= '(a)') '# Data content:'
            write(unit= file_unit_trp, fmt= '(a, tr4, a)') '#', 'Ray-tracing results for session:'
            write(unit= file_unit_trp, fmt= '(a, tr8, a)') '#', rdlog_epoch % session_name
            ! notes on output format:
            ! i0.1 outputs the minimal field width with at least one digit printed
            write(unit= file_unit_trp, fmt= '(a, tr4, a)') '#', 'Total number of observations:'
            write(unit= file_unit_trp, fmt= '(a, tr8, i0.1)') '#', nr_obs
            write(unit= file_unit_trp, fmt= '(a, tr4, a)') '#', 'Total number of different observing stations:'
            write(unit= file_unit_trp, fmt= '(a, tr8, i0.1)') '#', nr_stat
            write(unit= file_unit_trp, fmt= '(a)') '#'
            
            ! write information about the parameterization of the processing
            write(unit= file_unit_trp, fmt= '(a)') '# Chosen ray-tracing options:'
            write(unit= file_unit_trp, fmt= '(a, tr4, a)') '#', 'Assignment of observations to the epochs of the numerical weather model:'
            write(unit= file_unit_trp, fmt= '(a, tr8, a)') '#', epolog_mode
            write(unit= file_unit_trp, fmt= '(a, tr4, a)') '#', 'Method of vertical interpolation (horizontal treatment):'
            write(unit= file_unit_trp, fmt= '(a, tr8, a)') '#', interpolation_method
            write(unit= file_unit_trp, fmt= '(a, tr4, a)') '#', 'Method of numerical weather model usage:'
            write(unit= file_unit_trp, fmt= '(a, tr8, a)') '#', 'Full grid information used.'
            write(unit= file_unit_trp, fmt= '(a, tr4, a)') '#', 'Ray-tracing approach:'
            write(unit= file_unit_trp, fmt= '(a, tr8, a)') '#', raytr_method
            write(unit= file_unit_trp, fmt= '(a, tr4, a)') '#', 'Station coordinate catalogue-file:'
            write(unit= file_unit_trp, fmt= '(a, tr8, a)') '#', filename_stat_info
            write(unit= file_unit_trp, fmt= '(a)') '#'
            
            ! write information about the epochs and the horizontal resolution used during the processing
            write(unit= file_unit_trp, fmt= '(a)') '# Numerical weather model used:'
            
            ! write entries for all epochs
            ! note: size of "gridres" determines the number of epochs = gribfiles
            ! note: Implied loop for writing each line is used as it should be faster than a normal do loop.
            ! notes on output format:
            ! 1x outputs one blank
            ! i0.1 outputs the minimal field width with at least one digit printed
            ! ss supresses plus sign output for all following formats
            ! f5.3 --> form x.yyy
            write(unit= file_unit_trp, fmt= '(ss, a1, tr4, a, i0.1, a, a, 1x, a1, 1x, a, f5.3, a, f5.3, a)') ( '#', &
                                                                                        'Epoch ', ind_epoch, ': ', epoch_name(ind_epoch), '|', &
                                                                                        'Grib-file resolution: res_lat= ', gridres(ind_epoch) % res_lat, &
                                                                                        ' deg, res_lon= ', gridres(ind_epoch) % res_lon, ' deg', &
                                                                                        ind_epoch= 1, size(gridres) )
            
            write(unit= file_unit_trp, fmt= '(a)') '#'
            
            ! write information about the origin and use of the observational data
            write(unit= file_unit_trp, fmt= '(a)') '# Origin and use of the observational data:'
            write(unit= file_unit_trp, fmt= '(a, tr4, a)') '#', 'Observational data (scan number, source, observation time, site, azimuth, elevation)'
            write(unit= file_unit_trp, fmt= '(a, tr4, a)') '#', 'as well as pressure- and temperature-values provided in this file originate from the'
            write(unit= file_unit_trp, fmt= '(a, tr4, a)') '#', 'AzEl-file created by VieVS.'
            write(unit= file_unit_trp, fmt= '(a, tr4, a)') '#', 'The reported pressure- and temperature-values have not been used for the ray-tracing.'
            write(unit= file_unit_trp, fmt= '(a, tr4, a)') '#', 'Ray-tracing used temperature- and pressure-values only from the numerical weather'
            write(unit= file_unit_trp, fmt= '(a, tr4, a)') '#', 'model and the standard atmosphere model.'
            write(unit= file_unit_trp, fmt= '(a)') '#'
            write(unit= file_unit_trp, fmt= '(a)') '#'
            
            ! write detailed file content and format rules of the trp-file format
            write(unit= file_unit_trp, fmt= '(a)') '# ---------------------------------------------------------------------------------------'
            write(unit= file_unit_trp, fmt= '(a)') '# ----------------------------------- File description ----------------------------------'
            write(unit= file_unit_trp, fmt= '(a)') '# ---------------------------------------------------------------------------------------'
            write(unit= file_unit_trp, fmt= '(a)') '#'
            write(unit= file_unit_trp, fmt= '(a)') '# 1) Header record -- the first record line should have a signature:'
            write(unit= file_unit_trp, fmt= '(a)') '#    TROPO_PATH_DELAY  Exchange format  v N.N_additional  Format version of YYYY.MM.DD'
            write(unit= file_unit_trp, fmt= '(a)') '#    stating the version and date of the format'
            write(unit= file_unit_trp, fmt= '(a)') '#'
            write(unit= file_unit_trp, fmt= '(a)') '#    The header record allows to distinguish a valid file in the TROPO_PATH_DELAY'
            write(unit= file_unit_trp, fmt= '(a)') '#    format from files in other formats.'
            write(unit= file_unit_trp, fmt= '(a)') '#'
            write(unit= file_unit_trp, fmt= '(a)') '# 2) E-records section'
            write(unit= file_unit_trp, fmt= '(a)') '#    E-record  has letter E in the first field.'
            write(unit= file_unit_trp, fmt= '(a)') '#              It specifies the experiment name.'
            write(unit= file_unit_trp, fmt= '(a)') '#'
            write(unit= file_unit_trp, fmt= '(a)') '# 3) H-records section'
            write(unit= file_unit_trp, fmt= '(a)') '#    H-record  has letter H in the first field.'
            write(unit= file_unit_trp, fmt= '(a)') '#              It specifies the secondary experiment name.'
            write(unit= file_unit_trp, fmt= '(a)') '#'
            write(unit= file_unit_trp, fmt= '(a)') '# 4) M-records section'
            write(unit= file_unit_trp, fmt= '(a)') '#    M-record  has letter M in the first field.'
            write(unit= file_unit_trp, fmt= '(a)') '#              The M record keeps the model identifier.'
            write(unit= file_unit_trp, fmt= '(a)') '#'
            write(unit= file_unit_trp, fmt= '(a)') '# 5) U-records section'
            write(unit= file_unit_trp, fmt= '(a)') '#    U-record  has letter U in the first field.'
            write(unit= file_unit_trp, fmt= '(a)') '#              The U record instructs analysis software how to use the data.'
            write(unit= file_unit_trp, fmt= '(a)') '#              Analysis software may override this recommendation.'
            write(unit= file_unit_trp, fmt= '(a)') '#'
            write(unit= file_unit_trp, fmt= '(a)') '# 6) S-records section'
            write(unit= file_unit_trp, fmt= '(a)') '#    S-record  has letter S in the first field.'
            write(unit= file_unit_trp, fmt= '(a)') '#              It defines the name and coordinates of each site. A site cannot be'
            write(unit= file_unit_trp, fmt= '(a)') '#              defined more than once.'
            write(unit= file_unit_trp, fmt= '(a)') '#'
            write(unit= file_unit_trp, fmt= '(a)') '# 7) O-records section'
            write(unit= file_unit_trp, fmt= '(a)') '#    O-record  has letter O in the first field.'
            write(unit= file_unit_trp, fmt= '(a)') '#              It defines the circumstances of the observation:'
            write(unit= file_unit_trp, fmt= '(a)') '#              the time tag, site ID, azimuth and elevation in the local'
            write(unit= file_unit_trp, fmt= '(a)') '#              topocentric coordinate system, surface atmospheric pressure and'
            write(unit= file_unit_trp, fmt= '(a)') '#              surface air temperature, slant and zenith path delay and'
            write(unit= file_unit_trp, fmt= '(a)') '#              wet mapping factor.'
            write(unit= file_unit_trp, fmt= '(a)') '#'
            write(unit= file_unit_trp, fmt= '(a)') '# Records, which start from # character, are considered as comments.'
            write(unit= file_unit_trp, fmt= '(a)') '#'
            write(unit= file_unit_trp, fmt= '(a)') '# Format of an S-record:'
            write(unit= file_unit_trp, fmt= '(a)') '# ----------------------'
            write(unit= file_unit_trp, fmt= '(a)') '#'
            write(unit= file_unit_trp, fmt= '(a)') '# Field  1:1   A1    -- Record ID. Should be letter S.'
            write(unit= file_unit_trp, fmt= '(a)') '# field  2:3   a2       delimiter: blanks.'
            write(unit= file_unit_trp, fmt= '(a)') '# Field  4:11  A8    -- 8-letter long site identifier.'
            write(unit= file_unit_trp, fmt= '(a)') '#                       May contain any characters with decimal codes 32-255, but'
            write(unit= file_unit_trp, fmt= '(a)') '#                       blanks are allowed only at the end of the site identifier.'
            write(unit= file_unit_trp, fmt= '(a)') '#                       This site identifier should be unique among S-records.'
            write(unit= file_unit_trp, fmt= '(a)') '#                       This field should not necessarily have a special meaning.'
            write(unit= file_unit_trp, fmt= '(a)') '# field 12:13  a2       delimiter: blanks.'
            write(unit= file_unit_trp, fmt= '(a)') '# Field 14:26  F13.4 -- X site coordinate in a crust fixed reference frame.'
            write(unit= file_unit_trp, fmt= '(a)') '#                       Unit: meter'
            write(unit= file_unit_trp, fmt= '(a)') '# field 27:27  a1       delimiter: blank.'
            write(unit= file_unit_trp, fmt= '(a)') '# Field 28:40  F13.4 -- Y site coordinate in a crust fixed reference frame.'
            write(unit= file_unit_trp, fmt= '(a)') '#                       Unit: meter'
            write(unit= file_unit_trp, fmt= '(a)') '# field 41:41  a1       delimiter: blank.'
            write(unit= file_unit_trp, fmt= '(a)') '# Field 42:54  F13.4 -- Z site coordinate in a crust fixed reference frame.'
            write(unit= file_unit_trp, fmt= '(a)') '#                       Unit: meter'
            write(unit= file_unit_trp, fmt= '(a)') '# field 55:56  a2       delimiter: blanks.'
            write(unit= file_unit_trp, fmt= '(a)') '# Field 57:64  F8.4  -- Geodetic latitude of site,'
            write(unit= file_unit_trp, fmt= '(a)') '#                       interval: [-90; 90], positive to north.'
            write(unit= file_unit_trp, fmt= '(a)') '#                       Unit: degree'
            write(unit= file_unit_trp, fmt= '(a)') '# field 65:65  a1       delimiter: blank.'
            write(unit= file_unit_trp, fmt= '(a)') '# Field 66:73  F8.4  -- Geodetic longitude of site,'
            write(unit= file_unit_trp, fmt= '(a)') '#                       interval: [0; 360[, increasing towards east.'
            write(unit= file_unit_trp, fmt= '(a)') '#                       Unit: degree'
            write(unit= file_unit_trp, fmt= '(a)') '# field 74:74  a1       delimiter: blank.'
            write(unit= file_unit_trp, fmt= '(a)') '# Field 75:81  F7.2  -- Height of site above the reference ellipsoid.'
            write(unit= file_unit_trp, fmt= '(a)') '#                       Unit: meter'
            write(unit= file_unit_trp, fmt= '(a)') '#'
            write(unit= file_unit_trp, fmt= '(a)') '# Format of an O-record:'
            write(unit= file_unit_trp, fmt= '(a)') '# ----------------------'
            write(unit= file_unit_trp, fmt= '(a)') '#'
            write(unit= file_unit_trp, fmt= '(a)') '# Field  1:1    A1      -- Record ID. Should be letter O.'
            write(unit= file_unit_trp, fmt= '(a)') '# field  2:3    a2         delimiter: blanks.'
            write(unit= file_unit_trp, fmt= '(a)') '# Field  4:8    I5      -- Scan number'
            write(unit= file_unit_trp, fmt= '(a)') '# field  9:12   a4         delimiter: blanks.'
            write(unit= file_unit_trp, fmt= '(a)') '# Field 13:20   A8      -- Source name'
            write(unit= file_unit_trp, fmt= '(a)') '# field 21:25   a5         delimiter: blanks.'
            write(unit= file_unit_trp, fmt= '(a)') '# Field 26:46   A21     -- Calendar date in TAI in YYYY.MM.DD-hh:mm:ss.s format'
            write(unit= file_unit_trp, fmt= '(a)') '#                          Calendar date format:'
            write(unit= file_unit_trp, fmt= '(a)') '#                          field 26:29 -- I4   Year'
            write(unit= file_unit_trp, fmt= '(a)') '#                          field 30:30 -- a1   delimiter: .'
            write(unit= file_unit_trp, fmt= '(a)') '#                          field 31:32 -- I2   Month'
            write(unit= file_unit_trp, fmt= '(a)') '#                          field 33:33 -- a1   delimiter: .'
            write(unit= file_unit_trp, fmt= '(a)') '#                          field 34:35 -- I2   Day'
            write(unit= file_unit_trp, fmt= '(a)') '#                          field 36:36 -- a1   delimiter: -'
            write(unit= file_unit_trp, fmt= '(a)') '#                          field 37:38 -- I2   Hour'
            write(unit= file_unit_trp, fmt= '(a)') '#                          field 39:39 -- a1   delimiter: :'
            write(unit= file_unit_trp, fmt= '(a)') '#                          field 40:41 -- I2   Minute'
            write(unit= file_unit_trp, fmt= '(a)') '#                          field 42:42 -- a1   delimiter: :'
            write(unit= file_unit_trp, fmt= '(a)') '#                          field 43:46 -- F4.1 Seconds'
            write(unit= file_unit_trp, fmt= '(a)') '# field 47:48   a2         delimiter: blanks.'
            write(unit= file_unit_trp, fmt= '(a)') '# Field 49:56   A8      -- Site identifier'
            write(unit= file_unit_trp, fmt= '(a)') '# field 57:58   a2         delimiter: blanks.'
            write(unit= file_unit_trp, fmt= '(a)') '# Field 59:67   F9.5    -- Azimuth'
            write(unit= file_unit_trp, fmt= '(a)') '#                          Unit: degree'
            write(unit= file_unit_trp, fmt= '(a)') '# field 68:68   a1         delimiter: blank.'
            write(unit= file_unit_trp, fmt= '(a)') '# Field 69:76   F8.5    -- Elevation angle (outgoing)'
            write(unit= file_unit_trp, fmt= '(a)') '#                          Unit: degree'
            write(unit= file_unit_trp, fmt= '(a)') '# field 77:78   a2         delimiter: blanks.'
            write(unit= file_unit_trp, fmt= '(a)') '# Field 79:84   F6.1    -- Atmospheric surface pressure'
            write(unit= file_unit_trp, fmt= '(a)') '#                          Unit: hPa (mbar)'
            write(unit= file_unit_trp, fmt= '(a)') '# field 85:85   a1         delimiter: blank.'
            write(unit= file_unit_trp, fmt= '(a)') '# Field 86:90   F5.1    -- Surface air temperature'
            write(unit= file_unit_trp, fmt= '(a)') '#                          Unit: degree Celsius'
            write(unit= file_unit_trp, fmt= '(a)') '# field 91:92   a2         delimiter: blanks.'
            write(unit= file_unit_trp, fmt= '(a)') '# Field 93:107  1PD15.7 -- Slant total delay'
            write(unit= file_unit_trp, fmt= '(a)') '#                          Unit: seconds of time'
            write(unit= file_unit_trp, fmt= '(a)') '# field 108:108 a1         delimiter: blank.'
            write(unit= file_unit_trp, fmt= '(a)') '# Field 109:123 1PD15.7 -- Wet mapping factor'
            write(unit= file_unit_trp, fmt= '(a)') '#                          Tilting of the wet part of the atmosphere in the form'
            write(unit= file_unit_trp, fmt= '(a)') '#                          of a wet mapping function:'
            write(unit= file_unit_trp, fmt= '(a)') '#                          [wet_delay_along_the_path / wet_zenith_delay]'
            write(unit= file_unit_trp, fmt= '(a)') '# field 124:124 a1         delimiter: blank.'
            write(unit= file_unit_trp, fmt= '(a)') '# Field 125:139 1PD15.7 -- Hydrostatic zenith delay'
            write(unit= file_unit_trp, fmt= '(a)') '#                          Unit: seconds of time'
            write(unit= file_unit_trp, fmt= '(a)') '# field 140:140 a1         delimiter: blank.'
            write(unit= file_unit_trp, fmt= '(a)') '# Field 141:155 1PD15.7 -- Wet zenith delay'
            write(unit= file_unit_trp, fmt= '(a)') '#                          Unit: seconds of time'
            write(unit= file_unit_trp, fmt= '(a)') '#'
            write(unit= file_unit_trp, fmt= '(a)') '#######################################################################################################'
            write(unit= file_unit_trp, fmt= '(a)') '#'
            
            ! write E-records section
            ! notes on output format:
            ! 2x outputs two blanks
            write(unit= file_unit_trp, fmt= '(a1, 2x, a1, a)') 'E', '$', rdlog_epoch % session_name
            
            ! write H-records section
            ! notes on output format:
            ! 2x outputs two blanks
            write(unit= file_unit_trp, fmt= '(a1, 2x, a1, a)') 'H', '$', rdlog_epoch % session_name
            write(unit= file_unit_trp, fmt= '(a)') '#'
            
            ! write M-records section
            write(unit= file_unit_trp, fmt= '(6a, 1x, 2a)') 'M  Ray-tracing results from RADIATE program (version: ', &
                                                    RADIATE_version, &
                                                    ', subversion: ', &
                                                    RADIATE_subversion, &
                                                    ', created: ', &
                                                    date_formatted, &
                                                    time_formatted, &
                                                    ') developed by Armin Hofmeister, Technische Universität Wien.'
            write(unit= file_unit_trp, fmt= '(a)') '#'
            
            ! write U-records section
            write(unit= file_unit_trp, fmt= '(a)') 'U  NONE'
            write(unit= file_unit_trp, fmt= '(a)') '#'
            
            ! write S-records section
            ! write header for S-records
            write(unit= file_unit_trp, fmt= header_stat_fmt) '#', &
                                                             'Site|', &
                                                             'X-coordinate|', &
                                                             'Y-coordinate|', &
                                                             'Z-coordinate|', &
                                                             'Latitude|', &
                                                             'Longit.|', &
                                                             'Height'
            
            ! write S-records with specified format
            
            ! Attention: The station names and sources names are always written to the trp-file using a width of 8 characters. In case of station and sources names
            !            that are longer, the names will be truncated.
            !            The length of 8 characters is defined by the trp-format.
            
            ! write S-records for all observing stations in the session
            ! note: Implied loop for writing each line is used as it should be faster than a normal do loop.
            write(unit= file_unit_trp, fmt= data_stat_fmt) ( 'S', &
                                                             observing_stations_session(ind_stat) % station, &
                                                             x_stat(ind_stat), &
                                                             y_stat(ind_stat), &
                                                             z_stat(ind_stat), &
                                                             observing_stations_session(ind_stat) % lat_ell, &
                                                             observing_stations_session(ind_stat) % lon_ell, &
                                                             observing_stations_session(ind_stat) % h_ell, ind_stat= 1, nr_stat )
            
            write(unit= file_unit_trp, fmt= '(a)') '#'
            write(unit= file_unit_trp, fmt= '(a)') '#'
            
            
            ! write O-records section
            ! write header for O-records
            write(unit= file_unit_trp, fmt= header_obs_fmt) '#', &
                                                            'Scan|', &
                                                            'Source|', &
                                                            'Time tag in TAI|', &
                                                            'Site|', &
                                                            'Azimuth|', &
                                                            'Elev.|', &
                                                            'Pres.|', &
                                                            'Temp.|', &
                                                            'Slant total del|', &
                                                            'Wet map. factor|', &
                                                            'Hydr zenith del|', &
                                                            'Wet zenith del.'
            
            ! write O-records with specified format
            ! note: Implied loop for writing each line is used as it should be faster than a normal do loop.
            write(unit= file_unit_trp, fmt= data_obs_fmt) ( 'O', &
                                                            rdlog_epoch % observations(ind_obs) % scannr, &
                                                            rdlog_epoch % observations(ind_obs) % source, &
                                                            rdlog_epoch % observations(ind_obs) % year, &
                                                            '.', &
                                                            date_obs(ind_obs) % month, &
                                                            '.', &
                                                            date_obs(ind_obs) % day, &
                                                            '-', &
                                                            rdlog_epoch % observations(ind_obs) % hour, &
                                                            ':', &
                                                            rdlog_epoch % observations(ind_obs) % min, &
                                                            ':', &
                                                            obs_sec_str(ind_obs), & ! print the preprocessed string of the observation seconds that has no sign, but possible leading zeros
                                                            rdlog_epoch % observations(ind_obs) % station, &
                                                            azimuth_deg(ind_obs), &
                                                            elevation_deg(ind_obs), &
                                                            rdlog_epoch % observations(ind_obs) % pres, &
                                                            rdlog_epoch % observations(ind_obs) % temp, &
                                                            slant_total_delay_s(ind_obs), &
                                                            rdlog_epoch % delay(ind_obs) % mf_w, &
                                                            zenith_h_delay_s(ind_obs), &
                                                            zenith_w_delay_s(ind_obs), ind_obs= 1, nr_obs )
                
            ! rewrite file header record at last line of file
            write(unit= file_unit_trp, fmt= '(a)') fileheader
        
        end if
        
        ! close the created .trp-file
        ! note: close() for an unit that has not been opened, e.g. due to an error, has no effect, e.g. as file could not be created
        close(unit= file_unit_trp)

    end subroutine create_trp_global
    
end module module_create_trp_global