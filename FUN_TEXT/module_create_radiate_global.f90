! module_create_radiate_global.f90
    

!****************************************************************************
!
! PURPOSE:  Module containing subroutine for the creation of the .radiate-file.
!
!           This subroutine carries creates a .radiate-file containing all important results of the ray-tracing for the session.
!           The original AZEL-file information describing each observation is enhanced by the ray-tracing results.
!
!           The sorting order of the observations in the output is preserved in the way as it has been defined in subroutine "combine_and_sort_rd":
!               1.) chronological mjd
!               2.) if same mjd then alphabetical order of station name
!               3.) if same mjd and same station name then alphabetical order of source names
!
!           Note: Created as module to avoid explicit inferface for the subroutine as needed because of assumed-shape dummy argument in subroutine.
!
!
! INPUT:
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
!           % grib_epoch_mjd... note: deallocated at this stage:
!                                     prior used as vector containing the mjd of the grib-file epoch which has been used to calculate
!                                     the delay of each observation
!
!        save_path.............. specifies path to the saving directory for the .radiate-file
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
! 05.05.2015: create the Fortran-file based on the Matlab-file "write_rtr_calc_sessionwise_global.m"
! 06.05.2015: programming
! 11.05.2015: add comments
! 12.05.2015: programming
! 13.05.2015: programming
! 02.06.2015: correct comments
! 03.06.2015: change to implied loop for writing the data lines
! 09.06.2015: programming
! 11.06.2015: programming
! 17.12.2015: changes due to azel file extension by water vapour pressure
!             changes due to the output of the total mapping factor and the meteorological data at the station position from the NWM to the .radiate-file
! 12.01.2016: add comments
! 13.01.2016: correct comments
!
!****************************************************************************

module module_create_radiate_global

contains
    
    subroutine create_radiate_global( rdlog_epoch, &
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
        use module_type_definitions, only: rdlog_epoch_type, gridres_type
        use module_constants, only: len_statname, len_souname, C2K
        
    
        !----------------------------------------------------------------------------
    
        ! Variable definitions
        implicit none
    
    
        ! dummy arguments
        !----------------
    
        ! INPUT
    
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
        
        ! define variable for storing the field width for writing the header fields for station names and the source names
        ! note: In order to avoid a character string that is to short for string representation of the number,
        !       the string length is set to the size of the station or source name string +1 for the "|" in the header.
        !       Later trim() will be used to get rid of the trailing blanks, so that the field width is correctly written as character without blanks.
        !       In this way complicated calculation of the needed string length is avoided by the use of an acceptable setting of a too large length.
        character(len=len_statname+1) :: len_header_statname_str
        character(len=len_souname+1) :: len_header_souname_str
        
        ! define variable for storing the format for writing the data description header-lines
        character(len=:), allocatable :: descr_fmt
        
        ! define variable for storing the format for writing the data of the observations and their ray-tracing results
        character(len=:), allocatable :: data_fmt
        
        ! define variable for storing the converted temperature in [°C] at the station position interpolated from the numerical weather model
        ! note: The variable will be allocated later according to the number of observations.
        double precision, dimension(:), allocatable :: T_stat_out
        
        ! define variable for storing the current date
        character(len=8) :: date_value
        character(len=10) :: date_formatted
    
        ! define variable for storing the current time
        character(len=10) :: time_value
        character(len=12) :: time_formatted
        
        ! variable for storing the unit of the file
        integer, parameter :: file_unit_radiate= 1
    
        ! variable for storing the opening status of the file
        integer :: open_status
        
        ! define loop variable
        integer :: ind_epoch
        
        ! define loop variable
        integer :: ind_obs
        
        
        ! CONSTANTS
        
        ! see "module_constants" for "len_statname", "len_souname" and "C2K"
        
        
        !----------------------------------------------------------------------------
        
        !============================================================================
        ! Body of the subroutine
        !============================================================================
        
        ! Create .radiate-file
        ! The .radiate-file contains information and results of the ray-tracing processing for all observations of the session.
        
        ! get the current date and time
        call DATE_AND_TIME(date= date_value, time= time_value)
        
        ! format the date and time
        date_formatted= date_value(1:4) // '-' // date_value(5:6) // '-' // date_value(7:8)
        time_formatted= time_value(1:2) // ':' // time_value(3:4) // ':' // time_value(5:)
        
        ! create new .radiate-ascii-file for storing the data
        ! note: an old file is replaced by the new one as status is set to 'replace'
        open(unit= file_unit_radiate, file= save_path // rdlog_epoch % session_name // '.radiate', action= 'write', status= 'replace', iostat= open_status)
        
        ! Check if textfile can't be opened (iostat /= 0 --> error)
        if (open_status /= 0) then
            
            ! report warning message
            write(unit= *, fmt= '(tr4, a, a, a, /)') 'Warning: Problem with creating .radiate-file: "', save_path // rdlog_epoch % session_name // '.radiate"', '! No file created!'    
        
        else
            
            ! get the total number of observations in the specific session = size of the combined (and time interpolated) epologs = rdlog_epoch
            ! note: use size() of substructure "% observations"
            nr_obs= size(rdlog_epoch % observations)
        
            ! define field width for writing the header fields for station names and the source names
            ! note: write to character-variable to use it in format specification
            !       + 1 as "|" is printed at the end
            ! notes on output format:
            ! i0.1 outputs the minimal field width with at least one digit printed
            ! ss supresses plus sign output for all following formats
            write(unit= len_header_statname_str, fmt= '(ss, i0.1)') (len_statname + 1)
            write(unit= len_header_souname_str, fmt= '(ss, i0.1)') (len_souname + 1)
        
            ! define format of output for writing the data description header-lines
            ! notes on output format:
            ! a[w] outputs blanks to the left in case the length w is specified longer than needed for the string
            ! note: use variables to specify field width for header of station names and sources
            descr_fmt= '(a1, a6, a12, a5, a4, a3, a3, a6, a' // trim(len_header_statname_str) // ', a21, a21, a' // trim(len_header_souname_str) // ', a7, a8, a7, 6(a10), 2(a13), a10, 3(a11), a7, a8, a6)'
        
            ! define format for writing the data of the observations and their ray-tracing results
            ! notes on output format:
            ! 1x outputs one blank
            ! ss supresses plus sign output for all following formats
            ! note: for the station and source name use "a" without a field width
            !       to output the whole string length, which should be sized according to
            !       the constants "len_statname" and "len_souname". The header has been sized according to theses values.
            ! note: output precisions: for delays and bending: 0.1 mm, for elevation angles: 1e-7 rad (to reach precision of 0.1 mm in delays)
            !                          mapping factors: 1e-5 (to reach precision of 0.1 mm in delays)
            ! attention: Due the use of implied loop for writing the observation and ray-tracing data for each line,
            !            it is not possible to use the repeated format-syntax, e.g. 6(f9.4, 1x), within the format specification
            !            as this leads to an error in writing the lines.
            data_fmt='(1x, ss, i5, 1x, f11.5, 1x, i4, 1x, i3, 1x, i2, 1x, i2, 1x, f5.2, 1x, a, 1x, f20.15, 1x, f20.15, 1x, a, 1x, f6.2, 1x, f7.2, 1x, f6.2, 1x, f9.4, 1x, f9.4, 1x, f9.4, 1x, f9.4, 1x, f9.4, 1x, f9.4, 1x, f12.7, 1x, f12.7, 1x, f9.4, 1x, f10.5, 1x, f10.5, 1x, f10.5, 1x, f6.2, 1x, f7.2, 1x, f6.2)'
            
            
            ! convert the temperature in the "meteo_stat_out"-structure interpolated at the station position from the numerical weather model from [K] to [°C]
            ! in order to be consistent with the unit of the temperature at the station from the AzEl-file
            ! note: see "module_constants" for the temperature scale offset "C2K"
            !       [°C] = [K] - 273.15 = [°C] = [K] - C2K
            ! note: A different variable than the input variable is needed since an input variable can not be redefined.
            
            ! allocate the variable
            ! note: The size is according to the number of observations.
            allocate( T_stat_out(nr_obs) )
            
            ! assign the temperature in [°C]
            T_stat_out = rdlog_epoch % meteo_stat_out(:) % T - C2K
            
            
            !-----------------------------------------------------------------------------------------------
            
            ! write to .radiate-file
            
            ! write information about the data content of the file
            write(unit= file_unit_radiate, fmt= '(a)') '%*****************************************************************************************************************'
            
            ! write header
			write(unit= file_unit_radiate, fmt= '(a)') '% RADIATE format v 2.0'
			write(unit= file_unit_radiate, fmt= '(a)') '%'
			
			! write file description
			write(unit= file_unit_radiate, fmt= '(a)') '% This file contains ray-tracing results.'
			write(unit= file_unit_radiate, fmt= '(a)') '%'
			write(unit= file_unit_radiate, fmt= '(a)') '% Created by: Armin Hofmeister, Technische Universität Wien'
			write(unit= file_unit_radiate, fmt= '(a)') '%'
            
            ! write information about the program version
            write(unit= file_unit_radiate, fmt= '(a)') '% Created with RADIATE program developed by Armin Hofmeister'
            write(unit= file_unit_radiate, fmt= '(a, tr4, a, a)') '%', 'Version: ', RADIATE_version
            write(unit= file_unit_radiate, fmt= '(a, tr4, a, a)') '%', 'Subversion: ', RADIATE_subversion
            write(unit= file_unit_radiate, fmt= '(a)') '%'
			
            ! write time information about creation of the file
            write(unit= file_unit_radiate, fmt= '(a)') '% File created on:'
            write(unit= file_unit_radiate, fmt= '(a, tr4, a, 1x, a1, 1x, a)') '%',  date_formatted, '|', time_formatted
            write(unit= file_unit_radiate, fmt= '(a)') '%'
            
            ! write information about the data content of the file
            write(unit= file_unit_radiate, fmt= '(a)') '% Data content:'
            write(unit= file_unit_radiate, fmt= '(a, tr4, a)') '%', 'Ray-tracing results for session:'
            write(unit= file_unit_radiate, fmt= '(a, tr8, a)') '%', rdlog_epoch % session_name
            ! notes on output format:
            ! i0.1 outputs the minimal field width with at least one digit printed
            write(unit= file_unit_radiate, fmt= '(a, tr4, a)') '%', 'Total number of observations:'
            write(unit= file_unit_radiate, fmt= '(a, tr8, i0.1)') '%', nr_obs
            write(unit= file_unit_radiate, fmt= '(a)') '%'
            
            ! write information about the parameterization of the processing
            write(unit= file_unit_radiate, fmt= '(a)') '% Chosen ray-tracing options:'
            write(unit= file_unit_radiate, fmt= '(a, tr4, a)') '%', 'Assignment of observations to the epochs of the numerical weather model:'
            write(unit= file_unit_radiate, fmt= '(a, tr8, a)') '%', epolog_mode
            write(unit= file_unit_radiate, fmt= '(a, tr4, a)') '%', 'Method of vertical interpolation (horizontal treatment):'
            write(unit= file_unit_radiate, fmt= '(a, tr8, a)') '%', interpolation_method
            write(unit= file_unit_radiate, fmt= '(a, tr4, a)') '%', 'Method of numerical weather model usage:'
            write(unit= file_unit_radiate, fmt= '(a, tr8, a)') '%', 'Full grid information used'
            write(unit= file_unit_radiate, fmt= '(a, tr4, a)') '%', 'Ray-tracing approach:'
            write(unit= file_unit_radiate, fmt= '(a, tr8, a)') '%', raytr_method
            write(unit= file_unit_radiate, fmt= '(a, tr4, a)') '%', 'Station coordinate catalogue-file:'
            write(unit= file_unit_radiate, fmt= '(a, tr8, a)') '%', filename_stat_info
            write(unit= file_unit_radiate, fmt= '(a)') '%'
            
            ! write information about the epochs and the horizontal resolution used during the processing
            write(unit= file_unit_radiate, fmt= '(a)') '% Numerical weather model used:'
            
            ! write entries for all epochs
            ! note: size of "gridres" determines the number of epochs = gribfiles
            ! note: Implied loop for writing each line is used as it should be faster than a normal do loop.
            ! notes on output format:
            ! 1x outputs one blank
            ! i0.1 outputs the minimal field width with at least one digit printed
            ! ss supresses plus sign output for all following formats
            ! f5.3 --> form x.yyy
            write(unit= file_unit_radiate, fmt= '(ss, a1, tr4, a, i0.1, a, a, 1x, a1, 1x, a, f5.3, a, f5.3, a)') ( '%', &
                                                                                        'Epoch ', ind_epoch, ': ', epoch_name(ind_epoch), '|', &
                                                                                        'Grib-file resolution: res_lat= ', gridres(ind_epoch) % res_lat, &
                                                                                        ' deg, res_lon= ', gridres(ind_epoch) % res_lon, ' deg', &
                                                                                        ind_epoch= 1, size(gridres) )
            
            write(unit= file_unit_radiate, fmt= '(a)') '%'
            
            ! write information about the origin and use of the observational data
            write(unit= file_unit_radiate, fmt= '(a)') '% Origin and use of the observational data:'
            write(unit= file_unit_radiate, fmt= '(a, tr4, a)') '%', 'Observational data (scan number, observation time, station, azimuth, elevation, source)'
            write(unit= file_unit_radiate, fmt= '(a, tr4, a)') '%', 'as well as temperature-, pressure- and water vapour pressure-values provided in this file'
            write(unit= file_unit_radiate, fmt= '(a, tr4, a)') '%', 'originate from the AzEl-file created by VieVS.'
            write(unit= file_unit_radiate, fmt= '(a, tr4, a)') '%', 'The reported temperature-, pressure- and water vapour pressure-values from the numerical'
            write(unit= file_unit_radiate, fmt= '(a, tr4, a)') '%', 'weather model at the station positions have been determined from the numerical weather'
            write(unit= file_unit_radiate, fmt= '(a, tr4, a)') '%', 'model by interpolation within the ray-tracing application.'
            write(unit= file_unit_radiate, fmt= '(a, tr4, a)') '%', 'Ray-tracing used temperature-, pressure- and water vapour-values only from the numerical'
            write(unit= file_unit_radiate, fmt= '(a, tr4, a)') '%', 'weather model.'
            write(unit= file_unit_radiate, fmt= '(a)') '%'
            write(unit= file_unit_radiate, fmt= '(a)') '%'
            
            ! write information about the data in the different columns
            write(unit= file_unit_radiate, fmt= '(a)') '% Data columns:'
            write(unit= file_unit_radiate, fmt= '(a, tr4, a)') '%', '1  .... scannumber'
            write(unit= file_unit_radiate, fmt= '(a, tr4, a)') '%', '2  .... mjd'
            write(unit= file_unit_radiate, fmt= '(a, tr4, a)') '%', '3  .... year'
            write(unit= file_unit_radiate, fmt= '(a, tr4, a)') '%', '4  .... day of year'
            write(unit= file_unit_radiate, fmt= '(a, tr4, a)') '%', '5  .... hour'
            write(unit= file_unit_radiate, fmt= '(a, tr4, a)') '%', '6  .... min'
            write(unit= file_unit_radiate, fmt= '(a, tr4, a)') '%', '7  .... sec'
            write(unit= file_unit_radiate, fmt= '(a, tr4, a)') '%', '8  .... station'
            write(unit= file_unit_radiate, fmt= '(a, tr4, a)') '%', '9  .... azimuth in [rad]'
            write(unit= file_unit_radiate, fmt= '(a, tr4, a)') '%', '10 .... outgoing elevation angle (calculated, theoretical) in [rad]'
            write(unit= file_unit_radiate, fmt= '(a, tr4, a)') '%', '11 .... source'
            write(unit= file_unit_radiate, fmt= '(a, tr4, a)') '%', '12 .... temperature at station in [°C]'
            write(unit= file_unit_radiate, fmt= '(a, tr4, a)') '%', '13 .... pressure at station in [hPa]'
            write(unit= file_unit_radiate, fmt= '(a, tr4, a)') '%', '14 .... water vapour pressure at station in [hPa]'
            write(unit= file_unit_radiate, fmt= '(a, tr4, a)') '%', '15 .... zenith total delay in [m]'
            write(unit= file_unit_radiate, fmt= '(a, tr4, a)') '%', '16 .... zenith hydrostatic delay in [m]'
            write(unit= file_unit_radiate, fmt= '(a, tr4, a)') '%', '17 .... zenith wet delay in [m]'
            write(unit= file_unit_radiate, fmt= '(a, tr4, a)') '%', '18 .... slant total delay including geometric bending effect in [m]'
            write(unit= file_unit_radiate, fmt= '(a, tr4, a)') '%', '19 .... slant hydrostatic delay including geometric bending effect in [m]'
            write(unit= file_unit_radiate, fmt= '(a, tr4, a)') '%', '20 .... slant wet delay in [m]'
            write(unit= file_unit_radiate, fmt= '(a, tr4, a)') '%', '21 .... elevation angle at station in [rad]'
            write(unit= file_unit_radiate, fmt= '(a, tr4, a)') '%', '22 .... outgoing elevation angle from ray-tracing in [rad]'
            write(unit= file_unit_radiate, fmt= '(a, tr4, a)') '%', '23 .... geometric bending effect in [m]'
            write(unit= file_unit_radiate, fmt= '(a, tr4, a)') '%', '24 .... total mapping factor (includes treatment of geometric bending effect) [total_delay_along_the_path / total_zenith_delay]'
            write(unit= file_unit_radiate, fmt= '(a, tr4, a)') '%', '25 .... hydrostatic mapping factor (includes treatment of geometric bending effect) [hydrostatic_delay_along_the_path / hydrostatic_zenith_delay]'
            write(unit= file_unit_radiate, fmt= '(a, tr4, a)') '%', '26 .... wet mapping factor [wet_delay_along_the_path / wet_zenith_delay]'
            write(unit= file_unit_radiate, fmt= '(a, tr4, a)') '%', '27 .... temperature at the station position interpolated from the numerical weather model in [°C]'
            write(unit= file_unit_radiate, fmt= '(a, tr4, a)') '%', '28 .... pressure at the station position interpolated from the numerical weather model in [hPa]'
            write(unit= file_unit_radiate, fmt= '(a, tr4, a)') '%', '29 .... water vapour pressure at the station position interpolated from the numerical weather model in [hPa]'
            write(unit= file_unit_radiate, fmt= '(a)') '%*****************************************************************************************************************'
            write(unit= file_unit_radiate, fmt= '(a)') '%'
            write(unit= file_unit_radiate, fmt= '(a)') '%'
            write(unit= file_unit_radiate, fmt= descr_fmt) '%', &
                                                           '1|', &
                                                           '2|', &
                                                           '3|', &
                                                           '4|', &
                                                           '5|', &
                                                           '6|', &
                                                           '7|', &
                                                           '8|', &
                                                           '9|', &
                                                           '10|', &
                                                           '11|', &
                                                           '12|', &
                                                           '13|', &
                                                           '14|', &
                                                           '15|', &
                                                           '16|', &
                                                           '17|', &
                                                           '18|', &
                                                           '19|', &
                                                           '20|', &
                                                           '21|', &
                                                           '22|', &
                                                           '23|', &
                                                           '24|', &
                                                           '25|', &
                                                           '26|', &
                                                           '27|', &
                                                           '28|', &
                                                           '29'
            
            write(unit= file_unit_radiate, fmt= descr_fmt) '%', &    
                                                           'scan|', &
                                                           'mjd|', &
                                                           'year|', &
                                                           'doy|', &
                                                           'h|', &
                                                           'm|', &
                                                           'sec|', &
                                                           'station|', &
                                                           'azimuth [rad]|', &
                                                           'calc. elev. [rad]|', &
                                                           'source|', &
                                                           'T [°C]|', &
                                                           'p [hPa]|', &
                                                           'w[hPa]|', &
                                                           'ztd [m]|', &
                                                           'zhd [m]|', &
                                                           'zwd [m]|', &
                                                           'std [m]|', &
                                                           'shd [m]|', &
                                                           'swd [m]|', &
                                                           'e at st[rad]|', &
                                                           'e rtrd [rad]|', &
                                                           'bend. [m]|', &
                                                           'total mf|', &
                                                           'hydr. mf|', &
                                                           'wet mf|', &
                                                           'T [°C]|', &
                                                           'p [hPa]|', &
                                                           'w[hPa]'
            write(unit= file_unit_radiate, fmt= '(a)') '%'
            
            ! write observation data and ray-tracing results for each observation to the file with specified format
            ! note: Implied loop for writing each line is used as it should be faster than a normal do loop.
            write(unit= file_unit_radiate, fmt= data_fmt) ( rdlog_epoch % observations(ind_obs) % scannr, &
                                                            rdlog_epoch % observations(ind_obs) % mjd, &
                                                            rdlog_epoch % observations(ind_obs) % year, &
                                                            rdlog_epoch % observations(ind_obs) % doy, &
                                                            rdlog_epoch % observations(ind_obs) % hour, &
                                                            rdlog_epoch % observations(ind_obs) % min, &
                                                            rdlog_epoch % observations(ind_obs) % sec, &
                                                            rdlog_epoch % observations(ind_obs) % station, &
                                                            rdlog_epoch % observations(ind_obs) % az, &
                                                            rdlog_epoch % observations(ind_obs) % elev, &
                                                            rdlog_epoch % observations(ind_obs) % source, &
                                                            rdlog_epoch % observations(ind_obs) % temp, &
                                                            rdlog_epoch % observations(ind_obs) % pres, &
                                                            rdlog_epoch % observations(ind_obs) % wvpr, &
                                                            rdlog_epoch % delay(ind_obs) % dz_total, &
                                                            rdlog_epoch % delay(ind_obs) % dz_h, &
                                                            rdlog_epoch % delay(ind_obs) % dz_w, &
                                                            rdlog_epoch % delay(ind_obs) % ds_total_geom, &
                                                            rdlog_epoch % delay(ind_obs) % ds_h_geom, &
                                                            rdlog_epoch % delay(ind_obs) % ds_w, &
                                                            rdlog_epoch % delay(ind_obs) % e_stat, &
                                                            rdlog_epoch % delay(ind_obs) % e_outgoing_rt, &
                                                            rdlog_epoch % delay(ind_obs) % dgeo, &
                                                            rdlog_epoch % delay(ind_obs) % mf_total_geom, &
                                                            rdlog_epoch % delay(ind_obs) % mf_h_geom, &
                                                            rdlog_epoch % delay(ind_obs) % mf_w, &
                                                            T_stat_out(ind_obs), &
                                                            rdlog_epoch % meteo_stat_out(ind_obs) % p, &
                                                            rdlog_epoch % meteo_stat_out(ind_obs) % wvpr, ind_obs= 1, nr_obs )
            
        end if
        
        ! close the created .radiate-file
        ! note: close() for an unit that has not been opened, e.g. due to an error, has no effect, e.g. as file could not be created
        close(unit= file_unit_radiate)
        
    end subroutine create_radiate_global
    
end module module_create_radiate_global