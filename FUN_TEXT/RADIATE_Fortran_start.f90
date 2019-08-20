! RADIATE_Fortran_start.f90 
!
! FUNCTIONS:
! RADIATE_Fortran_start - Entry point of console application.
!

!****************************************************************************
!
! PROGRAM: RADIATE_Fortran_start
!
! PURPOSE:  Entry point for the console application.
!   
!           Ray-tracing program for one input AZEL-file (one session)
!
!               
!
!
! CALL: .exe with command line parameters:
!           1) AZEL-filename
!           2) Station catalogue-filename
!           3) additional parameters
!    
!    
!           extra) possible: output of file containing the print command infos, e.g. >Diary.txt
!                  Note: This parameter will not be counted as command argument!
!
!----------------------------------------------------------------------------
! 
! Fortran-file created by Armin Hofmeister
!   based on Matlab scripts created by Armin Hofmeister
! 
!----------------------------------------------------------------------------
! History:
! 
! 03.11.2014: create the Fortran-file based on the Matlab-file "RADIATE_global_start.m"
! 04.11.2014: programming
! 05.11.2014: programming
! 06.11.2014: programming
! 02.12.2014: test epolog mode detection
! 20.01.2015: test interpolation mode detection
! 28.01.2015: test interpolation mode detection
! 03.02.2015: comments
! 19.02.2015: programming: set program modes via command-line arguments
! 26.02.2015: add message that the arguments in the command-line are case-sensitive
! 20.04.2015: correct comments and messages
! 06.05.2015: add argument for command-line to specify if session index- and epolog-files
!             should be deleted after their usage (after loading)
!             improve error messaging
! 31.08.2015: change \ to / in path definitions for Windows and Linux compatibility
!             delete unused variable declaration "command_line"
! 10.09.2015: add comments
!             change directory names to upper case only

! Changes by Daniel Landskron:
! 19.11.2017: also enabled for individual internal creation of uniform azel-files by adding a 9th input argument
!
!Changes by Janina Boisits:
! 26.02.2018: add option 'microwave' or 'optical' to parameters
! 08.03.2018: add option 'wavelength' to info message
!
!****************************************************************************

program RADIATE_Fortran_start

    ! Define modules to be used
    
    ! Module for type definitions
    use module_type_definitions, only : parameters_type
    
    !----------------------------------------------------------------------------
    
    ! Variable definitions
    implicit none
    
    ! define loop variable
    integer :: i
    
    
    
    ! define structure for parameter settings for the calculation process
    type(parameters_type) :: parameters
    
    ! define variable for storing the number of command line arguments
    integer :: nr_args
    
    ! define variable for number of mandatory input arguments in command-line
    integer, parameter :: nr_mand_arg= 2
    
    ! define variable for (temporarily) storing the command arguments
    ! note: len must be specified as allocatable does not work in call of get_command_argument()
    character(len=1000) :: arg
    
    ! define variable for number of erroneous additional arguments
    integer :: arg_error
    
    
    !============================================================================
    ! Body of the program
    !============================================================================
    
    
    ! Get command-line arguments invoking the program
    
    ! determine the number of command-line arguments
    ! note: command-line arguments for input and output are not recognized, e.g. >../DATA/DIARY/Diary.txt is not recognized as argument, but work!
    ! determine the number of arguments
    nr_args=command_argument_count() ! attention: output files like >file.out don't count as an argument
    
    ! initialize variable for check on erroneous additional arguments
    arg_error= 0
    
    ! at least two arguments must be present representing the AZEL file name and the station data file
    if (nr_args < nr_mand_arg) then
        ! write error message and stop the program
        write(unit= *, fmt= '(a)') 'Error: Not enough input arguments provided! First argument in command must be the AZEL filename and second argument must be the station data filename! Program stopped!'
        write(unit= *, fmt= '(tr4, a)') 'First argument in command must be the AZEL filename!'
        write(unit= *, fmt= '(tr4, a)') 'Second argument in command must be the station data filename!'
        write(unit= *, fmt= '(a, /)') 'Program stopped!'
        ! stop the program
        stop
        
    ! in case at least two arguments are specified    
    else
        
        ! get the filename of the AZEL-file from the command-line input argument
        ! Attention: This must be the first argument!
        call get_command_argument(1, arg) ! attention: passing the argument directly to the structure field is not successful probably because of allocatable type of target variable
        !call getarg(1, arg) ! alternative method
        
        
        ! store the filename in the target variable
        ! note: use trim() to remove blanks as "arg" is initialized too long
        parameters % load_filename_azel= trim(arg)
        
        
        ! get the filename of the station data-file from the command-line input argument
        ! Attention: This must be the second argument!
        call get_command_argument(2, arg) ! attention: passing the argument directly to the structure field is not successful probably because of allocatable type of target variable
        
        ! store the filename in the target variable
        ! note: use trim() to remove blanks as "arg" is initialized too long
        parameters % load_filename_stat_info= trim(arg)
        
        
        ! Define default ray-tracing options
        ! Note: In case that not all possible settings are made through the user by specifying command-line arguments.
        !       If more than the mandatory arguments are set the default values for the paramters may be overwritten by provided command-line arguments!
        
        !***************************************************************************************************
        !**************************** Define DEFAULT PARAMETERS*********************************************
    
        ! set default vertical interpolation method
        ! options: 'gridwise', 'profilewise'
        parameters % interpolation_method= 'profilewise'
    
        ! set default ray-tracing method (approach for defining intersection points and delays)
        ! options: 'pwl'
        ! not realized for Fortran: 'ref_pwl', 'Thayer'
        parameters % raytr_method= 'pwl'
    
        ! set deafult time interpolation mode = mode of assigning the observations to epochs (select method of creating epologs)
        ! options: 'one_epoch_per_obs': description: each observation is assigned to only the one grib-file epoch where the obs. is in the interval [epoch-3h, epoch+3h[
        !          'two_epochs_per_obs': description: each observation is assigned to two grib-file epochs; note: also observations at the exact epoch time are
        !                                             covered twice, although the second epolog entry will not alter the time interpolation result!
        parameters % epolog_mode= 'two_epochs_per_obs'
    
        ! set default if .trp-file should be created
        ! .TRUE. .... create trp-file for each session
        ! .FALSE. ... don't create trp-files
        parameters % create_trp= .TRUE.
    
        ! set default if the error log should be saved
        ! .TRUE. .... save reported errors in a .err-file
        ! .FALSE. ... don't save reported errors in a .err-file
        parameters % save_errorlog= .TRUE.
        
        ! set default if a clean up of created session index- and epolog-files should be done
        ! .TRUE. .... clean up of created session index- and epolog-files
        ! .FALSE. ... no clean up of created session index- and epolog-files
        parameters % cleanup= .FALSE.
        
        ! set default if an Azel file shall be read, or created individually and uniformly
        ! .TRUE. .... existing azel file is read
        ! .FALSE. ... a uniform azel file is created for individual settings stored in INPUT/AZEL_spec.txt
        parameters % readAzel= .TRUE.
    
        ! set default wavelength range
        ! options: 'microwave': description: refractivity field computed for microwave observations
        !          'optical': description: refractivity field computed for optical observations
        parameters % wavelength= 'microwave'
    
        !***************************************************************************************************
        !***************************************************************************************************
        
        
        ! check if more than the mandatory two arguments are passed
        ! Note: The previously defined default settings are overwritten in case the specific parameter is found in the command-line!
        !       If a specific mode is specified more than once, the last argument is the one that is taken into account!
        if (nr_args > nr_mand_arg) then
            
            ! loop over all remaining command-line arguments
            do i= nr_mand_arg + 1, nr_args
        
                ! get specific command-line argument
                call get_command_argument(i, arg) ! attention: passing the argumentv directly to the structure field is not successful probably because of allocatable type of target variable
                !call getarg(i, arg) ! alternative method
        
                ! determine what argument it is
                ! note: use trim() to remove blanks as arg is initialized too long
                ! attention: the arguments are case-sensitive
                select case (trim(arg))
                    
                    ! possible cases for vertical interpolation method
                    ! options: 'profilewise', 'gridwise'
                    
                    ! case of "profilewise"
                    case ('-profilewise')
                        ! set vertical interpolation method
                        parameters % interpolation_method= 'profilewise'
                        
                    ! case of "gridwise"
                    case ('-gridwise')
                        ! set vertical interpolation method
                        parameters % interpolation_method= 'gridwise'
                        
                    !----------------------------------------------------------
                    
                    ! possible cases for ray-tracing method
                    ! options: 'pwl', 'ref_pwl', 'Thayer'
                    
                    ! case of "pwl"
                    case ('-pwl')
                        ! set ray-tracing method
                        parameters % raytr_method= 'pwl'
                        
                    ! case of "ref_pwl"
                    case ('-ref_pwl')
                        ! set ray-tracing method
                        parameters % raytr_method= 'ref_pwl'
                        
                    ! case of "Thayer"
                    case ('-Thayer')
                        ! set ray-tracing method
                        parameters % raytr_method= 'Thayer'
                        
                    !----------------------------------------------------------
                    
                    ! possible cases for time interpolation mode = mode of assigning the observations to epochs (select method of creating epologs)
                    ! options: 'one_epoch_per_obs': description: each observation is assigned to only the one grib-file epoch where the obs. is in the interval [epoch-3h, epoch+3h[
                    !          'two_epochs_per_obs': description: each observation is assigned to two grib-file epochs; note: also observations at the exact epoch time are
                    !                                covered twice, although the second epolog entry will not alter the time interpolation result!
                    
                    ! case of "one_epoch_per_obs"
                    case ('-one_epoch_per_obs')
                        ! set time interpolation mode = mode of assigning the observations to epochs
                        parameters % epolog_mode= 'one_epoch_per_obs'
                        
                    ! case of "two_epochs_per_obs"
                    case ('-two_epochs_per_obs')
                        ! set time interpolation mode = mode of assigning the observations to epochs
                        parameters % epolog_mode= 'two_epochs_per_obs'
                    
                    !----------------------------------------------------------
                    
                    ! possible cases for trp-creation
                    ! options: 'trp', 'notrp'
                    
                    ! case of "trp"
                    case ('-trp')
                        ! set mode for trp-creation
                        parameters % create_trp= .TRUE.
                        
                    ! case of "notrp"
                    case ('-notrp')
                        ! set mode for trp-creation
                        parameters % create_trp= .FALSE.
                        
                    !----------------------------------------------------------
                    
                    ! possible cases for errorlog-creation
                    ! options: 'errorlog', 'noerrorlog'
                    
                    ! case of "errorlog"
                    case ('-errorlog')
                        ! set mode for errorlog-creation
                        parameters % save_errorlog= .TRUE.
                        
                    ! case of "noerrorlog"
                    case ('-noerrorlog')
                        ! set mode for errorlog-creation
                        parameters % save_errorlog= .FALSE.
                    
                    !----------------------------------------------------------    
                    
                    ! possible cases for cleaning up created session index- and epolog-files
                    ! options: 'cleanup', 'nocleanup'
                    
                    ! case of "errorlog"
                    case ('-cleanup')
                        ! set mode for cleaning up created session index- and epolog-files
                        parameters % cleanup= .TRUE.
                        
                    ! case of "noerrorlog"
                    case ('-nocleanup')
                        ! set mode for cleaning up created session index- and epolog-files
                        parameters % cleanup= .FALSE.
                    
                    !----------------------------------------------------------
                                        
                    ! possible cases for reading an azel file or creating it individually
                    ! options: 'readAzel', '-createUniAzel'
                    
                    ! case of "errorlog"
                    case ('-readAzel')
                        ! set mode for cleaning up created session index- and epolog-files
                        parameters % readAzel= .TRUE.
                        
                    ! case of "noerrorlog"
                    case ('-createUniAzel')
                        ! set mode for cleaning up created session index- and epolog-files
                        parameters % readAzel= .FALSE.
                    
                    !----------------------------------------------------------
                                        
                    ! possible cases for wavelength range of observations
                    ! options: '-microwave', '-optical'
                    
                    ! case of "microwave"
                    case ('-microwave')
                        ! set mode for computing refractivity field for microwave observations
                        parameters % wavelength= 'microwave'
                        
                    ! case of "optical"
                    case ('-optical')
                        ! set mode for computing refractivity field for optical observations
                        parameters % wavelength= 'optical'
                    
                    !----------------------------------------------------------
                    
                    ! case of unknown argument
                    case default ! case if other cases fail
                        
                        ! raise number of erroneous additional arguments
                        arg_error= arg_error + 1
                    
                        ! report error message of unknown argument
                        ! notes on output format:
                        ! i0.1 outputs the minimal field width with at least one digit printed
                        write(unit= *, fmt= '(a, i0.1, a, a)') 'Error ', arg_error,': Unknown argument: ', trim(arg)
                        
                end select
                    
            end do ! end do-loop over additional arguments
            
            ! check if any of the additional argument options were defined wrong
            if ( arg_error > 0 ) then
                
                ! report number of erroneous additional arguments
                ! notes on output format:
                        ! i0.1 outputs the minimal field width with at least one digit printed
                write(unit= *, fmt= '(/, a, i0.1, a, /)') 'Encountered ', arg_error, ' erroneous additional argument(s)!'
            
                ! report about case sensitive arguments
                write(unit= *, fmt= '(a)') 'Attention:'
                write(unit= *, fmt= '(a, /)') 'The arguments are case-sensitive, so be sure to use the right cases!'
                        
                ! report info messages about possible arguments
                write(unit= *, fmt= '(a)') 'Note:'
                write(unit= *, fmt= '(a)') 'The first two arguments are mandatory and their input order is fixed, but the additional arguments can be passed in any order.'
                write(unit= *, fmt= '(tr4, a)') 'First argument in command must be the AZEL filename!'
                write(unit= *, fmt= '(tr4, a)') 'Second argument in command must be the station data filename!'
                write(unit= *, fmt= '(a)') 'Multiple arguments for the same program mode will lead to the usage of the last specified option for the specific mode.'
                write(unit= *, fmt= '(a, /)') 'In case one or more additional arguments are not specified then the default is taken for the unspecified mode.'
                        
                write(unit= *, fmt= '(/, a)') 'List of possible additional arguments (choose one per category):'
                        
                write(unit= *, fmt= '(/, tr4, a)') 'Options for vertical interpolation method:'
                write(unit= *, fmt= '(tr8, a)') 'a) "-profilewise" (Default)'
                write(unit= *, fmt= '(tr8, a)') 'b) "-gridwise"'
                        
                write(unit= *, fmt= '(/, tr4, a)') 'Options for ray-tracing method:'
                write(unit= *, fmt= '(tr8, a)') 'a) "-pwl" (Default)'
                write(unit= *, fmt= '(tr8, a)') 'b) "-ref_pwl"'
                write(unit= *, fmt= '(tr8, a)') 'c) "-Thayer"'
                        
                write(unit= *, fmt= '(/, tr4, a)') 'Options for time interpolation mode = mode of assigning the observations to epochs:'
                write(unit= *, fmt= '(tr8, a)') 'a) "-one_epoch_per_obs"'
                write(unit= *, fmt= '(tr8, a)') 'b) "-two_epochs_per_obs" (Default)'
                        
                write(unit= *, fmt= '(/, tr4, a)') 'Options for wavelength range of observations:'
                write(unit= *, fmt= '(tr8, a)') 'a) "-microwave" (Default)'
                write(unit= *, fmt= '(tr8, a)') 'b) "-optical"'
                        
                write(unit= *, fmt= '(/, tr4, a)') 'Options for trp-file creation:'
                write(unit= *, fmt= '(tr8, a)') 'a) "-trp" (Default)'
                write(unit= *, fmt= '(tr8, a)') 'b) "-notrp"'
                        
                write(unit= *, fmt= '(/, tr4, a)') 'Options for errorlog-file creation:'
                write(unit= *, fmt= '(tr8, a)') 'a) "-errorlog" (Default)'
                write(unit= *, fmt= '(tr8, a)') 'b) "-noerrorlog"'
                        
                write(unit= *, fmt= '(/, tr4, a)') 'Options for cleaning up created session index- and epolog-files:'
                write(unit= *, fmt= '(tr8, a)') 'a) "-cleanup"'
                write(unit= *, fmt= '(tr8, a)') 'b) "-nocleanup" (Default)'
                
                write(unit= *, fmt= '(/, tr4, a)') 'Options for reading azel file or individually creating them:'
                write(unit= *, fmt= '(tr8, a)') 'a) "-readAzel" (Default)'
                write(unit= *, fmt= '(tr8, a)') 'b) "-createUniAzel"'
                        
                write(unit= *, fmt= '(/, a, /)') 'Please redefine the argument(s). Program stopped!'
                ! stop the program
                stop
                
            end if ! end if any additional argument is wrong
            
        end if ! end if more than the mandatory two arguments are passed
    
    end if ! end if check at least mandatory arguments are passed
    
    !---------------------------------------------------------------------------------------------------
    
    
    ! Define paths to the directories containing the specific files and specify filename for station data file
    
    !***************************************************************************************************
    !**************************** Define PARAMETERS ****************************************************
    
    ! Choose path and filename, which should be used to load the VLBI station data.
    ! Supported files for the function *import_VLBI_stations* is "vlbi.ell" or textfiles with equivalent
    ! layout.
    parameters % load_path_stat_info= '../DATA/STATIONS/'
    
    ! define path to input-files (AZEL)
    parameters % load_path_azel= '../DATA/AZEL/'
    
    ! define path to the needed grib-files
    parameters % load_path_grib= '../DATA/GRIB/'
    
    ! define path to the directory where the files with the geoid undulations in different
    ! latitude/longitude resolution are stored
    parameters % load_path_undulation='../DATA/UNDULATIONS/';
    
    ! define path and filename of the textfile containing the specifications for the uniform azel file creation
	if ( parameters % load_filename_stat_info == 'slr.ell' ) then
        parameters % load_path_filename_indAzelSpec= '../DATA/INPUT/AZEL_spec_SLR.txt'
	else
		parameters % load_path_filename_indAzelSpec= '../DATA/INPUT/AZEL_spec.txt'
	end if
    
    
    !***************************************************************************************************
    !***************************************************************************************************
    
    !---------------------------------------------------------------------------------------------------
    
    ! Process
    
    ! call the subroutine *RayTrace_main* to do the ray-tracing processing
    call RayTrace_main_global( parameters )
    
    
    
end program RADIATE_Fortran_start
