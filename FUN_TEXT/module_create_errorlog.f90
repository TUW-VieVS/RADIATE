! module_create_errorlog.f90
    

!****************************************************************************
!
! PURPOSE:  Module containing subroutine for the creation of an error-logfile.
!
!           This subroutine carries creates an .err-file containing information about different errors that occured during
!           the processing of the ray-tracing.
!
!           Note: Created as module to avoid explicit inferface for the subroutine as needed because of assumed-shape dummy argument in subroutine.
!
!
! INPUT:
!        errorlog.... structure containing information about occured errors
!                     in the following variables:
!        
!           % error_nr...... variable for storing the error number
!           % error_type.... variable for storing the error type definition
!           % error_descr... variable for storing the error description
!
!        save_path.............. specifies path to the saving directory for the errorlog
!        sessionname............ name of the session for which the errorlog is created
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
! 30.04.2015: create the Fortran-file based on the Matlab-file "create_errorlog.m"
! 04.05.2015: programming
! 05.05.2015: comments, add info about grib-file resolutions
! 06.05.2015: add comments
! 07.05.2015: add comments, programming
! 13.05.2015: programming
! 03.06.2015: change to implied loop for for writing the data lines
! 11.06.2015: correct comments
!
!****************************************************************************

module module_create_errorlog

contains
    
    subroutine create_errorlog( errorlog, &
                                save_path, &
                                sessionname, &
                                RADIATE_version, &
                                RADIATE_subversion, &
                                epolog_mode, &
                                interpolation_method, &
                                raytr_method, &
                                filename_stat_info, &
                                epoch_name, &
                                gridres )
    
        ! Define modules to be used
        use module_type_definitions, only: errorlog_type, gridres_type
        
    
        !----------------------------------------------------------------------------
    
        ! Variable definitions
        implicit none
    
    
        ! dummy arguments
        !----------------
    
        ! INPUT
    
        ! define input variable for the error-structure
        ! note: dimension is defined through the number of occurred errors
        type(errorlog_type), dimension(:), intent(in) :: errorlog
        
        ! define input variable for the sessionname
        character(len=*), intent(in) :: sessionname
        
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
        
        ! define variable for storing the total number of errors
        integer :: nr_err
        
        ! variable for storing the unit of the file
        integer, parameter :: file_unit_err= 1
    
        ! variable for storing the opening status of the file
        integer :: open_status
        
        ! define variable for storing the current date
        character(len=8) :: date_value
        character(len=10) :: date_formatted
    
        ! define variable for storing the current time
        character(len=10) :: time_value
        character(len=12) :: time_formatted
        
        ! define loop variable
        integer :: ind_epoch
        
        ! define loop variable
        integer :: ind_err
        
        ! define format for writing the errors
        ! notes on output format:
        ! 1x outputs one blank
        ! ss supresses plus sign output for all following formats
        character(len=*), parameter :: data_fmt='(a1, 1x, ss, i5, 1x, a1, 1x, a, 1x, a1, 1x, a)' ! note: ss supresses plus sign output for all following formats
        
        
        !----------------------------------------------------------------------------
        
        !============================================================================
        ! Body of the subroutine
        !============================================================================
        
        ! Create error-logfile
        ! The error-logfile contains information about errors that occured during the processing.
        
        ! get number of errors that have occured
        nr_err= size(errorlog)
        
        ! check if any error has occured
        ! note: this check has also been done in the subroutine "RayTrace_main_global" to check wether it is necessary to create an error-logfile at all
        if ( nr_err > 0 ) then
            
            ! get the current date and time
            call DATE_AND_TIME(date= date_value, time= time_value)
    
            ! format the date and time
            date_formatted= date_value(1:4) // '-' // date_value(5:6) // '-' // date_value(7:8)
            time_formatted= time_value(1:2) // ':' // time_value(3:4) // ':' // time_value(5:)
            
            ! create new .err-ascii-file for storing the errorlog
            ! note: an old file is replaced by the new one as status is set to 'replace'
            open(unit= file_unit_err, file= save_path // sessionname // '.err', action= 'write', status= 'replace', iostat= open_status)
            
            ! Check if textfile can't be opened (iostat /= 0 --> error)
            if (open_status /= 0) then
                ! report warning message
                write(unit= *, fmt= '(tr4, a, a, a, /)') 'Warning: Problem with creating error-logfile: "', save_path // sessionname // '.err"', '! No file created!'
                
            else
                
                ! report process
                ! notes on output format:
                ! i0.1 outputs the minimal field width with at least one digit printed
                write(unit= *, fmt= '(tr4, a, i0.1, a)') 'Error-logfile for ', nr_err, ' encountered errors will be created.'
                
                
                ! write to errorlog-file
                
                ! write header with information about session-name, start-time of ray-tracing, parameterization and total number of errors
                write(unit= file_unit_err, fmt= '(a)') '%******************************************************************************************************************'
                
                ! "S" for session name
                ! notes on output format:
                ! 1x outputs one blank
                ! first character in the line containing the session name is set to "S" as to faciliate the finding of the session name later when the file is read in
                write(unit= file_unit_err, fmt= '(a)') '% error-logfile for session:'
                write(unit= file_unit_err, fmt= '(a1, 1x, a)') 'S', sessionname
                
                ! "N" for number of errors
                ! notes on output format:
                ! 1x outputs one blank
                ! i0.1 outputs the minimal field width with at least one digit printed
                ! ss supresses plus sign output for all following formats
                ! first character in the line containing the number of observations set to "N" as to faciliate the finding of the number later when the file is read in
                write(unit= file_unit_err, fmt= '(a)') '% total number of errors in the ray-tracing of this session:'
                write(unit= file_unit_err, fmt= '(a1, 1x, ss, i0.1)') 'N', nr_err
                
                ! "P" for processing information about start time of processing and some parameterizations
                ! notes on output format:
                ! 1x outputs one blank
                ! i0.1 outputs the minimal field width with at least one digit printed
                ! ss supresses plus sign output for all following formats
                ! first character in the line containing the number of observations set to "N" as to faciliate the finding of the number later when the file is read in
                write(unit= file_unit_err, fmt= '(a)') '% information about the processing:'
                ! write time information about creation of the file
                write(unit= file_unit_err, fmt= '(a)') '% date and time of .err-file-creation:'
                write(unit= file_unit_err, fmt= '(a1, 1x, a, 1x, a1, 1x, a)') 'P', date_formatted, '|', time_formatted
                ! write time information about processing parameterization
                write(unit= file_unit_err, fmt= '(a)') '% parameterization of the processing:'
                write(unit= file_unit_err, fmt= '(a1, 1x, a, 5( 1x, a1, 1x, a ))') 'P', &
                                                                                    RADIATE_version, '|', &
                                                                                    RADIATE_subversion, '|', &
                                                                                    epolog_mode, '|', &
                                                                                    interpolation_method, '|', &
                                                                                    raytr_method, '|', &
                                                                                    filename_stat_info
                ! write information about the epochs and the horizontal resolution used during the processing
                write(unit= file_unit_err, fmt= '(a)') '% Grib-file epoch(s) and resolutions:'
                
                ! write entries for all epochs
                ! note: size of "gridres" determines the number of epochs = gribfiles
                ! note: Implied loop for writing each line is used as it should be faster than a normal do loop.
                ! notes on output format:
                ! 1x outputs one blank
                ! i0.1 outputs the minimal field width with at least one digit printed
                ! ss supresses plus sign output for all following formats
                ! f5.3 --> form x.yyy
                write(unit= file_unit_err, fmt= '(ss, a1, 1x, a, i0.1, a, a, 1x, a1, 1x, a, f5.3, a, f5.3, a)') ( 'P', &
                                                                                            'Epoch ', ind_epoch, ': ', epoch_name(ind_epoch), '|', &
                                                                                            'Grib-file resolution: res_lat= ', gridres(ind_epoch) % res_lat, &
                                                                                            ' deg, res_lon= ', gridres(ind_epoch) % res_lon, ' deg', &
                                                                                            ind_epoch= 1, size(gridres) )
                
                ! write column description
                write(unit= file_unit_err, fmt= '(a)') '%'
                write(unit= file_unit_err, fmt= '(a)') '% Columns:'
                ! notes on output format:
                ! tr4 tabs 4 positions to the right
                write(unit= file_unit_err, fmt= '(a, tr4, a)') '%', '0  .... data specifier sign (S... session, N... number of errors, P... processing details, D... error data)'
                write(unit= file_unit_err, fmt= '(a, tr4, a)') '%', '1  .... error number'
                write(unit= file_unit_err, fmt= '(a, tr4, a)') '%', '2  .... error type definition'
                write(unit= file_unit_err, fmt= '(a, tr4, a)') '%', '3  .... error description'
                write(unit= file_unit_err, fmt= '(a)') '%******************************************************************************************************************'
                write(unit= file_unit_err, fmt= '(a)') '%'
                
                ! write errors to the file
                ! note: Implied loop for writing each error entry is used as it should be faster than a normal do loop.
                !       Implicit writing of all structure variables is not supported
                !       due to allocatable attribute of some variables in the structure.
                    
                ! "D" for data of error
                write(unit= file_unit_err, fmt= data_fmt) ('D', &
                                                            errorlog(ind_err) % error_nr, &
                                                            '|', &
                                                            errorlog(ind_err) % error_type, &
                                                            '|', &
                                                            errorlog(ind_err) % error_descr, ind_err= 1, nr_err)
                
            end if
        
            ! close the created textfile
            ! note: close() for an unit that has not been opened, e.g. due to an error, has no effect, e.g. as file could not be created
            close(unit= file_unit_err)
            
        ! in case of no error entry in "errorlog"
        else
            
            ! report process
            write(unit= *, fmt= '(tr4, a)') 'No errors occurred during ray-tracing of the session. Therefore no .err-file is created. A possibly existing .err-file of a former processing run will be deleted.'
            
            ! delete a possible former .err-file for the same session
            ! note: in case of an old processing with occurred errors an .err-file that still exists needs to be deleted
            ! try to open a possibly existing old .err-file for the same session
            open(unit= file_unit_err, file= save_path // sessionname // '.err', status= 'old', iostat= open_status)
            
            ! check if there is an old .err-file
            ! note: if "iostat" is not 0 then there is no old file with the same filename
            if (open_status == 0) then
                
                ! close the file again by using "delete" as status resulting in a dletion of the old file
                ! note: "delete" only works if the file has not been opened as open(readonly)
                close(unit= file_unit_err, status= 'delete')
                
            end if
            
        end if ! end of if in case at least one error has occurred during the processing
        
    end subroutine create_errorlog
    
end module module_create_errorlog