! module_resize_errorlog.f90

!****************************************************************************
!
! PURPOSE:  Module containing subroutine to resize an existing "errorlog"-structure.
!
!           Note: Created as module to avoid explicit inferface for the subroutine as needed because of allocatable dummy argument in subroutine.
!            
!           The output "errorlog" is a copy of the original "errorlog", but with extended size.
! 
!
! INPUT and OUTPUT:
!         errorlog....... structure of type errorlog_type, which will be resized
!
!
! INPUT:
!         new_size.... size used to resize the errorlog
!
!
!----------------------------------------------------------------------------
! 
! Fortran-file created by Armin Hofmeister
! 
!----------------------------------------------------------------------------
! History:
! 
! 09.04.2015: create the Fortran-file
! 13.05.2015: correct comments
!
!****************************************************************************


module module_resize_errorlog

contains
    
    subroutine resize_errorlog( errorlog, new_size )

        ! Define modules to be used
    
        use module_type_definitions, only: errorlog_type
    
        !----------------------------------------------------------------------------
    
        ! Variable definitions
        implicit none
    
    
        ! dummy arguments
        !----------------
    
        ! INPUT and OUTPUT
    
        ! intent(in out) for the following variable:
    
        ! define input and output variable for the array that should be resized
        type(errorlog_type), dimension(:), allocatable, intent(in out) :: errorlog
    
    
        ! INPUT
        
        ! intent(in) for the following variable:
    
        ! define input variable for storing the value for resizing the errorlog
        integer, intent(in) :: new_size
    
    
        ! local variables
        !----------------
    
        ! define temporal variable for storing the values of "errorlog"
        type(errorlog_type), dimension(:), allocatable :: temp_errorlog
    
    
    
        !----------------------------------------------------------------------------
    
        !============================================================================
        ! Body of the subroutine
        !============================================================================
    
        ! allocate the temporal variable with the new size
        allocate( temp_errorlog(new_size) )
    
        ! assign the values stored in "errorlog" to the temporal variable
        ! note: indexing (1:size(errorlog)) is necessary even if errolog is size 0 (indexing 1:0 is possible)
        temp_errorlog(1:size(errorlog))= errorlog
    
        ! deallocate the original "errorlog"
        deallocate( errorlog )
    
        ! allocate the new "errorlog" with new size
        allocate( errorlog(new_size) )
    
        ! assign the values from the temporal variable to the output variable
        errorlog= temp_errorlog
    
        ! deallocate the temporal variable
        deallocate( temp_errorlog )
    
    end subroutine resize_errorlog

end module module_resize_errorlog