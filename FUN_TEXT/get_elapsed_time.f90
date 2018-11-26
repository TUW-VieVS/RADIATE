! get_elapsed_time.f90

!****************************************************************************
!
! PURPOSE:  Subroutine to calculate the elapsed time
!            
!           "get_elapsed_time" calculates the elapsed time (wall time) between an input event determined through
!           the call of "system_clock" and the call of this subroutine.
!
! 
!
! INPUT:
!         icount_start...   number of clock counts at start time
!         count_rate.....   processor count rate
! 
! 
! OUTPUT:
!         icount_current...   number of clock counts at current time (call of this subroutine)
!         elapsed_time.....   time in [s] between start time and call of this subroutine
!
!----------------------------------------------------------------------------
! 
! Fortran-file created by Armin Hofmeister
! 
!----------------------------------------------------------------------------
! History:
! 
! 10.11.2014: create the Fortran-file
! 13.11.2014: comments
!
!****************************************************************************
    
    
subroutine get_elapsed_time( icount_start, count_rate, icount_current, elapsed_time )

    ! Define modules to be used
    
    !----------------------------------------------------------------------------
    
    ! Variable definitions
    implicit none
    
    
    ! dummy arguments
    !----------------
    
    ! define variable for the start time count and count rate
    integer(8), intent(in) :: icount_start, count_rate ! use kind=8 to ensure that count_max is not reached and therefore resetting of the count to zero is avoided
    
    ! define variable for the current count count
    integer(8), intent(out) :: icount_current ! use same kind (=8) as above to ensure that count_max is not reached and therefore resetting of the count to zero is avoided
    
    ! define variable for storing the calculated elapsed time in [s]
    real, intent(out) :: elapsed_time
    
    
    !----------------------------------------------------------------------------
    
    !============================================================================
    ! Body of the subroutine
    !============================================================================


    ! get the current count
    call system_clock(icount_current)
    
    ! caculate the elapsed time in [s]
    elapsed_time= real(icount_current-icount_start) / real(count_rate)
    
    
    end subroutine get_elapsed_time