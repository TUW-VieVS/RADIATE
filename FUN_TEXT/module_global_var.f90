! module_global_var.f90 

!****************************************************************************
!
! PURPOSE:  Module for defining global variables
!
!----------------------------------------------------------------------------
! 
! Fortran-file created by Armin Hofmeister
! 
!----------------------------------------------------------------------------
! History:
! 
! 03.12.2014: create the Fortran-file
! 08.01.2015: add variables for start time and rate for getting elapsed
! 29.01.2015: comments
! 20.04.2015: correct comments
!
!****************************************************************************
    
module module_global_var

    ! Variable definitions
    implicit none
    
    save ! stellt sicher, dass der Inhalt in den Speicherplaetzen zwischen den einzelnen Einbindevorgaengen in den einzelnen Programmeinheiten unveraendert bleibt,
         ! d.h. der Wert beim nächsten Aufruf in einer Subroutine ist jener wie zuletzt die Subroutine abgeschlossen wurde und die Werte gehen nicht verloren.
         ! Aber Achtung: Die Werte können während der Subroutine verändert werden!
    
    ! define variables for system clock values used for timing of the program
    ! variables for the start time count and count rate and intermediate count number needed for "get_elapsed_time" subroutine whenever called and values will only be assigned in the first assignment 
    integer(8) :: icount_start, count_rate, icount_interm ! use kind=8 to ensure that count_max is not reached and therefore resetting of the count to zero is avoided
    ! variable for storing the calculated elapsed time
    real :: elapsed_time
    
    ! define variable for storing the gridded undulation data
    double precision, dimension(:,:), allocatable :: undulation_grid_global
    
end module module_global_var