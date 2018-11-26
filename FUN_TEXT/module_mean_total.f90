! module_mean_total.f90

!****************************************************************************
!
! PURPOSE:  Module inheriting subroutines of "mean_total" for different classes of variables (double precision, ...) to calculate the total mean value of all input values of a variable.
!            
!           "mean_total" is a subroutine to calculate the arithmeticmean of all input values from a single variable.
!           Note: In case the input variable is empty an exception may follow.
!                 For double precision the result for the mean is set to NaN!
!
!           Via an integrated interface the correct subroutine is selected
!           when calling "mean_total". Depending on the dummy arguments of the
!           individual subroutines this is decided.
!
!
! INPUT:
!         a... input variable with values
! 
! 
! OUTPUT:
!         m... arithmetic mean of all values from "a"
!
!----------------------------------------------------------------------------
! 
! Fortran-file created by Armin Hofmeister
! 
!----------------------------------------------------------------------------
! History:
! 
! 26.01.2015: create the Fortran-file
!
!****************************************************************************

    
module module_mean_total2D
    
    ! Define modules to be used
    use, intrinsic :: ieee_arithmetic ! intrinsic module to e.g. specify NaN and Inf
    
    ! declaration part
    implicit none
    public :: mean_total2D
    private :: mean_total2D_double
    
    
    ! interface for decision of correct subroutine
    interface mean_total2D
        module procedure mean_total2D_double
    end interface
    
    
contains
    
    subroutine mean_total2D_double( a, &
                                    m )

        ! Define modules to be used
    
        !----------------------------------------------------------------------------
    
        ! Variable definitions
        implicit none
    
    
        ! dummy arguments
        !----------------
    
        ! INPUT
    
        ! define variable for storing the input variable used for calculating the mean
        double precision, dimension(:, :), intent(in) :: a
    
    
        ! OUTPUT
        
        ! define variable for storing the mean value
        double precision, intent(out) :: m
        
    
        ! local variables
        !----------------
        
        ! define variable for storing the size of the input variable
        integer :: s
        
        !----------------------------------------------------------------------------
    
        !============================================================================
        ! Body of the subroutine
        !============================================================================
        
        ! determine total number of elements in input variable
        s= size(a) ! note size() without dim specified, delivers total number of elements
        
        ! check if total number of elements in "a" is not 0
        if (s /= 0) then
            
            ! calculate the mean of all input values
            ! note: sum() without dimension specification delivers sum of all elements
            m= sum(a) / s
        
        ! prevent divide by 0 error and set mean to NaN
        else
            
            ! determine the NaN-value as signaling NaN (producing exceptions when used in calculation)
            m= ieee_value(0.0d0, ieee_signaling_nan) ! 0.0d0 specifies double precision kind
        
        end if
        
        
    end subroutine mean_total2D_double
                                  
end module module_mean_total2D