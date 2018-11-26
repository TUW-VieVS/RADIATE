! module_meshgrid2D.f90

!****************************************************************************
!
! PURPOSE:  Module inheriting subroutines of "meshgrid2D" for different classes of variables (double precision, ...) to create a cartesian grid in 2D space.
!            
!           "meshgrid2D" is a subroutine to create a cartesian grid in 2D space using the input of two
!           vectors specifying the values of the grid for each dimension.
!
!           Via an integrated interface the correct subroutine is selected
!           when calling "meshgrid2D". Depending on the dummy arguments of the
!           individual subroutines this is decided.
!
!           Adapted from MATLAB meshgrid description:
!           meshgrid2D(x_vec,y_vec, x_grid, y_grid) replicates the grid vectors xgv and ygv to produce the coordinates of a rectangular grid (X, Y).
!           The grid vector x_vec is replicated size(y_vec) times to form the columns (as whole vector) of X.
!           The grid vector y_vec is replicated size(x_vec) times to form the rows (as whole vector) of Y.
!
!
! INPUT:
!         x_vec... input vector for x-values
!         y_vec... input vector for y-values
! 
! 
! OUTPUT:
!         x_grid... grid of x-values
!         y_grid... grid of y-values
!
!----------------------------------------------------------------------------
! 
! Fortran-file created by Armin Hofmeister
!   based on Matlab script
! 
!----------------------------------------------------------------------------
! History:
! 
! 12.01.2015: create the Fortran-file based on the Matlab-intrinsic subroutine "meshgrid.m"
! 13.01.2015: programming
! 20.01.2015: comments
! 26.01.2015: comments
!
!****************************************************************************

    
module module_meshgrid2D
    
    ! Define modules to be used
    
    
    ! declaration part
    implicit none
    public :: meshgrid2D
    private :: meshgrid2D_double
    
    
    ! interface for decision of correct subroutine
    interface meshgrid2D
        module procedure meshgrid2D_double
    end interface
    
    
contains
    
    subroutine meshgrid2D_double( x_vec, &
                                  y_vec, &
                                  x_grid, &
                                  y_grid )

        ! Define modules to be used
    
        !----------------------------------------------------------------------------
    
        ! Variable definitions
        implicit none
    
    
        ! dummy arguments
        !----------------
    
        ! INPUT
    
        ! define variable for storing the input vector for x-values
        double precision, dimension(:), intent(in) :: x_vec
    
        ! define variable for storing the input vector for y-values
        double precision, dimension(:), intent(in) :: y_vec
    
    
        ! OUTPUT
    
        ! Notes: Dimension specification at initialization of the variable does not work as variables get correct sizes, but undefined addresses of memory location!
        !        Solution: allocate output variables in calling function!
        !        This behaviour may be due to the fact that there is no explicit interface and/or because the subroutine is part of a module.
        
        ! define variable for storing grid of x-values
        double precision, dimension(:, :), intent(out) :: x_grid
    
        ! define variable for storing grid of x-values
        double precision, dimension(:, :), intent(out) :: y_grid
        
    
        ! local variables
        !----------------
        
        
        !----------------------------------------------------------------------------
    
        !============================================================================
        ! Body of the subroutine
        !============================================================================
        
        
        ! create the grid of the x-values
        ! note: spread adds a defined number of copies of an array to a defined dimension 
        x_grid= spread(x_vec, 1, size(y_vec)) ! note: dim= 1 to spread as new columns
        
        ! create the grid of the y-values
        y_grid= spread(y_vec, 2, size(x_vec)) ! note: dim= 2 to spread as new rows
        
        
    end subroutine meshgrid2D_double
                                  
end module module_meshgrid2D