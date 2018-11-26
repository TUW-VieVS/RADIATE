! module_lagrange_int.f90
    

!****************************************************************************
!
! PURPOSE:  Module containing subroutine for lagrange interpolation.
!
!           This subroutine carries out a lagrange interpolation using given values y at locations x to
!           calculate the value y_int at the desired location x_int.
!           This version supports the ouput of the Lagrange basis polynomials L. With them it is possible to calculate
!           for other y-values (of other domains) the y_int at the same x_int.
!
!           Advantage of the output of the Lagrange basis polynomials L: It is not necessary to calculate the L again if the x_int location stays the same!
!
!           Input values of x need not to be sorted. They must just be in the order of according y-values.
!
!           Equations taken from scriptum "Mathematische Methoden der Geowissenschaften" (version 2008) pages 2-3, equations (1.5 - 1.6).
!
!           Note: Created as module to avoid explicit inferface for the subroutine as needed because of assumed-shape dummy arguments in subroutine.
!
!
! INPUT:
!        x_int... location where y_int should be determined;
!                 x_int must be a scalar
!        x....... vector containing the locations for which values y are known;
!                 x must contain at least 2 locations
!        y....... vector containing the values at the locations x;
!                 length of y must be equal to length of x
!
!
! OUTPUT:
!        y_int... determined value at the location x_int
!
!         optional output: in case output parameters are present at the call:
!
!         optional: L_out... Lagrange basis polynomials for the interpolation of y_int at x_int
!                            note: Lagrange basis polynomials are independent from y-values and only depend on the x-values,
!                                  so they can be used for any y-values at the same x-values in order to calculate y_int.
!
!----------------------------------------------------------------------------
! 
! Fortran-file created by Armin Hofmeister
!   based on Matlab scripts created by Armin Hofmeister
! 
!----------------------------------------------------------------------------
! History:
! 
! 27.04.2015: create the Fortran-file based on the Matlab-file "lagrange_int_mult_y.m"
! 30.04.2015: comments
! 13.05.2015: correct comments
! 13.01.2016: correct comments
!
!****************************************************************************

module module_lagrange_int

contains
    
    subroutine lagrange_int( x_int, &
                             x, &
                             y, &
                             y_int, &
                             L_out ) ! optional
    
        ! Define modules to be used
        
    
        !----------------------------------------------------------------------------
    
        ! Variable definitions
        implicit none
    
    
        ! dummy arguments
        !----------------
    
        ! INPUT
    
        ! define input variable for storing location where y_int should be determined
        ! note: x_int must be a scalar
        double precision, intent(in) :: x_int
        
        ! define input variable for storing vector containing the locations for which values y are known
        ! note: x must contain at least 2 locations --> dimension(:)
        double precision, dimension(:), intent(in) :: x
        
        ! define input variable for storing vector containing the values at the locations x;
        ! note: The the size of y must be equal to the size of x!
        double precision, dimension(:), intent(in) :: y
        
        
        ! OUTPUT
    
        ! define output variable for storing the determined value at the location x_int
        ! note: y_int is a scalar
        double precision, intent(out) :: y_int
        
        ! OPTIONAL OUTPUT
        
        ! optional output variable for the Lagrange basis polynomials
        double precision, dimension(:), intent(out), optional :: L_out
        
    
        ! local variables
        !----------------
        
        ! define variable for storing the polynomial order = number of known pairs of (x, y)
        integer :: n
        
        ! define variable for storing the Lagrange basis polynomials
        double precision, dimension(:), allocatable :: L
        
        ! define variable for loop index
        integer :: i, j
        
        !----------------------------------------------------------------------------
        
        !============================================================================
        ! Body of the subroutine
        !============================================================================


        ! Check if x_input is a scalar
        ! note: This check can be omitted as through definition the variable x_int must is a scalar!
        
        
        ! Define polynomial order
        
        ! get number of known values
        ! note: size() without dimension specification returns the total number total elements in the array and
        !       x should only have one dimension
        n= size(x)
        
        
        ! Check if number of elements in x fit to number of according elements in y
        ! note: As x and y have only one dimension it is sufficient to check using size() without dimension specification
        if ( n /= size(y) ) then
            
            ! report error message
            write(unit= *, fmt= '(a)') 'Error in lagrange_int_mult_y: Number of elements in y must be equal to number of elements in x! Program stopped!'
            ! stop the program
            stop
            
        end if
        
        
        ! Lagrange Interpolation
        
        ! if x and y are scalars no lagrange interpolation can be done
        if (n < 2) then
            
            ! report error message
            write(unit= *, fmt= '(a)') 'Error in lagrange_int_mult_y: Lagrange interpolation is not supported if number of known pairs of values is smaller than 2! Program stopped!'
            ! stop the program
            stop
            
        else
            
            ! Calculate Lagrange basis polynomials for each order
            
            ! calculate terms for each order
            
            ! allocate Lagrange basis polynomials for each order
            allocate( L(n) )
            
            ! initialize Lagrange basis polynomials for each order
            !note: initialize all L entries with ones due to the following multiplication for each order [*(x_int-x(j))/(x(i)-x(j))]
            L= 1
        
            ! loop over polynomial order to create each Lagrange basis polynomials
            do i= 1, n
                
                ! loop over polynomial order to sequentially create one complete Lagrange basis polynomial
                do j= 1, n
                    
                    ! only calculate the new step if indices of i and j are different
                    if ( i /= j) then
                        L(i)= L(i) * ( x_int - x(j) ) / ( x(i) - x(j) )
                    end if
                    
                end do
                
            end do
            
            ! calculate interpolated value y_int at x_int
            ! note: dot_product() does the following: sum(vector_a * vector_b), where vector_a * vector_b is the element-wise multiplication
            y_int= dot_product(L, y)
            
            ! determine if optional output arguments are needed
            if ( present(L_out) ) then
                L_out= L
            end if
            
        end if
            
    end subroutine lagrange_int
    
end module module_lagrange_int