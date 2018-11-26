! module_qsort.f90 

!****************************************************************************
!
! PURPOSE:  Module inheriting subroutines for sorting either double precision, or integer
!           values or character strings in ascending (alphabetical order).
!           The used sorting algorithm quicksort is implemented from
!           program code provided by rosettacode.org.
!           Quicksort is an unstable sorting algorithm. This means a previous order is not preserved.
!
!           Via an integrated interface the correct subroutine is selected
!           when calling "qsort". Depending on the dummy arguments of the
!           individual subroutines this is decided.
!
!
! INPUT and OUTPUT of the call "qsort(A, nlist)"
!
! INPUT:
!         A...   structure containing vector "value" for sorting and "order" for sorting indices
!         nlist... size of vector which will be sorted
! 
! 
! OUTPUT:
!         A...   same as input, but sorted using the values vector:
!                   % order... indices after sorting usabel for values_sorted= values_unsorted(order)
!                   % value... sorted version of input for "values"
!
!----------------------------------------------------------------------------
! 
! Fortran-file created by Armin Hofmeister
!   based on program code provided by rosettacode.org (quicksort-algorithm)
! 
!----------------------------------------------------------------------------
! History:
! 
! 06.11.2014: create the Fortran-file
! 11.11.2014: programming
! 12.11.2014: programming, adding sorting subroutine for integer
! 13.11.2014: comments
! 20.04.2015: correct comments
! 21.04.2015: add comments
! 14.03.2016: comment
!
!****************************************************************************
    
    
module module_qsort
    
    ! Define modules to be used
    use module_sort_type_definitions ! module for type definitions used in sorting subroutines
    
    
    ! declaration part
    implicit none
    public :: qsort
    private :: qsort_integer_ascend, qsort_double_ascend, qsort_char_alpha
    
    
    ! interface for decision of correct subroutine
    interface qsort
        module procedure qsort_integer_ascend, qsort_double_ascend, qsort_char_alpha
    end interface
    
    
contains
    
    
    !----------------------------------------------------------------------------------------------
    ! subroutine for sorting integer values
    !----------------------------------------------------------------------------------------------
    
    ! recursive subroutine for sorting integer values in ascending order
    recursive subroutine qsort_integer_ascend(A, nlist)
        
        ! declaration part
        implicit none
    
        ! dummy arguments
        integer, intent(in) :: nlist
        type(sort_integer_type), dimension(nlist), intent(in out) :: A ! integer containing type
        
        ! local variables
        integer :: left, right
        real :: random ! random number
        integer :: pivot ! pivot element of same type as A % value
        type(sort_integer_type) :: temp ! temporary struture of same type as "A"
        integer :: marker
        
        ! Body of the subroutine
        
        ! if the list has more than one entry
        if (nlist > 1) then
            
            ! get a random number [0, 1[
            call random_number(random)
            
            ! determine the pivot element using the random number
            pivot= A(int(random*real(nlist-1))+1) % value ! random pivot (not best performance, but avoids worst-case scenario of number of operations for sorting)
            
            left = 0
            right = nlist + 1
            
            do while (left < right)
                
                right = right - 1
                
                do while (A(right) % value > pivot)
                    right = right - 1
                end do
                
                left = left + 1
                
                do while (A(left) % value < pivot)
                    left = left + 1
                end do
                
                if (left < right) then
                    
                    temp= A(left)
                    A(left) = A(right)
                    A(right) = temp
                    
                end if
                
            end do
            
            if (left == right) then
                marker = left + 1
            else
                marker = left
            end if
            
            call qsort_integer_ascend(A(:marker-1), marker-1)
            call qsort_integer_ascend(A(marker:), nlist-marker+1)
            
        end if
        
    end subroutine qsort_integer_ascend
    
    
    !----------------------------------------------------------------------------------------------
    ! subroutine for sorting double precision values
    !----------------------------------------------------------------------------------------------
    
    ! recursive subroutine for sorting double precision values in ascending order
    recursive subroutine qsort_double_ascend(A, nlist)
        
        ! declaration part
        implicit none
    
        ! dummy arguments
        integer, intent(in) :: nlist
        type(sort_double_type), dimension(nlist), intent(in out) :: A ! double containing type
        
        ! local variables
        integer :: left, right
        real :: random ! random number
        double precision :: pivot ! pivot element of same type as A % value
        type(sort_double_type) :: temp ! temporary struture of same type as "A"
        integer :: marker
        
        ! Body of the subroutine
        
        if (nlist > 1) then
            
            call random_number(random)
            pivot= A(int(random*real(nlist-1))+1) % value ! random pivot (not best performance, but avoids worst-case scenario of number of operations for sorting)
            left = 0
            right = nlist + 1
            
            do while (left < right)
                
                right = right - 1
                
                do while (A(right) % value > pivot)
                    right = right - 1
                end do
                
                left = left + 1
                
                do while (A(left) % value < pivot)
                    left = left + 1
                end do
                
                if (left < right) then
                    
                    temp= A(left)
                    A(left) = A(right)
                    A(right) = temp
                    
                end if
                
            end do
            
            if (left == right) then
                marker = left + 1
            else
                marker = left
            end if
            
            call qsort_double_ascend(A(:marker-1), marker-1)
            call qsort_double_ascend(A(marker:), nlist-marker+1)
            
        end if
        
    end subroutine qsort_double_ascend
    
    
    !----------------------------------------------------------------------------------------------
    ! subroutine for sorting character strings
    !----------------------------------------------------------------------------------------------
    
    ! recursive subroutine for sorting character strings alphabetical
    recursive subroutine qsort_char_alpha(A, nlist)
        
        ! declaration part
        implicit none
    
        ! dummy arguments
        integer, intent(in) :: nlist
        type(sort_char_type), dimension(nlist), intent(in out) :: A ! character string containing type
        
        ! local variables
        integer :: left, right
        real :: random ! random number
        character(len=:), allocatable :: pivot ! pivot element of same type as A % value
        type(sort_char_type) :: temp ! temporary struture of same type as "A"
        integer :: marker
        
        ! Body of the subroutine
        
        if (nlist > 1) then
            
            call random_number(random)
            pivot= A(int(random*real(nlist-1))+1) % value ! random pivot (not best performance, but avoids worst-case scenario of number of operations for sorting)
            left = 0
            right = nlist + 1
            
            do while (left < right)
                
                right = right - 1
                
                do while (A(right) % value > pivot)
                    right = right - 1
                end do
                
                left = left + 1
                
                do while (A(left) % value < pivot)
                    left = left + 1
                end do
                
                if (left < right) then
                    
                    temp= A(left)
                    A(left) = A(right)
                    A(right) = temp
                    
                end if
                
            end do
            
            if (left == right) then
                marker = left + 1
            else
                marker = left
            end if
            
            call qsort_char_alpha(A(:marker-1), marker-1)
            call qsort_char_alpha(A(marker:), nlist-marker+1)
            
        end if
        
    end subroutine qsort_char_alpha
    
end module module_qsort