! module_msort.f90 

!****************************************************************************
!
! PURPOSE:  Module inheriting subroutines for sorting either double precision or integer
!           values or character strings in ascending (alphabetical order).
!           The used sorting algorithm mergesort is implemented from
!           program code provided by rosettacode.org.
!           Mergesort is a stable sorting algorithm. This means a previous order is preserved if possible.
!
!           Via an integrated interface the correct subroutine is selected
!           when calling "msort". Depending on the dummy arguments of the
!           individual subroutines this is decided.
!
!
! INPUT and OUTPUT of the call "msort(A, nlist)"
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
!   based on program code provided by rosettacode.org (mergesort-algorithm)
! 
!----------------------------------------------------------------------------
! History:
! 
! 12.11.2014: create the Fortran-file
! 13.11.2014: comments
! 20.04.2015: correct comments
! 21.04.2015: add comments
! 14.03.2016: comment
!
!****************************************************************************
    
    
module module_msort
    
    ! Define modules to be used
    use module_sort_type_definitions ! module for type definitions used in sorting subroutines
    
    
    ! declaration part
    implicit none
    public :: msort
    private :: merge_integer, msort_integer_ascend, merge_double, msort_double_ascend, merge_char, msort_char_alpha
    
    
    ! interface for decision of correct subroutine
    interface msort
        module procedure merge_integer, msort_integer_ascend, merge_double, msort_double_ascend, merge_char, msort_char_alpha
    end interface
   
    
contains
    
    
    !----------------------------------------------------------------------------------------------
    ! subroutines for sorting integer values
    !----------------------------------------------------------------------------------------------
    
    ! subroutine for merging lists after splitting the original list into sublists (that will be sorted)
    subroutine merge_integer(A, nA, B, nB, C, nC)
        
        ! declaration part
        implicit none
        
        ! dummy arguments
        integer, intent(in) :: nA, nB, nC ! normal usage: nA + nB = nC
        type(sort_integer_type), dimension(nA), intent(in out) :: A ! integer containing type
        type(sort_integer_type), dimension(nB), intent(in) :: B ! integer containing type, B overlays C (nA + 1:nC)
        type(sort_integer_type), dimension(nC), intent(in out) :: C ! integer containing type
        
        ! local variables
        integer :: i, j, k
        
        ! Body of the subroutine
        
        ! initialize indices
        i = 1
        j = 1
        k = 1
        
        ! determine using list A in which position the entries of the list B have to be placed in order that the merged entries of A and B result in sorted list C
        ! in this step the stability of the sorting algorithm is preserved
        ! as entries of list B will only be assigned before entries of list A if they are smaller than entries in list A
        do while (i <= nA .AND. j <= nB)
            if (A(i) % value <= B(j) % value) then
                C(k) = A(i)
                i = i + 1
            else
                C(k) = B(j)
                j = j + 1
            end if
            
            k = k + 1
            
        end do
        
        ! in case list B is fully merged, but there are still entries remaining in list A to be assigned to list C
        ! assign the remaining values of list A to the output list C for sorted output
        do while (i <= nA)
            C(k) = A(i)
            i = i + 1
            k = k + 1
        end do
        
        return
        
    end subroutine merge_integer
    
    
    ! recursive subroutine for mergesort algorithm
    ! subroutine is called recursively until one or two elements are left in the reduced sorting list so that they are sorted or just need to be swapped
    ! and then get merged again to whole list
    recursive subroutine msort_integer_ascend(A, n, T)
        
        ! declaration part
        implicit none
        
        ! dummy arguments
        integer, intent(in) :: n
        type(sort_integer_type), dimension(n), intent(in out) :: A ! integer containing type
        type(sort_integer_type), dimension((n+1)/2), intent(out) :: T ! integer containing type
        
        ! local variables
        type(sort_integer_type) :: V ! integer containing type
        integer :: nA, nB
        
        ! Body of the subroutine
        
        ! no sorting in case of one element
        if (n < 2) return
        
        ! sorting in case of two elements
        if (n == 2) then
            ! swap elements if necessary
            if (A(1) % value > A(2) % value) then
                V = A(1)
                A(1) = A(2)
                A(2) = V
            end if
            return
        end if
        
        ! half the total list that should be sorted
        nA = (n + 1) / 2 ! first half
        nB = n - nA ! second half
        
        ! sort the first half of the list
        call msort_integer_ascend(A, nA, T)
        ! sort the second half of the list
        call msort_integer_ascend(A(nA + 1), nB, T)
        
        ! merge the values merge the list-before last with the last list
        ! only done in cases where the two lists don't already have the correct order
        if (A(nA) % value > A(nA + 1) % value) then
            T(1:nA) = A(1:nA)
            call merge_integer(T, nA, A(nA + 1), nB, A, n)
        end if
        
        return
        
    end subroutine msort_integer_ascend
    
    
    !----------------------------------------------------------------------------------------------
    ! subroutines for sorting double precision values
    !----------------------------------------------------------------------------------------------
    
    ! subroutine for merging lists after splitting the original list into sublists (that will be sorted)
    subroutine merge_double(A, nA, B, nB, C, nC)
        
        ! declaration part
        implicit none
        
        ! dummy arguments
        integer, intent(in) :: nA, nB, nC ! normal usage: nA + nB = nC
        type(sort_double_type), dimension(nA), intent(in out) :: A ! double containing type
        type(sort_double_type), dimension(nB), intent(in) :: B ! double containing type, B overlays C (nA + 1:nC)
        type(sort_double_type), dimension(nC), intent(in out) :: C ! double containing type
        
        ! local variables
        integer :: i, j, k
        
        ! Body of the subroutine
        
        ! initialize indices
        i = 1
        j = 1
        k = 1
        
        ! determine using list A in which position the entries of the list B have to be placed in order that the merged entries of A and B result in sorted list C
        ! in this step the stability of the sorting algorithm is preserved
        ! as entries of list B will only be assigned before entries of list A if they are smaller than entries in list A
        do while (i <= nA .AND. j <= nB)
            if (A(i) % value <= B(j) % value) then
                C(k) = A(i)
                i = i + 1
            else
                C(k) = B(j)
                j = j + 1
            end if
            
            k = k + 1
            
        end do
        
        ! in case list B is fully merged, but there are still entries remaining in list A to be assigned to list C
        ! assign the remaining values of list A to the output list C for sorted output
        do while (i <= nA)
            C(k) = A(i)
            i = i + 1
            k = k + 1
        end do
        
        return
        
    end subroutine merge_double
    
    
    ! recursive subroutine for mergesort algorithm
    ! subroutine is called recursively until one or two elements are left in the reduced sorting list so that they are sorted or just need to be swapped
    ! and then get merged again to whole list
    recursive subroutine msort_double_ascend(A, n, T)
        
        ! declaration part
        implicit none
        
        ! dummy arguments
        integer, intent(in) :: n
        type(sort_double_type), dimension(n), intent(in out) :: A ! double containing type
        type(sort_double_type), dimension((n+1)/2), intent(out) :: T ! double containing type
        
        ! local variables
        type(sort_double_type) :: V ! double containing type
        integer :: nA, nB
        
        ! Body of the subroutine
        
        ! no sorting in case of one element
        if (n < 2) return
        
        ! sorting in case of two elements
        if (n == 2) then
            ! swap elements if necessary
            if (A(1) % value > A(2) % value) then
                V = A(1)
                A(1) = A(2)
                A(2) = V
            end if
            return
        end if
        
        ! half the total list that should be sorted
        nA = (n + 1) / 2 ! first half
        nB = n - nA ! second half
        
        ! sort the first half of the list
        call msort_double_ascend(A, nA, T)
        ! sort the second half of the list
        call msort_double_ascend(A(nA + 1), nB, T)
        
        ! merge the values merge the list-before last with the last list
        ! only done in cases where the two lists don't already have the correct order
        if (A(nA) % value > A(nA + 1) % value) then
            T(1:nA) = A(1:nA)
            call merge_double(T, nA, A(nA + 1), nB, A, n)
        end if
        
        return
        
    end subroutine msort_double_ascend
    
    
    !----------------------------------------------------------------------------------------------
    ! subroutines for sorting character strings
    !----------------------------------------------------------------------------------------------
    
    ! subroutine for merging lists after splitting the original list into sublists (that will be sorted)
    subroutine merge_char(A, nA, B, nB, C, nC)
        
        ! declaration part
        implicit none
        
        ! dummy arguments
        integer, intent(in) :: nA, nB, nC ! normal usage: nA + nB = nC
        type(sort_char_type), dimension(nA), intent(in out) :: A ! character string containing type
        type(sort_char_type), dimension(nB), intent(in) :: B ! character string containing type, B overlays C (nA + 1:nC)
        type(sort_char_type), dimension(nC), intent(in out) :: C ! character string containing type
        
        ! local variables
        integer :: i, j, k
        
        ! Body of the subroutine
        
        ! initialize indices
        i = 1
        j = 1
        k = 1
        
        ! determine using list A in which position the entries of the list B have to be placed in order that the merged entries of A and B result in sorted list C
        ! in this step the stability of the sorting algorithm is preserved
        ! as entries of list B will only be assigned before entries of list A if they are smaller than entries in list A
        do while (i <= nA .AND. j <= nB)
            if (A(i) % value <= B(j) % value) then
                C(k) = A(i)
                i = i + 1
            else
                C(k) = B(j)
                j = j + 1
            end if
            
            k = k + 1
            
        end do
        
        ! in case list B is fully merged, but there are still entries remaining in list A to be assigned to list C
        ! assign the remaining values of list A to the output list C for sorted output
        do while (i <= nA)
            C(k) = A(i)
            i = i + 1
            k = k + 1
        end do
        
        return
        
    end subroutine merge_char
    
    
    ! recursive subroutine for mergesort algorithm
    ! subroutine is called recursively until one or two elements are left in the reduced sorting list so that they are sorted or just need to be swapped
    ! and then get merged again to whole list
    recursive subroutine msort_char_alpha(A, n, T)
        
        ! declaration part
        implicit none
        
        ! dummy arguments
        integer, intent(in) :: n
        type(sort_char_type), dimension(n), intent(in out) :: A ! character string containing type
        type(sort_char_type), dimension((n+1)/2), intent(out) :: T ! character string containing type
        
        ! local variables
        type(sort_char_type) :: V ! character string containing type
        integer :: nA, nB
        
        ! Body of the subroutine
        
        ! no sorting in case of one element
        if (n < 2) return
        
        ! sorting in case of two elements
        if (n == 2) then
            ! swap elements if necessary
            if (A(1) % value > A(2) % value) then
                V = A(1)
                A(1) = A(2)
                A(2) = V
            end if
            return
        end if
        
        ! half the total list that should be sorted
        nA = (n + 1) / 2 ! first half
        nB = n - nA ! second half
        
        ! sort the first half of the list
        call msort_char_alpha(A, nA, T)
        ! sort the second half of the list
        call msort_char_alpha(A(nA + 1), nB, T)
        
        ! merge the values merge the list-before last with the last list
        ! only done in cases where the two lists don't already have the correct order
        if (A(nA) % value > A(nA + 1) % value) then
            T(1:nA) = A(1:nA)
            call merge_char(T, nA, A(nA + 1), nB, A, n)
        end if
        
        return
        
    end subroutine msort_char_alpha
    
    
end module module_msort