! module_sort_type_definitions.f90 

!****************************************************************************
!
! PURPOSE:  Module for defining types (structures) used in the sorting modules:
!           module_qsort
!           module_msort
!
!----------------------------------------------------------------------------
! 
! Fortran-file created by Armin Hofmeister
! 
!----------------------------------------------------------------------------
! History:
! 
! 12.11.2014: create the Fortran-file
! 13.11.2014: comments
!
!****************************************************************************
    
module module_sort_type_definitions

    ! define type for sorting of integer values
    type, public :: sort_integer_type
        integer :: order ! original order of unsorted data
        integer :: value ! values to be sorted
    end type sort_integer_type
    
    ! define type for sorting of double precision values
    type, public :: sort_double_type
        integer :: order ! original order of unsorted data
        double precision :: value ! values to be sorted
    end type sort_double_type
    
    ! define type for sorting of character strings
    type, public :: sort_char_type
        integer :: order ! original order of unsorted data
        character(len=:), allocatable :: value ! values to be sorted
    end type sort_char_type

end module module_sort_type_definitions