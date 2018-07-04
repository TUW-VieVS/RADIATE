! module_date_type_definition.f90 

!****************************************************************************
!
! PURPOSE:  Module for defining type (structure) used in the mjd2date and date2mjd conversions
!           subroutines:
!                          mjd2date
!                          date2mjd
!
!----------------------------------------------------------------------------
! 
! Fortran-file created by Armin Hofmeister
! 
!----------------------------------------------------------------------------
! History:
! 
! 13.11.2014: create the Fortran-file
! 03.12.2014: comments
!
!****************************************************************************
    
module module_date_type_definition

    ! define type for date representation
    type, public :: date_type
        integer :: year
        integer :: month
        integer :: day
        integer :: hour
        integer :: min
        double precision :: sec
    end type date_type
    
end module module_date_type_definition