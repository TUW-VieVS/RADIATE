! module_date2mjd.f90

!****************************************************************************
!
! PURPOSE:  Module containing elemental subroutine to convert civil date to modified Julian date
!           (elemental procedure enabling array input).
!
!           "date2mjd(date, mjd)" converts the given civil dates (in form of 
!           a "date" structure of dimension(n) containing year, month, day, hour, min, sec) to
!           modified Julian date(s) mjd dimension(n).
!    
!           Note: If the year is lower than 100, but year >= 80 it is assummed to be 19yy.
!                 If the year is lower than 100, but year < 80 it is assummed to be 20yy.
!
!           Attention: Input date-values will not be checked for correctness plausibility.
!                      Thus, be sure to avoid input of negative or not possible values as this leads to wrong mjd!
! 
!
! INPUT:
!         date...... structure (date_type) containing the converted mjd:
!                           year, month, day, hour and min as integers and
!                           sec as double precision value
! 
! OUTPUT:
!         mjd... modified julian date
!
!
!----------------------------------------------------------------------------
! 
! Fortran-file created by Armin Hofmeister
!   based on Matlab scripts created by AW (Version of 2007-01-22)
! 
!----------------------------------------------------------------------------
! History:
! 
! 15.01.2015: create the Fortran-file based on the Matlab-file "date2mjd.m"
! 19.01.2015: programming
! 20.01.2015: comments
! 29.01.2015: change to module as to save explicit interface block
!
!****************************************************************************

module module_date2mjd
    
contains
    
    elemental subroutine date2mjd(date, mjd)
    
        ! Define modules to be used
        use module_date_type_definition
    
    
        !----------------------------------------------------------------------------
    
        ! Variable definitions
        implicit none
    
    
        ! dummy arguments
        !----------------
    
        ! INPUT
        ! input variable for date
        type(date_type), intent(in) :: date
            
        ! OUTPUT
        ! output variable for mjd
        double precision, intent(out) :: mjd
    
    
        ! local arguments
        !----------------
    
        ! define date structure variable to store temporal date values which can be redefined
        ! as input variable for date is only intent(in) and can not be redefined if necessary!
        type(date_type) :: temp_date
    
        ! define variable for storing the fraction of day (unit = fractional day)
        double precision :: frcday
    
        ! define variable for storing the Julian date
        double precision :: jd
    
        !----------------------------------------------------------------------------
    
        !============================================================================
        ! Body of the subroutine
        !============================================================================
    
        ! assign input date-variable to temporal variable, whose values are not fixed to the input
        ! as it doesn't have intent(in) option restricting redefinement of values
        ! note: in the course of the subroutine to calculate the mjd some values of the input date may have to be redefined
        temp_date= date
    
        ! check if year is < 100
        if (temp_date % year < 100) then
            ! check if year is < 80
            if (temp_date % year < 80) then
                ! assume it is meant 20xx for a year < 80
                temp_date % year = 2000 + temp_date % year
            else
                ! assume it is meant 19xx for 80<= year < 100
                temp_date % year = 1900 + temp_date % year
            end if
        end if
    
    
        ! compute the mjd
    
        ! modify the month and year to fit the mjd-calculation formula
        ! check if the month is <= 2
        if (temp_date % month <= 2) then
            ! reduce year by 1
            temp_date % year = temp_date % year - 1
            ! raise month by 12
            temp_date % month = temp_date % month + 12
        end if
    
        ! compute fraction of day
        frcday = ( (temp_date % sec / 60 + temp_date % min) / 60 + temp_date % hour) / 24
    
        ! convert to Julian date
        ! note: aint() truncates a value to a whole number, this means in other words rounding towards zero (like fix in MATLAB).
        !       Without specifying a kind value aint() preserves the input kind, but input must be a kind of real!
        !       aint() is different from anint(), which rounds towards the nearest integer!
        jd = aint(365.25d0 * temp_date % year) + aint(30.6001d0 * (temp_date % month + 1)) + &
                  temp_date % day + frcday + 1720981.5d0
    
        ! convert to modified Julian day
        mjd = jd - 2400000.5d0
    
    end subroutine date2mjd
    
end module module_date2mjd
