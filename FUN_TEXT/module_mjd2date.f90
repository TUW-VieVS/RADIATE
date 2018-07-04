! module_mjd2date.f90

!****************************************************************************
!
! PURPOSE:  Module containing elemental subroutine to convert modified Julian date to civil date
!           (elemental procedure enabling array input).
!
!           "mjd2date(mjd, date)" converts the given modified Julian date(s) mjd size(n)
!           to civil dates returned as "date" structure size(n) containing year, month, day, hour, min, sec.
!           
!           Note, the accuracy of this function is 10^-4 sec if mjd is not
!           an integer vector. For high timing accuracy, treat the fractions of a day
!           separately i.e., keep an integer vector mjd and a double vector sid
!           (seconds wihtin day).
!           
!           Note: The algorithm has been taken from Hofmann-Wellenhof et al. (1991)
!                 and is valid only for epochs between Mar-1-1900 and Feb-28-2100.
!            
! 
!
! INPUT:
!         mjd... modified julian date
! 
! 
! OUTPUT:
!         date...... structure (date_type) containing the converted mjd:
!                           year, month, day, hour and min as integers and
!                           sec as double precision value
!
!----------------------------------------------------------------------------
! 
! Fortran-file created by Armin Hofmeister
!   based on Matlab scripts created by AW (2007-01-16 )
! 
!----------------------------------------------------------------------------
! History:
! 
! 13.11.2014: create the Fortran-file based on the Matlab-file "mjd2date.m"
! 26.11.2014: comments
! 20.01.2015: comments
! 29.01.2015: change to module as to save explicit interface block
! 11.05.2015: change mod() to modulo()
!
!****************************************************************************

module module_mjd2date
    
contains
    
    elemental subroutine mjd2date(mjd, date)
    
        ! Define modules to be used
        use module_date_type_definition
    
    
        !----------------------------------------------------------------------------
    
        ! Variable definitions
        implicit none
    
    
        ! dummy arguments
        !----------------
    
        ! input variable for mjd
        double precision, intent(in) :: mjd
    
        ! output variable for date
        type(date_type), intent(out) :: date
    
    
        ! local arguments
        !----------------
    
        ! define variable for storing the seconds within day
        double precision :: sid
    
    
        ! define auxiliary variables
        double precision :: mjd_floor
        double precision :: a, b, c, d, e
        double precision :: year, month, day, hour, min, sec
    
        !----------------------------------------------------------------------------
    
        !============================================================================
        ! Body of the subroutine
        !============================================================================
    
        ! Determine the sid
        ! Check if input mjd has no .0 fraction part, use 1.0 for modulo() in other to fit the data type
        ! note: modulo() in Fortran works like mod() in MATLAB (if non-zero result then the sign of the second argument is retained)
        if (modulo(mjd, 1.0d0) == 0) then
            sid = 0
            mjd_floor = floor(mjd)
        else
            mjd_floor = floor(mjd) ! mjd of each day's 00:00:00
            sid = (mjd - mjd_floor)*86400 ! seconds within day   
        end if
    
        ! helper quantities
        a = mjd_floor + 2400001 
        b = a + 1537
        c = floor( (b - 122.1d0) / 365.25d0 )
        d = floor( 365.25d0 * c )
        e = floor( (b - d) / 30.6001d0 )
    
        ! assign year, month and day to output structure
        ! note: int() is not necessary, but clarifies the assignment
        day = b - d - floor(30.6001d0 * e) ! days
        month = e - 1 - 12 * floor(e / 14.0d0) ! months
        year = c - 4715 - floor((7 + month) / 10) ! years
    
    
        ! determine values of hour, min, sec
    
        ! in case an integer or real mjd without fractional part was the input
        if (sid == 0) then
           hour = 0
           min = 0
           sec = 0
        ! in case there was a fractional part in the input mjd
        else
            hour  = floor(sid / 3600)
            min = floor((sid - hour * 3600) / 60)
            sec = sid - (hour * 60 + min) * 60
        
            ! round seconds to 10^-4 secs
            sec = anint(sec * 1d4) * 1d-4 ! note: anint() rounds to nearest integer and returns a real value of same kind as input, but input must be of a kind of real; use d for exponent to preserve double precision
        
            ! check if rounding lead to an increase of the seconds to 60.xxx
            if (sec >= 60) then
                sec = 0 ! set them to 0
                min = min + 1 ! increase minutes by 1
            
                ! check if now the minutes have been raised to 60
                if (min == 60) then
                    min = 0 ! set them 0
                    hour = hour + 1 ! increase hours by 1
                
                    ! raise of day, month, year ...... not necessary as algorithm calculated sid as seconds within day, so day should be correct as initially calculated
                end if
            
            end if

        end if
    
    
        ! assign year, month and day to output structure
        ! attention: here the values will be converted to integer for year, month, day, hour, min (but: these vales should not have a fractional part other than .0)
        ! note: int() is not necessary, but clarifies the assignment
        date % year = int(year)
        date % month = int(month)
        date % day = int(day)
        date % hour = int(hour)
        date % min = int(min)
    
        date % sec = sec
    
    end subroutine mjd2date
    
end module module_mjd2date