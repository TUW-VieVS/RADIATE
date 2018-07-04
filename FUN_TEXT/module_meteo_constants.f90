! module_meteo_constants.f90

!****************************************************************************
!
! PURPOSE:  Module for the definition of meteorological constants.
!
!----------------------------------------------------------------------------
! 
! Fortran-file created by Armin Hofmeister
! 
!----------------------------------------------------------------------------
! History:
! 
! 29.01.2015: create the Fortran-file
! 11.06.2015: comments
!
!****************************************************************************

module module_meteo_constants
    
    ! Variable definitions
    implicit none
    
    save ! Ensures that the content of the memory spaces between different inclusions to different program units remains unchanged, i.e. the value of the
         ! variable in the next call will be the same as when the subroutine was left at the last call and the value does not get lost after the subroutine.
         ! Attention: The value of a variable can be changed during the subroutine!
    
    
    ! CONSTANTS
    
    ! Define constants and variables concerning meteorological fundamentals
    
    ! universal gas constant
    ! see scriptum Atmospheric Effects in Geodesy 2012, equation (1.23) on page 7
    double precision, parameter :: R = 8314.510d0 ! in [kg*m^2/(kmol*K*s^2)] = [J/(kmol*K)]
    
    ! molar weight of water vapour and dry constituents
    ! see scriptum Atmospheric Effects in Geodesy 2012, equations (1.35) and (1.36) on page 8
    double precision, parameter :: Mw = 18.01528d0 ! in [g/mol] = [kg/kmol]
    double precision, parameter :: Md = 28.9644d0 ! in [g/mol] = [kg/kmol]
    
    ! specific gas constants
    ! see scriptum Atmospheric Effects in Geodesy 2012, page 7 between eq. (1.27) and (1.28)
    double precision, parameter :: Rw = R/Mw ! in [J/(K*kg)] = [m^2/(K*s^2)]
    double precision, parameter :: Rd = R/Md ! in [J/(K*kg)] = [m^2/(K*s^2)]
    
    ! refractivity coefficients
    ! "Rueger - best average" values
    ! see scriptum Atmospheric Effects in Geodesy 2012, equation (2.6) on page 20
    double precision, parameter :: k1 = 77.6890d0 ! in [K/hPa]
    double precision, parameter :: k2 = 71.2952d0 ! in [K/hPa]
    double precision, parameter :: k3 = 375463d0 ! in [K^2/hPa]
    
    ! see scriptum Atmospheric Effects in Geodesy 2012, equation (2.22) on page 25
    double precision, parameter :: k2s = k2 - k1 * (Mw / Md) ! in [K/hPa]

end module module_meteo_constants