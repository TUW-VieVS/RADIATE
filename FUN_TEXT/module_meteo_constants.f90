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
! Changes by Janina Boisits:
! 07.03.2018: add constants for optical refractivity computation
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
    
    ! standard pressure and temperature of dry air
    ! see Mendes & Pavlis 'High-accuracy zenith delay prediction at optical wavelengths'
    double precision, parameter :: Pd = 1013.25d0 ! in [hPa]
    double precision, parameter :: Td_k = 288.15d0 ! in [K]
    double precision, parameter :: td_deg = 15d0 ! in [°C]
    
    ! standart pressure and temperature of water vapour
    ! see Mendes & Pavlis 'High-accuracy zenith delay prediction at optical wavelengths'
    double precision, parameter :: Pw = 13.33d0 ! in [hPa]
    double precision, parameter :: Tw_k = 293.15d0 ! in [K]
    double precision, parameter :: tw_deg = 20d0 ! in [°C]
    
    ! coefficients for computing compressibility factors
    ! see Mendes & Pavlis 'High-accuracy zenith delay prediction at optical wavelengths'
    double precision, parameter :: a0 = 1.58123d-4 ! in [K/hPa]
    double precision, parameter :: a1 = -2.9331d-6 ! in [1/hPa]
    double precision, parameter :: a2 = 1.1043d-8 ! in [K/(K*hPa)]
    double precision, parameter :: b0 = 5.707d-4 ! in [K/hPa]
    double precision, parameter :: b1 = -2.051d-6 ! in [1/hPa]
    double precision, parameter :: c0 = 1.9898d-2 ! in [K/hPa]
    double precision, parameter :: c1 = -2.376d-4 ! in [1/hPa]
    double precision, parameter :: d0 = 1.83d-7 ! in [K^2/hPa^2]
    double precision, parameter :: e0 = -0.765d-4 ! in [K^2/hPa^2]
    
    ! standard compressibility factor of dry air and water vapour
    ! see Mendes & Pavlis 'High-accuracy zenith delay prediction at optical wavelengths'
    double precision, parameter :: Zd = 1 - (Pd / Td_k) * (a0 + a1*td_deg + a2*td_deg**2) + (Pd / Td_k)**2 * d0
    double precision, parameter :: Zw = 1 - (Pw / Tw_k) * (a0 + a1*tw_deg + a2*tw_deg**2 + (b0 + b1*tw_deg) + (c0 + c1*tw_deg)) + (Pw / Tw_k)**2 * (d0 + e0)
    
    ! standard density of water vapour
    ! see Mendes & Pavlis 'High-accuracy zenith delay prediction at optical wavelengths'
    double precision, parameter :: rho_ws = (Pw * Mw) / (Zw * R * Tw_k) ! in [(hPa*kg)/J]
    
    ! coefficients of group refractive index of dry air
    ! see Mendes & Pavlis 'High-accuracy zenith delay prediction at optical wavelengths'
    double precision, parameter :: k0_opt = 238.0185d0 ! in [mym^-2]
    double precision, parameter :: k1_opt = 5792105d0 ! in [mym^-2]
    double precision, parameter :: k2_opt = 57.362d0 ! in [mym^-2]
    double precision, parameter :: k3_opt = 167917d0 ! in [mym^-2]
    double precision, parameter :: xc = 375d0 ! IAG recommendations
    double precision, parameter :: Cco2 = 1 + 0.534d0 * 1d-6 * (xc - 450d0)
    !double precision, parameter :: Ngaxs = 1d-2 * (k1 * (k0+sigma**2) / (k0-sigma**2)**2 + k3 * (k2+sigma**2) / (k2-sigma**2)**2) * Cco2
    
    ! coefficients of group refractive index of water vapour
    ! see Mendes & Pavlis 'High-accuracy zenith delay prediction at optical wavelengths'
    double precision, parameter :: w0 = 295.235d0 ! in []
    double precision, parameter :: w1 = 2.6422d0 ! in [mym^2]
    double precision, parameter :: w2 = -0.032380d0 ! in [mym^4]
    double precision, parameter :: w3 = 0.004028d0 ! in [mym^6]
    double precision, parameter :: cf = 1.022d0
    !double precision, parameter :: Ngws = 1d-2 * cf * (w0 + 3 *w1*sigma**2 + 5*w2*sigma**4 + 7*w3*sigma**6)


end module module_meteo_constants