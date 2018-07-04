! module_gph2horth.f90
    

!****************************************************************************
!
! PURPOSE:  Module containing elemental subroutine to converts geopotential (=dynamic) height [m] to orthometric height [m]
!           (elemental procedure enabling array input).
!            
!           This subroutine converts geopotential (=dynamic) height [m] to orthometric height [m]
!           using the formula by meteorologist (Kraus, 2004), whose formula inherits an approximation approach
!           (in order to avoid iterative calculation).
! 
!
! INPUT:
!        gph....... geopotential (=dynamic) height(s) in [m]
!        lat_deg... latitude value in [°]
! 
! OUTPUT:
!        h_orth... orthometric height in [m]
!
!----------------------------------------------------------------------------
! 
! Fortran-file created by Armin Hofmeister
!   based on Matlab script "gph2horth.m" created by Armin Hofmeister, which is an adapted version the Matlab script "gpm2z.m" created by Johannes Böhm
! 
!----------------------------------------------------------------------------
! History:
! 
! 21.01.2015: create the Fortran-file based on the Matlab-file "gph2horth.m"
! 29.01.2015: change to module as to save explicit interface block
! 12.02.2015: use constant from mudule for "deg2rad"
! 11.05.2015: add comments
! 02.06.2015: comments
! 01.10.2015: comments
!
!****************************************************************************

module module_gph2horth

contains
    
    elemental subroutine gph2horth( gph, &
                                    lat_deg, &
                                    h_orth )

        ! Define modules to be used
        use module_constants, only: deg2rad
    
    
        !----------------------------------------------------------------------------
    
        ! Variable definitions
        implicit none
    
    
        ! dummy arguments
        !----------------
    
        ! INPUT
    
        ! variable for input of geopotential (=dynamic) height(s) in [m]
        double precision, intent(in) :: gph
    
        ! variable for input of latitude value in [°]
        double precision, intent(in) :: lat_deg
    
    
        ! OUTPUT
    
        ! variable for output of orthometric height in [m]
        double precision, intent(out) :: h_orth
    
    
        ! local variables
        !----------------
    
        ! variable for storing latitude value in [rad]
        double precision :: lat_rad
        
        
        ! CONSTANTS
        
        ! see "module_constants" for "deg2rad"
    
    
        !----------------------------------------------------------------------------
    
        !============================================================================
        ! Body of the subroutine
        !============================================================================
    
        ! transform latitude values from [°] to [rad]
        ! note: see "module_constants" for constant "deg2rad"
        ! note: extra variable for storing latitude in [rad] is needed as intent(in) variable can not be redefined
        lat_rad= lat_deg * deg2rad
    
        ! transformation of geopotential (=dynamic) height [m] to orthometric height [m]
        ! using formula by meteorologists (kraus)
        ! phi in [rad]
        h_orth = 1 / (2 * 1.57d-7) - sqrt( ( 1 / (2 * 1.57d-7) )**2 - gph / ( 1 - 0.0026373d0 * cos( 2 * lat_rad ) + 0.0000059d0 * (cos( 2 * lat_rad ) )**2 ) / 1.57d-7 )
        
    end subroutine gph2horth
    
end module module_gph2horth