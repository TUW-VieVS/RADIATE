! module_mjd2doy.f90

!****************************************************************************
!
! PURPOSE:  Module containing elemental subroutine to convert modified Julian date to doy.
!
! 
!
! INPUT:
!         mjd...... modified Julian date
! 
! OUTPUT:
!         doy...... day of year
!
!
!----------------------------------------------------------------------------
! 
! Fortran-file created by Daniel Landskron on the basis of mjd2date.m from VieVS
! 
!----------------------------------------------------------------------------
! History:
! 
! 19.11.2017: create the Fortran-file
!
!****************************************************************************

module module_mjd2doy
    
contains
    
    elemental subroutine mjd2doy(mjd, doy)
    
    
        !----------------------------------------------------------------------------
    
        ! Variable definitions
        implicit none
    
    
        ! dummy arguments
        !----------------
    
        ! INPUT
        ! input variable for date
        double precision, intent(in) :: mjd
            
        ! OUTPUT
        ! output variable for mjd
        integer, intent(out) :: doy
    
    
        ! local arguments
        !----------------
    
        ! argument needed for the conversion
        integer :: mjd_2000
        integer :: mjd_int
        integer :: ic0
        double precision :: ic
        double precision :: ir
        double precision :: iy1
        
    
        !----------------------------------------------------------------------------
    
        !============================================================================
        ! Body of the subroutine
        !============================================================================
    
        mjd_2000 = 51544  ! mjd for year=2000

        mjd_int = int(mjd)

        ic0 = mjd_int-mjd_2000
        ic = ic0/1461;
        
        ic = int(ic)
        
        if (ic0<0) then
            ic = ic-1
        end if

        ir = ic0-ic*1461

        if (ir<366) then
            doy = ir+1
        end if
 
        if (ir >= 366) then
            ir = ir-366
            iy1 = ir/365
            iy1 = int(iy1)
            doy = ir-iy1*365+1
        end if
        
    end subroutine mjd2doy
    
end module module_mjd2doy