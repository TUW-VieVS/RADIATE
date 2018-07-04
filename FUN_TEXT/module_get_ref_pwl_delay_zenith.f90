! module_get_ref_pwl_delay_zenith.f90


!****************************************************************************
!
! PURPOSE:  Module containing subroutine for support of refined piece-wise linear ray-tracing (calculation of zenith refractive indices and delays)
!           
!           Note: Created as module to avoid explicit inferface for the subroutine as needed because of assumed-shape arrays in subroutine.
!           
!           "get_ref_pwl_delay_zenith" is a subroutine to calculate the mean total, hydrostatic and wet zenith refractive indices and also
!           the total, hydrostatic and wet zenith delays for a specific zenith path between the last and the current intersection point.
!           
!           The equations are taken from Hobiger et al. 2008, Fast and accurate ray-tracing algorithms for
!           real-time space geodetic applications using numerical weather models, equations (30 - 31)
!				  
!
!
! INPUT and OUTPUT: 
!
!         n_m_z......... vector of mean zenith total RI calculated with n_z_lower and n_z_upper
!         n_m_h_z....... vector of mean zenith hydrostatic RI calculated with n_h_z_lower and n_h_z_upper
!         n_m_w_z....... vector of mean zenith wet RI calculated with n_w_z_lower and n_w_z_upper
!         dz_total...... vector of zenith total delay in [m] for each level
!         dz_h.......... vector of zenith hydrostatic delay in [m] for each level
!         dz_w.......... vector of zenith wet delay in [m] for each level
!
!
! INPUT:
!         n_z_lower..... value for zenith total RI of the intersection point in the lower level
!         n_h_z_lower... value for zenith hydrostatic RI of the intersection point in the lower level
!         n_w_z_lower... value for zenith wet RI of the intersection point in the lower level
!         n_z_upper..... value for zenith total RI of the intersection point in the upper level
!         n_h_z_upper... value for zenith hydrostatic RI of the intersection point in the upper level
!         n_w_z_upper... value for zenith wet RI of the intersection point in the upper level
!         level......... value for index of the current intermediate level = upper level in base data set
!         dh............ vector of zenith distance between two consecutive levels
!
!----------------------------------------------------------------------------
! 
! Fortran-file created by Armin Hofmeister
!   based on Fortran scripts created by Armin Hofmeister
! 
!----------------------------------------------------------------------------
! History:
! 
! 07.09.2015: create the Fortran-file based on the Fortran-file "module_get_ref_pwl_delay.f90"
!             adapt subroutine for zenith delay calculations only
! 09.09.2015: correct comments
!
!****************************************************************************

module module_get_ref_pwl_delay_zenith

contains
    
    subroutine get_ref_pwl_delay_zenith( n_m_z, n_m_h_z, n_m_w_z, & ! whole vectors
                                         dz_total, dz_h, dz_w, & ! whole vectors
                                         n_z_lower, n_h_z_lower, n_w_z_lower, & ! scalars
                                         n_z_upper, n_h_z_upper, n_w_z_upper, & ! scalars
                                         level, & ! scalar
                                         dh ) ! whole vector
    
    
    
        
        ! Define modules to be used
        use module_constants, only: epsilon_refr_ind
    
    
        !----------------------------------------------------------------------------
    
        ! Variable definitions
        implicit none
    
    
        ! dummy arguments
        !----------------
    
        ! INPUT and OUTPUT
    
        ! intent(in out) for the following variables:
        
        ! define input and output variables for storing the mean zenith total RI calculated with n_lower and n_upper (total hydrostatic, wet)
        ! note: size = total number of height levels
        double precision, dimension(:), intent(in out) :: n_m_z, n_m_h_z, n_m_w_z ! whole vectors
        
        ! define input and output variables for storing the zenith delays in between two levels
        ! note: size = total number of height levels -1 as the delays are in between two levels
        double precision, dimension(:), intent(in out) :: dz_total, dz_h, dz_w ! whole vectors
        
        
        ! INPUT
        
        ! intent(in) for the following variables:
        
        ! define input variables for storing zenith RI of the intersection point in the lower level (total hydrostatic, wet)
        double precision, intent(in) :: n_z_lower, n_h_z_lower, n_w_z_lower ! scalars
        
        ! define input variables for storing zenith RI of the intersection point in the upper level (total hydrostatic, wet)
        double precision, intent(in) :: n_z_upper, n_h_z_upper, n_w_z_upper ! scalars
        
        ! define input variable for storing the index of the current intermediate level = upper level in base data set
        integer, intent(in) :: level ! scalar
        
        ! define input variable for storing the zenith distance between two consecutive levels
        ! note: size = total number of height levels -1 as the zenith distance is between two consecutive levels
        double precision, dimension(:), intent(in) :: dh ! whole vector
        
        
        ! local variables
        !----------------
        
        
        ! CONSTANTS
        
        ! see "module_constants" for "epsilon_refr_ind"
        
    
        !----------------------------------------------------------------------------
    
        !============================================================================
        ! Body of the subroutine
        !============================================================================
        
        
        ! Define constant "epsilon_refr_ind" for checking if refractive index is almost equal to 1
        ! see "module_constants" for "epsilon_refr_ind"
        
        
        ! Calculate the zenith total delay
        
        ! check if refractivities n_z_lower and n_z_upper are not (almost) equal to 1 --> would cause NaN
        ! results in further calculations
        if ( (abs(n_z_lower - 1) >= epsilon_refr_ind) .AND. (abs(n_z_upper - 1) >= epsilon_refr_ind) ) then ! if true than normal calculation doesn't cause NaN results
            
            ! calculate mean total zenith refractive index at the next (="level") height (value will be used for
            ! determining the delay)
            ! see Hobiger et al. 2008, equation (30) on page 8
            n_m_z(level)= 1 + (n_z_lower - 1) * sqrt( (n_z_upper - 1) / (n_z_lower - 1) )
            
            ! calculate zenith total delay (slant or zenith, depending on value s)
            dz_total(level-1)= ( n_m_z(level) - n_m_z(level-1) ) / log( (n_m_z(level) - 1) / (n_m_z(level-1) - 1) ) * dh(level-1) ! in [m]
            
        else
            
            ! set mean zenith total refractive index at the next (="level") height to 1
            n_m_z(level)= 1
            
            ! set zenith total delay to 0
            dz_total(level-1)= 0 ! in [m]
            
        end if
        
        
        ! Calculate the zenith hydrostatic delay
        
        ! check if refractivities n_h_z_lower and n_h_z_upper are not (almost) equal to 1 --> would cause NaN
        ! results in further calculations
        if ( (abs(n_h_z_lower - 1) >= epsilon_refr_ind) .AND. (abs(n_h_z_upper - 1) >= epsilon_refr_ind) ) then ! if true than normal calculation doesn't cause NaN results
            
            ! calculate mean zenith hydrostatic refractive index at the next (="level") height (value will be used for
            ! determining the delay)
            ! see Hobiger et al. 2008, equation (30) on page 8
            n_m_h_z(level)= 1 + (n_h_z_lower - 1) * sqrt( (n_h_z_upper - 1) / (n_h_z_lower - 1) )
            
            ! calculate zenith hydrostatic delay (slant or zenith, depending on value s)
            dz_h(level-1)= ( n_m_h_z(level) - n_m_h_z(level-1) ) / log( (n_m_h_z(level) - 1) / (n_m_h_z(level-1) - 1) ) * dh(level-1) ! in [m]
            
        else
            
            ! set mean zenith hydrostatic refractive index at the next (="level") height to 1
            n_m_h_z(level)= 1
            
            ! set zenith hydrostatic delay to 0
            dz_h(level-1)= 0 ! in [m]
            
        end if
        
        
        ! Calculate the zenith wet delay
        
        ! check if refractivities n_w_z_lower and n_w_z_upper are not (almost) equal to 1 --> would cause NaN
        ! results in further calculations
        if ( (abs(n_w_z_lower - 1) >= epsilon_refr_ind) .AND. (abs(n_w_z_upper - 1) >= epsilon_refr_ind) ) then ! if true than normal calculation doesn't cause NaN results
            
            ! calculate mean zenith wet refractive index at the next (="level") height (value will be used for
            ! determining the delay)
            ! see Hobiger et al. 2008, equation (30) on page 8
            n_m_w_z(level)= 1 + (n_w_z_lower - 1) * sqrt( (n_w_z_upper - 1) / (n_w_z_lower - 1) )
            
            ! calculate zenith wet delay (slant or zenith, depending on value s)
            dz_w(level-1)= ( n_m_w_z(level) - n_m_w_z(level-1) ) / log( (n_m_w_z(level) - 1) / (n_m_w_z(level-1) - 1) ) * dh(level-1) ! in [m]
            
        else
            
            ! set mean zenith wet refractive index at the next (="level") height to 1
            n_m_w_z(level)= 1
            
            ! set zenith wet delay to 0
            dz_w(level-1)= 0 ! in [m]
            
        end if
        
    end subroutine get_ref_pwl_delay_zenith
    
end module module_get_ref_pwl_delay_zenith