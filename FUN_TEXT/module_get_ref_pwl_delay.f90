! module_get_ref_pwl_delay.f90


!****************************************************************************
!
! PURPOSE:  Module containing subroutine for support of refined piece-wise linear ray-tracing (calculation of slant refractive indices and delays)
!           
!           Note: Created as module to avoid explicit inferface for the subroutine as needed because of assumed-shape arrays in subroutine.
!           
!           "get_ref_pwl_delay" is a subroutine to calculate the mean total upper and lower slant refractive indices needed for refined piece-wise linear
!           ray-tracing.
!           Furthermore the mean total, hydrostatic and wet slant refractive indices and also the total, hydrostatic and wet slant delays for a specific
!           ray path between the last and the current intersection point are calculated.
!           
!           The equations are taken from Hobiger et al. 2008, Fast and accurate ray-tracing algorithms for
!           real-time space geodetic applications using numerical weather models, equations (27 - 31)
!           
!           Note: Don't use values of "dh" instead of the difference R(level) - R(level-1) for the calculations of n_m_lower and n_m_upper as
!                 for the first difference and therefore also the specific dh /= R(level) - R(level-1).
!                 This is because the first real and intermediate level is the station level!
!				  
!
!
! INPUT and OUTPUT: 
!
!         n_m........... vector of mean slant total RI calculated with n_lower and n_upper
!         n_m_h......... vector of mean slant hydrostatic RI calculated with n_h_lower and n_h_upper
!         n_m_w......... vector of mean slant wet RI calculated with n_w_lower and n_w_upper
!         ds_total...... vector of slant total delay in [m] for each level
!         ds_h.......... vector of slant hydrostatic delay in [m] for each level
!         ds_w.......... vector of slant wet delay in [m] for each level
!
!
! INPUT:
!         n_lower....... value for slant total RI of the intersection point in the lower level
!         n_h_lower..... value for slant hydrostatic RI of the intersection point in the lower level
!         n_w_lower..... value for slant wet RI of the intersection point in the lower level
!         n_upper....... value for slant total RI of the intersection point in the upper level
!         n_h_upper..... value for slant hydrostatic RI of the intersection point in the upper level
!         n_w_upper..... value for slant wet RI of the intersection point in the upper level
!         R............. vector of geocentric height of the levels (not intermediate levels!)
!         level......... value for index of the current intermediate level = upper level in base data set
!         s............. vector of slant distance between two consecutive levels
! 
! 
! OUTPUT:
!         n_m_lower... value for mean slant total refractive index (RI) for the lower level determined by vertical
!                      exponential decay
!         n_m_upper... value for mean slant total refractive index (RI) for the upper level determined by vertical
!                      exponential decay
!
!----------------------------------------------------------------------------
! 
! Fortran-file created by Armin Hofmeister
!   based on Matlab scripts created by Armin Hofmeister
! 
!----------------------------------------------------------------------------
! History:
! 
! 19.02.2015: create the Fortran-file based on the Matlab-file "get_ref_pwl_delay.m"
! 25.02.2015: programming
! 26.02.2015: programming
! 09.04.2015: comment
! 13.05.2015: correct comments
! 10.06.2015: comments
! 07.09.2015: move zenith delay part to an own subroutine
! 09.09.2015: correct comments
!
!****************************************************************************

module module_get_ref_pwl_delay

contains
    
    subroutine get_ref_pwl_delay( n_m, n_m_h, n_m_w, & ! whole vectors
                                  ds_total, ds_h, ds_w, & ! whole vectors
                                  n_lower, n_h_lower, n_w_lower, & ! scalars
                                  n_upper, n_h_upper, n_w_upper, & ! scalars
                                  R, & ! whole vector
                                  level, & ! scalar
                                  s, & ! whole vector
                                  n_m_lower, n_m_upper ) ! scalars
    
    
    
        
        ! Define modules to be used
        use module_constants, only: epsilon_refr_ind
    
    
        !----------------------------------------------------------------------------
    
        ! Variable definitions
        implicit none
    
    
        ! dummy arguments
        !----------------
    
        ! INPUT and OUTPUT
    
        ! intent(in out) for the following variables:
        
        ! define input and output variables for storing the mean slant total RI calculated with n_lower and n_upper (total hydrostatic, wet)
        ! note: size = total number of height levels
        double precision, dimension(:), intent(in out) :: n_m, n_m_h, n_m_w ! whole vectors
        
        ! define input and output variables for storing the slant delays in between two levels
        ! note: size = total number of height levels -1 as the delays are in between two levels
        double precision, dimension(:), intent(in out) :: ds_total, ds_h, ds_w ! whole vectors
        
        
        ! INPUT
        
        ! intent(in) for the following variables:
        
        ! define input variables for storing slant RI of the intersection point in the lower level (total hydrostatic, wet)
        double precision, intent(in) :: n_lower, n_h_lower, n_w_lower ! scalars
        
        ! define input variables for storing slant RI of the intersection point in the upper level (total hydrostatic, wet)
        double precision, intent(in) :: n_upper, n_h_upper, n_w_upper ! scalars
        
        ! define input variable for storing the geocentric height of the levels (not intermediate levels!)
        ! note: size = total number of height levels
        double precision, dimension(:), intent(in) :: R ! whole vector
        
        ! define input variable for storing the index of the current intermediate level = upper level in base data set
        integer, intent(in) :: level ! scalar
        
        ! define input variable for storing the slant distance between two consecutive levels
        ! note: size = total number of height levels -1 as the slant distance is between two consecutive levels
        double precision, dimension(:), intent(in) :: s ! whole vector
        
        
        ! OUTPUT
        
        ! intent(out) for the following variables:
        
        ! define ouput variables for storing the mean slant total refractive index (RI) for the lower/upper level determined by vertical
        ! exponential decay
        double precision, intent(out) :: n_m_lower, n_m_upper ! scalars
        
        
        ! local variables
        !----------------
        
        ! define variable for storing the auxiliary parameter C for vertical exponential decay assumption
        double precision :: C ! scalar
        
        
        ! CONSTANTS
        
        ! see "module_constants" for "epsilon_refr_ind"
        
    
        !----------------------------------------------------------------------------
    
        !============================================================================
        ! Body of the subroutine
        !============================================================================
        
        
        ! Define constant "epsilon_refr_ind" for checking if refractive index is almost equal to 1
        ! see "module_constants" for "epsilon_refr_ind"

    
        ! Calculate slant n_m_lower and slant n_m_upper needed for ray tracing and the slant total delay
        
        ! check if refractivities n_lower and n_upper are not (almost) equal to 1 --> would cause NaN
        ! results in further calculations
        if ( (abs(n_lower - 1) >= epsilon_refr_ind) .AND. (abs(n_upper - 1) >= epsilon_refr_ind) ) then ! if true than normal calculation doesn't cause NaN results
            
            ! calculate mean slant total refractive indices for the lower and upper level
            
            ! calculate auxiliary parameter C for vertical exponential assumption
            ! see Hobiger et al. 2008, equation (29) on page 8
            C= log( (n_upper - 1) / (n_lower - 1) ) / ( R(level) - R(level-1) )
            
            ! calculate mean slant total refractive index for the lower level
            ! see Hobiger et al. 2008, equation (27) on page 8, no factor of 2!!!
            n_m_lower= 1 + (n_lower - 1) / (C * ( R(level) - R(level-1) )) * (sqrt( (n_upper - 1) / (n_lower - 1) ) - sqrt( (n_lower - 1) / (n_upper - 1) ))
            
            ! calculate mean slant total refractive index for the upper level
            ! see Hobiger et al. 2008, equation (28) on page 8, no factor of 2!!!
            n_m_upper= 1 + (n_lower - 1) / (C * ( R(level) - R(level-1) )) * (sqrt( ((n_upper - 1) / (n_lower - 1))**3 ) - sqrt( (n_upper - 1) / (n_lower - 1) ))
        
            ! calculate mean slant total refractive index at the next (= "level") intersection point (value will be used for
            ! determining the delay)
            ! see Hobiger et al. 2008, equation (30) on page 8
            n_m(level)= 1 + (n_lower - 1) * sqrt( (n_upper - 1) / (n_lower - 1) )
        
            ! calculate slant total delay (slant or zenith, depending on value s)
            ds_total(level-1)= ( n_m(level) - n_m(level-1) ) / log( (n_m(level) - 1) / (n_m(level-1) - 1) ) * s(level-1) ! in [m]
        
        else
            
            ! set mean slant total refractive index for the lower level to 1
            n_m_lower= 1
            
            ! set mean slant total refractive index for the upper level to 1
            n_m_upper= 1
            
            ! set mean slant total refractive index at the next (="level") intersection point to 1
            n_m(level)= 1
            
            ! set slant total delay to 0
            ds_total(level-1)= 0 ! in [m]
            
        end if
        
        
        ! Calculate the slant hydrostatic delay
        
        ! check if refractivities n_h_lower and n_h_upper are not (almost) equal to 1 --> would cause NaN
        ! results in further calculations
        if ( (abs(n_h_lower - 1) >= epsilon_refr_ind) .AND. (abs(n_h_upper - 1) >= epsilon_refr_ind) ) then ! if true than normal calculation doesn't cause NaN results
            
            ! calculate mean slant hydrostatic refractive index at the next (="level") intersection point (value will be used for
            ! determining the delay)
            ! see Hobiger et al. 2008, equation (30) on page 8
            n_m_h(level)= 1 + (n_h_lower - 1) * sqrt( (n_h_upper - 1) / (n_h_lower - 1) )
            
            ! calculate slant hydrostatic delay (slant or zenith, depending on value s)
            ds_h(level-1)= ( n_m_h(level) - n_m_h(level-1) ) / log( (n_m_h(level) - 1) / (n_m_h(level-1) - 1) ) * s(level-1) ! in [m]
            
        else
            
            ! set mean slant hydrostatic refractive index at the next (="level") intersection point to 1
            n_m_h(level)= 1
            
            ! set slant slant hydrostatic delay to 0
            ds_h(level-1)= 0 ! in [m]
            
        end if
        
        
        ! Calculate the slant wet delay
        
        ! check if refractivities n_w_lower and n_w_upper are not (almost) equal to 1 --> would cause NaN
        ! results in further calculations
        if ( (abs(n_w_lower - 1) >= epsilon_refr_ind) .AND. (abs(n_w_upper - 1) >= epsilon_refr_ind) ) then ! if true than normal calculation doesn't cause NaN results
            
            ! calculate mean slant wet refractive index at the next (="level") intersection point (value will be used for
            ! determining the delay)
            ! see Hobiger et al. 2008, equation (30) on page 8
            n_m_w(level)= 1 + (n_w_lower - 1) * sqrt( (n_w_upper - 1) / (n_w_lower - 1) )
            
            ! calculate slant wet delay (slant or zenith, depending on value s)
            ds_w(level-1)= ( n_m_w(level) - n_m_w(level-1) ) / log( (n_m_w(level) - 1) / (n_m_w(level-1) - 1) ) * s(level-1) ! in [m]
            
        else
            
            ! set mean slant wet refractive index at the next (="level") intersection point to 1
            n_m_w(level)= 1
            
            ! set slant wet delay to 0
            ds_w(level-1)= 0 ! in [m]
            
        end if
        
    end subroutine get_ref_pwl_delay
    
end module module_get_ref_pwl_delay