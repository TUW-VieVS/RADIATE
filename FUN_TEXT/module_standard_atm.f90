! module_standard_atm.f90

!****************************************************************************
!
! PURPOSE:  Module to to determine temperature, pressure as a function of height
!           using a standard atmosphere model.
!
!           "module_standard_atm" is a module containing subroutines to determine temperature, pressure as a
!           function of height using a standard atmospher model. Additionally density and
!           speed of sound are optional output parameters.
!           The subroutines are based on the US Standard Atmosphere of 1976.
!
!           The only public and callable subroutine outside the module is "standard_atm", which is the main subroutine for the calculations.
!
!           All calculations are performed in metric units, but converted to english if the user chooses.
!           Note: Calculations assume constant gravitational acceleration.
!
!           Note: Although the input height is expected to be a geopotential height, the input of ellipsoidal or geometric height leads to
!                 only neglectable differences since this is just a global average model anyway.
!
! Inputs and outputs as needed by the subroutine "standard_atm":
!
! INPUT:
!         h... geopotential height in [m] or [ft]; possible input range is [0m, 84852m]
!
!         optional: unit... default (parameter not present): input and output is metric;
!                           specified: logical type: .TRUE. for english units in input and output, .FALSE. for metric units (SI-units) in input and output
! 
! 
! OUTPUT:
!         T... temperature in [K] (default) or [R] (Rankine-scale)
!         P... pressure in [hPa] (default) or [psi]
!
!         optional: Rho... density in [kg/m^3] (default) or [lbm/ft^3]
!         optional: A..... speed of sound in [m/s] (default) or [ft/s]
!
!
!----------------------------------------------------------------------------
! 
! Fortran-file created by Armin Hofmeister
!   based on Matlab scripts created by Richard Rieber (Version: 3/17/2006 )
!   with updates by Armin Hofmeister
! 
!----------------------------------------------------------------------------
! History:
! 
! 27.01.2015: create the Fortran-file based on the Matlab-file "standard_atm.m"
! 28.01.2015: corrections
! 20.04.2015: correct comments and messages
! 31.08.2015: delete unused variable declaration "T_in" in subroutine "IsoThermal"
! 28.09.2015: comments added, minor correction of value for T0 and g
!             fix error for interpolation in case of height in between last two interpolation levels
!
!****************************************************************************

module module_standard_atm

    ! Define modules to be used
    
    
    ! declaration part
    implicit none
    private ! everything that is not explicitely declared public is private and known only within this module
    
    ! Define other constants
        
    ! acceleration of gravity [m/s/s]
    double precision, parameter, private :: g = 9.80665d0
        
    ! ratio of specific heats
    double precision, parameter, private :: gamma = 1.4d0
        
    ! gas constant for air [J/kg-K]
    double precision, parameter, private :: R = 287d0
    
    ! set main subroutine to public
    public :: standard_atm
    private :: Gradient, IsoThermal
    
    
contains
    
    subroutine standard_atm( height, &
                             T, &
                             P, &
                             unit, &
                             Rho, &
                             A )
    
        ! Define modules to be used
    
    
        !----------------------------------------------------------------------------
    
        ! Variable definitions
        implicit none
    
    
        ! dummy arguments
        !----------------
    
        ! INPUT
    
        ! input variable for height
        double precision, intent(in) :: height
    
        ! optional input variable for unit specification
        logical, intent(in), optional :: unit
    
    
        ! OUTPUT
    
        ! output variable for temperature
        double precision, intent(out) :: T
    
        ! output variable for pressure
        double precision, intent(out) :: P
    
        ! optional output variable for density
        double precision, intent(out), optional :: Rho
    
        ! optional output variable for speed of sound
        double precision, intent(out), optional :: A
    
    
        ! local arguments
        !----------------
        
        ! CONSTANTS
        
        ! Define transformation constants between metric and english units
        
        ! transformation constant from [m] to [ft]
        double precision, parameter :: m2ft = 3.2808d0
        
        ! transformation constant from [K] to [R]
        double precision, parameter :: K2R = 1.8d0
        
        ! transformation constant from [pa] to [psi]
        double precision, parameter :: pa2psi = 1.450377377302092d-4
        
        ! transformation constant from [kg] to [lmb]
        double precision, parameter :: kg2lbm = 2.20462262184878d0
        
        
        ! Define constants of the standard atmosphere model US 1976
        
        ! altitudes [m]
        double precision, parameter :: Start = 0d0
        double precision, parameter :: H1 = 11000d0
        double precision, parameter :: H2 = 20000d0
        double precision, parameter :: H3 = 32000d0
        double precision, parameter :: H4 = 47000d0
        double precision, parameter :: H5 = 51000d0
        double precision, parameter :: H6 = 71000d0
        double precision, parameter :: H7 = 84852d0
        
        ! lapse rates [K/m]
        double precision, parameter :: L1 = -0.0065d0
        double precision, parameter :: L2 = 0d0
        double precision, parameter :: L3 = 0.001d0
        double precision, parameter :: L4 = 0.0028d0
        double precision, parameter :: L5 = 0d0
        double precision, parameter :: L6 = -0.0028d0
        double precision, parameter :: L7 = -0.002d0
        
        
        ! OTHER LOCAL VARIABLES
        
        ! variable for locally storing the unit specification (changeable local variable compared to "unit" with intent(in) specification)
        logical :: unit_used
        
        ! variable for storing the input height, which can be changed (by transformation) compared to the intent(in) variable "height"
        double precision :: height_used
        
        ! variable for temperature in [K]
        double precision :: T_temp
        
        ! variable for pressure in [Pa]
        double precision :: P_temp
        
        ! variable for density in [kg/m^3]
        double precision :: Rho_temp
    
        !----------------------------------------------------------------------------
    
        !============================================================================
        ! Body of the subroutine
        !============================================================================    
    
    
        ! Handling of optional input argument "unit", assign to local "unit_used" variable
        ! note: reassigning of the optional input arguments to a new variable is necessary as a dummy argument
        !       corresponding to an actual argument that is not present must not be given
        !       a value within the procedure
        
        ! check if "unit" is specified as input argument
        if (present(unit)) then
            ! set "unit_used" to input value (may be .FALSE.= metric or .TRUE.= english)
            unit_used= unit
        ! in case "unit" is not specified
        else
            ! use metric unit as default (unit= .FALSE.)
            unit_used= .FALSE.
        end if
    
    
        ! Determine if height is >0 as reqired by the standard atmosphere model
        
        if (height < 0) then
            ! report error message in case that input height is below 0 m
            write(unit= *, fmt= '(a)') 'Error: Input height to standard atmosphere subroutine is below 0 m. For this value the standard atmosphere is not defined. Program is stopped!'
            ! stop the program
            stop
        end if
        
        
        ! Check unit of input height and transform to [m] if it is [ft]
        
        ! transform input height to [m] if necessary
        ! note: As "height" is an input variable with intent(in) parameter, it can not be redefined,
        !       therefore the local variable "height_used" is used for storing the input height!
        
        ! check if input height is in [ft] as "unit_used" is set to .TRUE. by user, which means english units in input
        if (unit_used) then
            height_used = height / m2ft ! in [ft]
        else
            height_used = height
        end if
        
        
        ! set initial values
        
        ! temperature in [K]
        T_temp = 288.15d0
        
        ! pressure in [Pa]
        P_temp = 1.01325d5
        
        ! density in [kg/m^3]
        Rho_temp = 1.225d0
        
        
        ! Calculate temperature, pressure and density for desired input height using the US Standard Atmosphere of 1976
        
        if (height_used <= H1) then
            
            call Gradient( Start, height_used, T_temp, P_temp, Rho_temp, L1 )
            
        else if ( (height_used > H1) .AND. (height_used <= H2) ) then
            
            call Gradient( Start, H1, T_temp, P_temp, Rho_temp, L1 )
            call IsoThermal( H1, height_used, T_temp, P_temp, Rho_temp )
            
        else if ( (height_used > H2) .AND. (height_used <= H3) ) then
            
            call Gradient( Start, H1, T_temp, P_temp, Rho_temp, L1 )
            call IsoThermal( H1, H2, T_temp, P_temp, Rho_temp )
            call Gradient( H2, height_used, T_temp, P_temp, Rho_temp, L3 )
            
        else if ( (height_used > H3) .AND. (height_used <= H4) ) then
            
            call Gradient( Start, H1, T_temp, P_temp, Rho_temp, L1 )
            call IsoThermal( H1, H2, T_temp, P_temp, Rho_temp )
            call Gradient( H2, H3, T_temp, P_temp, Rho_temp, L3 )
            call Gradient( H3, height_used, T_temp, P_temp, Rho_temp, L4 )
            
        else if ( (height_used > H4) .AND. (height_used <= H5) ) then
            
            call Gradient( Start, H1, T_temp, P_temp, Rho_temp, L1 )
            call IsoThermal( H1, H2, T_temp, P_temp, Rho_temp )
            call Gradient( H2, H3, T_temp, P_temp, Rho_temp, L3 )
            call Gradient( H3, H4, T_temp, P_temp, Rho_temp, L4 )
            call IsoThermal( H4, height_used, T_temp, P_temp, Rho_temp )
            
        else if ( (height_used > H5) .AND. (height_used <= H6) ) then
            
            call Gradient( Start, H1, T_temp, P_temp, Rho_temp, L1 )
            call IsoThermal( H1, H2, T_temp, P_temp, Rho_temp )
            call Gradient( H2, H3, T_temp, P_temp, Rho_temp, L3 )
            call Gradient( H3, H4, T_temp, P_temp, Rho_temp, L4 )
            call IsoThermal( H4, H5, T_temp, P_temp, Rho_temp )
            call Gradient( H5, height_used, T_temp, P_temp, Rho_temp, L6 )
            
        else if ( (height_used > H6) .AND. (height_used <= H7) ) then
            
            call Gradient( Start, H1, T_temp, P_temp, Rho_temp, L1 )
            call IsoThermal( H1, H2, T_temp, P_temp, Rho_temp )
            call Gradient( H2, H3, T_temp, P_temp, Rho_temp, L3 )
            call Gradient( H3, H4, T_temp, P_temp, Rho_temp, L4 )
            call IsoThermal( H4, H5, T_temp, P_temp, Rho_temp )
            call Gradient( H5, H6, T_temp, P_temp, Rho_temp, L6 )
            call Gradient( H6, height_used, T_temp, P_temp, Rho_temp, L7 )
            
        else
            ! report error message in case that input height is above the upper limit (H7)
            write(unit= *, fmt= '(a, f10.3, a, f9.3, a)') 'Error: Input height of ', height_used, 'm exceeds upper limit of the standard atmosphere model of ', H7, ' m! Program is stopped!'
            ! stop the program
            stop
        end if
        
        
        ! Define (possible) output parameters
        
        ! check what units are desired for the output
        ! if "unit_used"= .FALSE. --> metric units (as calculated)
        if (.NOT. unit_used) then
            
            ! mandatory output parameters
            T = T_temp ! in [K]
            P = P_temp / 100 ! in [hPa]
            
            ! optional output parameters
            ! determine if output is desired
            if (present(Rho)) then
                Rho = Rho_temp ! in [kg/m^3]
            end if
            
            if (present(A)) then
                A = ( R * gamma * T_temp )**0.5 ! in [m/s]
            end if
        
        ! if "unit_used"= .TRUE. --> conversion to english units
        else
            
            ! mandatory output parameters
            T = T_temp * K2R ! in [R]
            P = P_temp * pa2psi ! in [psi]
            
            ! optional output parameters
            ! determine if output is desired
            if (present(Rho)) then
                Rho = Rho_temp * kg2lbm / (m2ft**3) ! in [lbm/ft^3]
            end if
        
            if (present(A)) then
                A = ( R * gamma * T_temp )**0.5 * m2ft ! in [ft/s]
            end if
            
        end if
        

    end subroutine standard_atm
    
    
    
    !--------------------------- auxiliary subroutines --------------------------------------------
    
    subroutine Gradient( Z0, Z1, T, P, Rho, Lapse )
        
        ! PURPOSE:  Gradient calculation (interpolation at interpolation height level)
    
        ! Define modules to be used
    
    
        !----------------------------------------------------------------------------
    
        ! Variable definitions
        implicit none
    
    
        ! dummy arguments
        !----------------
    
        ! INPUT
    
        ! input variable for lower height level
        double precision, intent(in) :: Z0
        
        ! input variable for interpolation height level
        double precision, intent(in) :: Z1
        
        
        ! input variable for lapse rate
        double precision, intent(in) :: Lapse
        
        
        ! INPUT and OUTPUT
        
        ! input/output variable for temperature
        double precision, intent(in out) :: T
        
        ! input/output variable for pressure
        double precision, intent(in out) :: P
        
        ! input/output variable for density
        double precision, intent(in out) :: Rho
    
    
        ! local arguments
        !----------------
    
        ! variable for temporal storage of input temperature
        double precision :: T_in
        
    
        !----------------------------------------------------------------------------
    
        !============================================================================
        ! Body of the subroutine
        !============================================================================
        
        ! assign input temperature to local variable (avoiding that input value is lost)
        T_in = T
        
        T = T_in + Lapse * (Z1 - Z0)
        P = P * ( T / T_in )**( -g / ( Lapse * R ) )
        !Rho = Rho * ( T / T_in)**( -( g / (Lapse * R ) ) + 1 )
        Rho = P / ( R * T )
    
    end subroutine Gradient
    
    
    subroutine IsoThermal( Z0, Z1, T, P, Rho )
        
        ! PURPOSE:  Isothermal calculation
    
        ! Define modules to be used
    
    
        !----------------------------------------------------------------------------
    
        ! Variable definitions
        implicit none
    
    
        ! dummy arguments
        !----------------
    
        ! INPUT
    
        ! input variable for lower height level
        double precision, intent(in) :: Z0
        
        ! input variable for interpolation height level
        double precision, intent(in) :: Z1
        
        
        ! INPUT and OUTPUT
        
        ! input/output variable for temperature
        double precision, intent(in out) :: T
        
        ! input/output variable for pressure
        double precision, intent(in out) :: P
        
        ! input/output variable for density
        double precision, intent(in out) :: Rho
        
        
        !----------------------------------------------------------------------------
    
        !============================================================================
        ! Body of the subroutine
        !============================================================================
        
        !note: this is not necessary to declare as input = output
        !T = T
        P = P * exp( -( g / ( R * T ) ) * ( Z1 - Z0 ) )        
        !Rho = Rho * exp( -( g / ( R * T ) ) * ( Z1 - Z0 ) )
        Rho = P / ( R * T )
    
    end subroutine IsoThermal
    
    
end module module_standard_atm