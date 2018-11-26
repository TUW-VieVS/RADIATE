! module_char_case.f90 

!****************************************************************************
!
! PURPOSE:  Module inheriting elemental subroutines for changing the case of all characters
!           in an input string (elemental procedure enabling array input).
!           Only literals a-z or A-Z are changed no special literals or numbers.
!           Subroutines for changing case rely on ASCII code, so this subroutines
!           should work on any computer not depending on the default character table.
!           
!           The used case subroutine algorithms are implemented from
!           program code provided by rosettacode.org.
!
!
!
! INPUT:
!         str...   character string (array) which will be changed in terms of case to upper or lower case (each character)
! 
! 
! OUTPUT:
!         str...   same as input, but changed in case depending on the called subroutine
!
!----------------------------------------------------------------------------
! 
! Fortran-file created by Armin Hofmeister
!   based on program code provided by rosettacode.org (string case)
! 
!----------------------------------------------------------------------------
! History:
! 
! 11.11.2014: create the Fortran-file
! 13.11.2014: comments
! 20.04.2015: correct comments
! 14.03.2016: correct comments
!
!****************************************************************************
    
    
module module_char_case
    
contains
    
    ! elemental subroutine for uppercase conversion (array input possible)
    elemental subroutine upper(str)
        
        ! declaration part
        implicit none
    
        ! dummy argument
        character(len=*), intent(in out) :: str
        
        ! local variable
        integer :: i
        
        ! loop over whole string
        do i = 1, len(str) 
            ! determine if current character is lower case
            select case (str(i:i))
                case('a':'z')
                    ! change case through getting correct ascii code of upper case character (-32)
                    str(i:i) = achar(iachar(str(i:i)) - 32)
            end select
        end do
    
    end subroutine upper


    ! elemental subroutine for lowercase conversion (array input possible)
    elemental subroutine lower(str)
        
        ! declaration part
        implicit none
        
        ! dummy argument
        character(len=*), intent(in out) :: str
        
        ! local variable
        integer :: i
    
        do i = 1, len(str)
            ! determine if current character is upper case
            select case (str(i:i))
                ! change case through getting correct ascii code of lower case character (+32)
                case('A':'Z')
                    str(i:i) = achar(iachar(str(i:i)) + 32)
            end select
        end do
    
    end subroutine lower
    
    
end module module_char_case