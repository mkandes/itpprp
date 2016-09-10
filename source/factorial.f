! =====================================================================
! function FACTORIAL
!
! factorial is an function that recursively computes the value of the 
! factorial of integer n.
! --------------------------------------------------------------------- 
      INTEGER RECURSIVE FUNCTION factorial ( n ) RESULT ( nFactorial )

      IMPLICIT NONE

      INTEGER, INTENT ( IN ) :: n

      nFactorial = 0

      IF ( n == 0 ) THEN

         nFactorial = 1

      ELSE IF ( n >= 1 .AND. n <= 12 ) THEN

         nFactorial = n * factorial ( n - 1 )

      ELSE

         nFactorial = -1
         WRITE ( UNIT = 6 , FMT = * ) 'math : factorial :: &
            & ERROR - n must be an integer greater than or equal to 0, &
            & but less than or equal to 12 because KIND ( n ) may be a &
            & 32-bit integer.'
         STOP

      END IF

      RETURN
      END FUNCTION
! =====================================================================
