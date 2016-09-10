! =====================================================================
! function ALAGUERRE
!
! alaguerre is an function that recursively computes the value of the 
! associated Laguerre polynomial of degree n and order k at the point x.
! --------------------------------------------------------------------- 
      DOUBLE PRECISION RECURSIVE FUNCTION alaguerre ( n , k , x ) RESULT ( alaguerreNK ) 

      IMPLICIT NONE

      INTEGER, INTENT ( IN ) :: n
      INTEGER, INTENT ( IN ) :: k

      DOUBLE PRECISION, INTENT ( IN ) :: x

      alaguerreNK = 0.0

      IF ( k > -1 ) THEN

         IF ( n == 0 ) THEN

            alaguerreNK = 1.0

         ELSE IF ( n == 1 ) THEN

            alaguerreNK = 1.0 + REAL ( k ) - x

         ELSE IF ( n >= 2 ) THEN

            alaguerreNK = ( ( 2.0 * REAL ( n ) + REAL ( k ) - 1.0 - x ) * alaguerre ( n - 1 , k , x ) - & 
               & ( REAL ( n ) + REAL ( k ) - 1.0 ) * alaguerre ( n - 2 , k , x ) ) / REAL ( n )

         ELSE

            ! ERROR

         END IF
            
      ELSE

         ! ERROR 

      END IF

      RETURN 
      END FUNCTION
! =====================================================================

