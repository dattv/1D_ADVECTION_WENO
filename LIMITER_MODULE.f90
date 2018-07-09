MODULE LIMITER_MODULE
    USE IR_Precision
    CONTAINS
!
!----------------------------------------------------------------------*
!
    SUBROUTINE MUSCLE(NCELL, U, UL, UR, LIMITE)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NCELL
    REAL(R8P), DIMENSION(0:NCELL), INTENT(IN) :: U
    REAL(R8P), DIMENSION(0:NCELL), INTENT(INOUT)  :: UL, UR
    INTEGER, INTENT(IN) :: LIMITE
    
    REAL(R8P) :: TOLLIM, OMEGA
    REAL(R8P) :: SIGN, DUPW, DLOC, DELTA, RATIO
    INTEGER :: I
    
    ! BODY
    TOLLIM = 1.E-5; OMEGA = 1._R8P/3._R8P
    
    DO I = 1, NCELL-1, 1
        DUPW = U(I) - U(I-1)
        DLOC = U(I+1) - U(I)
        
        IF(ABS(DUPW).LE.TOLLIM) DUPW=TOLLIM*SIGN(1.0D0,DUPW)
        IF(ABS(DLOC).LE.TOLLIM) DLOC=TOLLIM*SIGN(1.0D0,DLOC)
        
        DELTA = 0.5D0*(1.0D0+OMEGA)*DUPW + 0.5D0*(1.0D0-OMEGA)*DLOC
        RATIO = DUPW/DLOC 
        
!       Compute slope limiter functions. The subroutines carry
!       DELTA, multiply it by the slope limiter and  return
!       a limited DELTA to be used in the boundary extrapolation
!       step.
!       
!       Slope limiters used are:
!       
!       LIMITE = 1, Godunov's first order upwind method
!       LIMITE = 2, upwind second order method (non-monotone)
!       LIMITE = 3, upwind TVD, with SUPERBEE type limiter
!       LIMITE = 4, upwind TVD, with VAN LEER type limiter
!       LIMITE = 5, upwind TVD, with VAN ALBADA type limiter
!       LIMITE = 6, upwind TVD, with MINMOD type limiter
!       LIMITE = 7, upwind TVD, with MINMAX type limiter   
        
        IF(LIMITE.EQ.1) DELTA = 0.0
        IF(LIMITE.EQ.2) DELTA = DELTA
        IF(LIMITE.EQ.3) CALL SBSLIC(RATIO, OMEGA, DELTA)
        IF(LIMITE.EQ.4) CALL VLSLIC(RATIO, OMEGA, DELTA)
        IF(LIMITE.EQ.5) CALL VASLIC(RATIO, OMEGA, DELTA)
        IF(LIMITE.EQ.6) CALL MISLIC(RATIO, OMEGA, DELTA)
        IF(LIMITE.EQ.7) CALL MINMAX(DUPW, DLOC, DELTA)
        
        UL(I) = U(I) - 0.5D0*DELTA
        UR(I) = U(I) + 0.5D0*DELTA
    END DO        
    END SUBROUTINE MUSCLE
!
!----------------------------------------------------------------------*
!
      SUBROUTINE SBSLIC(R, OMEGA, DELTA)
!
!     Purpose: to compute a SUPERBEE type slope limiter DELTA
!
      IMPLICIT NONE
!
!     Declaration of variables
!
      REAL(R8P)  DELTA, DENOR, OMEGA, PHI, PHIR, R
!
      PHI             = 0.D0
      IF(R.GE.0.0)PHI = 2.D0*R
      IF(R.GE.0.5)PHI = 1.D0
!
      IF(R.GE.1.0)THEN
         DENOR = 1.D0 - OMEGA + (1.D0 + OMEGA)*R
         PHIR  = 2.D0/DENOR
         PHI   = MIN(PHIR, R)
         PHI   = MIN(PHI, 2.0D0)
      ENDIF
!
      DELTA = PHI*DELTA
!
      END
!
!----------------------------------------------------------------------*
!
      SUBROUTINE VLSLIC(R, OMEGA, DELTA)
!
!     Purpose: to compute a VAN LEER type slope limiter DELTA
!
      IMPLICIT NONE
!
!     Declaration of variables
!
      REAL(R8P)  DELTA, DENOR, OMEGA, PHI, PHIR, R
!
      PHI = 0.D0
!
      IF(R.GE.0.D0)THEN
         DENOR = 1.D0 - OMEGA + (1.D0 + OMEGA)*R
         PHIR  = 2.D0/DENOR
         PHI   = 2.D0*R/(1.D0 + R)
         PHI   = MIN(PHI, PHIR)
      ENDIF
!
      DELTA    = PHI*DELTA
!
      END
!
!----------------------------------------------------------------------*
!
      SUBROUTINE VASLIC(R, OMEGA, DELTA)
!
!     Purpose: to compute a VAN ALBADA type slope limiter DELTA
!
      IMPLICIT NONE
!
!     Declaration of variables
!
      REAL(R8P)  DELTA, DENOR, OMEGA, PHI, PHIR, R
!
      PHI = 0.0
!
      IF(R.GE.0.0)THEN
         DENOR = 1.D0 - OMEGA + (1.D0 + OMEGA)*R
         PHIR  = 2.D0/DENOR
         PHI   = R*(1.D0 + R)/(1.D0 + R*R)
         PHI   = MIN(PHI, PHIR)
      ENDIF
!
      DELTA    = PHI*DELTA
!
      END
!
!----------------------------------------------------------------------*
!
      SUBROUTINE MISLIC(R, OMEGA, DELTA)
!
!     Purpose: to compute a MINMOD type slope limiter DELTA
!
      IMPLICIT NONE
!
!     Declaration of variables
!
      REAL(R8P)  DELTA, DENOR, OMEGA, PHI, PHIR, R
!
      PHI             = 0.D0
      IF(R.GE.0.D0)PHI = R
!
      IF(R.GE.1.0)THEN
         DENOR = 2.D0*(1.D0 - OMEGA + (1.D0 + OMEGA)*R)
         PHIR  = 4.D0/DENOR
         PHI   = MIN(1.D0, PHIR)
      ENDIF
!
      DELTA    = PHI*DELTA
!
      END
!
!----------------------------------------------------------------------*
!
      SUBROUTINE MINMAX(DUPW, DLOC, DELTA)
!
!     Purpose: to compute a MINMAX type slope limiter DELTA.
!              This is the most diffusive of all limiters
!              for centred schemes
!
      IMPLICIT NONE
!
!     Declaration of variables
!
      REAL(R8P)  BETAL, BETAR, DELTA, DLOC, DUPW, SIGNO
!
      BETAL = 1.D0
      BETAR = 1.D0
      SIGNO = 0.5D0*(SIGN(1.D0,DUPW) + SIGN(1.D0,DLOC))
      DELTA = SIGNO*(MIN(BETAL*ABS(DUPW),BETAR*ABS(DLOC)))
!
      END
END MODULE LIMITER_MODULE    