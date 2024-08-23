C----------------------------------------------------------------------
C----------------------------------------------------------------------
C   acn : Front Heteroclinic continuation to Homotopy from Allen-Cahn-Nagumo to LLK
C----------------------------------------------------------------------
C----------------------------------------------------------------------
C
      SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP) 
C     ---------- ---- 
C 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      DIMENSION U(NDIM),PAR(*),F(NDIM),DFDU(NDIM,*),DFDP(NDIM,*)
C 
      PMU = PAR(1)
C	  A = PAR(2)
C	  PH = PAR(3)
	  V0 = PAR(4)
	  BT = PAR(5)
	  PLC = PAR(6)
	  PM = PAR(7)
	  PK = PAR(8)
	  PC = PAR(9)
C	  EP1 = 0.0172	!Notice that another EP1 below
C	  EP2 = 0.0348-0.0172	!Notice that another EP2 below
C	  SC = 1.0 ! time scale adjust parameter
C	  SE = 1.0 ! W scale parameter

      U1 = U(1)
      W1 = U(2)

      VU = V0*(TANH(BT*(U1-PLC))+PM)
	  FU = PK*U1 - VU - PC
	  VP = BT*V0/COSH(BT*(U1-PLC))/COSH(BT*(U1-PLC))
	  GU = 1.0-VP/PK+PMU

C      FW1 = PMU*W1-(U1-EP1)*(EP2+EP1-U1)*(U1-EP1-EP2*A)/EP2/EP2
	  FW2 = 6.0*FU/PK/U1/U1 + 3.0*GU*W1/U1
      F(1) = W1
      F(2) = FW2
C
      IF(IJAC.EQ.0)RETURN
C
      DFDU(1,1) = 0.0
      DFDU(1,2) = 1.0

C      F1U = A - 2.0*(1.0 + A)*(U1-EP1)/EP2
C      F1U = F1U + 3.0*(U1-EP1)*(U1-EP1)/EP2/EP2
	  V2P = 1.0/COSH(BT*(U1-PLC))/COSH(BT*(U1-PLC))
	  V2P = -2.0*BT*BT*V0*V2P*TANH(BT*(U1-PLC))
	  F2U = -12.0*FU/PK/U1/U1/U1
	  F2U = F2U+3.0*(2.0-2.0*VP/PK-GU*W1)/U1/U1
	  F2U = F2U-3.0*V2P*W1/U1/PK

      DFDU(2,1) = F2U
      DFDU(2,2) = 3.0*GU/U1
C
      IF(IJAC.EQ.1)RETURN
C
C No parameter derivatives are specified with this example
C
      RETURN
      END
C
      SUBROUTINE STPNT(NDIM,U,PAR,T)
!     ---------- ----- 

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM
      DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM),PAR(*)
      DOUBLE PRECISION, INTENT(IN) :: T

       PAR(1)= 0.0 !MU
C	   PAR(2) = 0.35	!A
C	   PAR(3) = 1.0	!phi = PH
	   PAR(4) = 0.0168	!V0 = V0
	   PAR(5) = 89.7	!beta = BT
	   PAR(6) = 0.025	!Lc = PLC
	   PAR(7) = 0.913	!M = PM
	   PAR(8) = 1.25	!k = PK
	   PAR(9) = 0.0155	!c=PC

       U(1)=0.0266751
       U(2)=0.0
	   
      END SUBROUTINE STPNT
	  
      SUBROUTINE BCND 
      END SUBROUTINE BCND

      SUBROUTINE ICND 
      END SUBROUTINE ICND

      SUBROUTINE FOPT 
      END SUBROUTINE FOPT

      SUBROUTINE PVLS
      END SUBROUTINE PVLS




