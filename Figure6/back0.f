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
	  A = PAR(2)
	  PH = PAR(3)
C	  PLAM = PAR(4)
	  PK = PAR(8)
	  PC = PAR(9)
	  
	  V0 = 0.0168	!V0 = V0
	  BT = 89.7	!beta = BT
	  PLC = 0.025	!Lc = PLC
	  PM = 0.913	!M = PM
	  
C	  PLAM = 2.0*PK/(PMU+1.0)
	  PAR(4) = 2.0*PAR(8)/(PAR(1)+1.0)
C	  PK=PLAM*(PMU+1.0)/2.0
	  
	  EP1 = 0.0172	!Notice that another EP1 below
	  EP2 = 0.0348-0.0172	!Notice that another EP2 below
	  SC = 0.02 ! time scale adjust parameter
	  SE = 1.0 ! W scale parameter

      U1 = U(1)
      W1 = U(2)

      VU = V0*(TANH(BT*(U1-PLC))+PM)
	  FU = PK*U1 - VU - PC
	  VP = BT*V0/COSH(BT*(U1-PLC))/COSH(BT*(U1-PLC))
	  GU = 1.0-VP/PK+PMU

      FW1 = PMU*W1-(U1-EP1)*(EP2+EP1-U1)*(U1-EP1-EP2*A)/EP2/EP2
	  FW2 = 6.0*FU/PK/U1/U1 + 3.0*GU*SE*W1/U1
      F(1) = ((1.0-PH)+PH*SC*SE)*W1
      F(2) = (1.0-PH)*FW1 + PH*SC/SE*FW2
C
      IF(IJAC.EQ.0)RETURN
C
      DFDU(1,1) = 0.0
      DFDU(1,2) = (1.0-PH)+PH*SC*SE

      F1U = A - 2.0*(1.0 + A)*(U1-EP1)/EP2
      F1U = F1U + 3.0*(U1-EP1)*(U1-EP1)/EP2/EP2
	  V2P = 1.0/COSH(BT*(U1-PLC))/COSH(BT*(U1-PLC))
	  V2P = -2.0*BT*BT*V0*V2P*TANH(BT*(U1-PLC))
	  F2U = -12.0*FU/PK/U1/U1/U1
	  F2U = F2U+3.0*(2.0-2.0*VP/PK-GU*SE*W1)/U1/U1
	  F2U = F2U-3.0*V2P*SE*W1/U1/PK

      DFDU(2,1) = (1.0-PH)*F1U + PH*SC/SE*F2U
      DFDU(2,2) = PMU*(1.0-PH) + 3.0*PH*SC*GU/U1
C
      IF(IJAC.EQ.1)RETURN
C
C No parameter derivatives are specified with this example
C
      RETURN
      END
C
      SUBROUTINE STPNT(NDIM,U,PAR,T)
C     ----------------
C
C Sets parameter values for homoclinic bifurcation analysis (IPS=9).
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION U(NDIM),PAR(*)
C
C COMMON block needed if IPS=9 (homoclinic bifurcations) :
      COMMON /BLHOM/ ITWIST,ISTART,IEQUIB,NFIXED,NPSI,NUNSTAB,NSTAB
C
C----------------------------------------------------------------------
C Problem parameters (only PAR(1-9) are available to the user) :
C
        PAR(1) = SQRT(2.0)*(0.5-0.35)	!PMU= C and MU for each eq.
		PAR(2) = 0.35	!A
		PAR(3) = 0.0	!phi = PH
		PAR(8) = 1.25	!k = PK
		PAR(9) = 0.0163	!c=PC
	  EP1 = 0.0172	!Notice that another EP1 below
	  EP2 = 0.0348-0.0172	!Notice that another EP2 below
C
        PAR(11)=  100.0D0         ! truncated time interval 
C----------------------------------------------------------------------
C If IEQUIB=1 then put initial equilibrium in PAR(11+i), i=1,...,NDIM :
C
        IF (IEQUIB.NE.0) THEN
          PAR(12) = EP1
          PAR(13) = 0.0
          PAR(14) = EP2+EP1
          PAR(15) = 0.0
        ENDIF
C----------------------------------------------------------------------
C IF ISTART=2 then put analytic homoclinic orbit here with T in the
C   interval [0,1]
C 
C test example (a=0,b=1)
C
C	  U(x)=EP2/(1 + Exp[x/Sqrt[2]])+EP1
      S=(T-0.5)*PAR(11)
      U(1) = EP2/(1.0+EXP(-S/SQRT(2.0)))+EP1
      U(2) = EP2*EXP(S/SQRT(2.0))/SQRT(2.0)
	  U(2) = U(2)/(1.0+EXP(S/SQRT(2.0)))/(1.0+EXP(S/SQRT(2.0)))
C
      RETURN
      END
C
C
C
      SUBROUTINE PVLS(NDIM,U,PAR)
C     ---------- ----
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      DIMENSION U(*),PAR(*)
C Homoclinic bifurcations COMMON block needed here :
      COMMON /BLHOM/ ITWIST,ISTART,IEQUIB,NFIXED,NPSI,NUNSTAB,NSTAB
C
C If IEQUIB=0 put analytic equilibrium in PAR(11+i), i=1,...,NDIM :
C
      DO I=1,NDIM
      PAR(11+i)=0
      ENDDO
	  EP1 = 0.0172	!Notice that another EP1 below
	  EP2 = 0.0348-0.0172	!Notice that another EP2 below

          PAR(12) = EP1
          PAR(13) = 0.0
          PAR(14) = EP2+EP1
          PAR(15) = 0.0
C
      RETURN
      END
C
      SUBROUTINE BCND
      RETURN
      END
C
      SUBROUTINE ICND
      RETURN
      END
C
      SUBROUTINE FOPT
      RETURN
      END




