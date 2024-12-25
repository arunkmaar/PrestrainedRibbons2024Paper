! JB
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!
! Setting up bvp for 1-perversion solutions for the case of preferred
! short wavelength deformation  
! 
!----------------------------------------------------------------------


SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NDIM, IJAC, ICP(*)
  DOUBLE PRECISION, INTENT(IN) :: U(NDIM), PAR(*)
  DOUBLE PRECISION, INTENT(OUT) :: F(NDIM)
  DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDIM,NDIM),DFDP(NDIM,*)

  !----------------------------------------------------------------------  
  ! tau = nondimensional twisting strain
  ! dTau = d tau / dS
  ! d2Tau = d^2 tau / dS^2
  ! d3TauByl = (1/l) d^3 tau / dS^3
  ! MByF0 = M / F0, where M is moment 
  ! epsilon = nondimensional applied strain
  ! l = nondimensional length 
  ! Su = S/l, we use Su as the independent variable in our Auto-07p equations
  ! mExt = the amplitude of applied external moment (non zero values used in guess step only, otherwise set zero)

  DOUBLE PRECISION theta, tau, dTau, d2Tau, d3TauByl, MByF0, epsilon, l, Su, mExt


  ! Constant and paramters
  l = 500./18.
  epsilon = PAR(1)
  mExt = PAR(2)


  ! Primary variables
  theta = U(1)
  tau = U(2)
  dTau = U(3)
  d2Tau = U(4)
  d3TauByl = U(5)
  MByF0 = U(6)
  Su = U(7)


  ! Equilibrium equations
  !F(1) =  d theta / d Su  = tau * l, and similarly others
  F(1) = tau* l
  F(2) = dTau* l
  F(3) = d2Tau* l
  F(4) = d3TauByl * l * l

  ! F(5) =  d^4 d3TauByl / d Su  = d^4 tau / d S, which
  ! is calculated from the consitutive relationship in section 3.4 in the paper

  F(5) = ((-4.3233604686883735*d2Tau)/&
       (1 + 0.02385932600124142*tau**2 + 1.0914887945929934e-6*tau**4 + &
       1.5225974561258426e-11*tau**6) + &
       MByF0/(1 + 0.02385932600124142*tau**2 + 1.0914887945929934e-6*tau**4 + &
       1.5225974561258426e-11*tau**6) - &
       (30.07824199468193*(55.135280077297786 + epsilon)*tau)/&
       (1 + 0.02385932600124142*tau**2 + 1.0914887945929934e-6*tau**4 + &
       1.5225974561258426e-11*tau**6) + &
       (-2.255868149601145*tau**3 + dTau**2*&
       (-0.18982761108352694*tau - 1.8084561083863596e-6*tau**3) + &
       d3TauByl*dTau*(-2.6510362223601573*tau - 0.00024255306546510963*tau**3 - &
       5.075324853752808e-9*tau**5) + &
       d2Tau**2*(-0.07157797800372427*tau - 6.548932767557961e-6*tau**3 - &
       1.370337710513258e-10*tau**5) + &
       dTau**4*(3.0655352092436086e-6*tau + 1.9015983192737368e-9*tau**3 + &
       2.806987695866122e-13*tau**5 + 1.2924421620474643e-17*tau**7) + &
       d2Tau*(-0.18982761108352694*tau**2 - 9.042280541931798e-7*tau**4 + &
       dTau**2*(0.06001169357741679 - 6.9667951166287056e-6*tau**2 + &
       1.4448190824359843e-9*tau**4 + 1.871325130577415e-13*tau**6 + &
       6.4622108102373216e-18*tau**8)))/&
       (1 + 0.02385932600124142*tau**2 + 1.0914887945929934e-6*tau**4 + &
       1.5225974561258426e-11*tau**6))

  ! F(6) = M'= distributed external moment; the equilibrium equation
  F(6) = 0.
  ! F(7) = d Su/ d Su = 1 
  F(7) = 1.

END SUBROUTINE FUNC

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

SUBROUTINE STPNT(NDIM,U,PAR,T)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NDIM
  DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM),PAR(*)
  DOUBLE PRECISION, INTENT(IN) :: T

  !----------------------------------------------------------------------
  ! Initialize the equation parameters
  ! epsilon = 0
  PAR(1) = 0.0
  ! external moment = 0
  PAR(2) = 0.0

  ! Initialize the solution
  U(1)  = 0.0
  U(2)  = 0.0
  U(3)  = 0.0
  U(4)  = 0.0
  U(5)  = 0.0
  U(6)  = 0.0
  U(7)  = T 


END SUBROUTINE STPNT

!----------------------------------------------------------------------
!----------------------------------------------------------------------

SUBROUTINE BCND(NDIM,PAR,ICP,NBC,U0,U1,FB,IJAC,DBC)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NDIM, ICP(*), NBC, IJAC
  DOUBLE PRECISION, INTENT(IN) :: PAR(*), U0(NDIM), U1(NDIM)
  DOUBLE PRECISION, INTENT(OUT) :: FB(NBC)
  DOUBLE PRECISION, INTENT(INOUT) :: DBC(NBC,*)

  !----------------------------------------------------------------------
  ! We use suffix 'At0' with the variable name to imply the value of variable at Su = 0, i.e. S = 0. 
  ! For instance thetaAt0  means theta at Su = 0
  DOUBLE PRECISION thetaAt0, tauAt0, dTauAt0, d2TauAt0, d3TauBylAt0, SuAt0

  ! We use suffix 'Atl' with the variable name to imply the value of variable at Su = 1, i.e. S = l.
  ! For instance thetaAtl  means theta at Su = 1
  DOUBLE PRECISION thetaAtl, tauAtl, dTauAtl, d2TauAtl, d3TauBylAtl


  thetaAt0 = U0(1)
  tauAt0 = U0(2)
  dTauAt0 = U0(3)
  d2TauAt0 = U0(4)
  d3TauBylAt0 = U0(5)
  SuAt0 = U0(7)

  thetaAtl = U1(1)
  tauAtl = U1(2)
  dTauAtl = U1(3)
  d2TauAtl = U1(4)
  d3TauBylAtl = U1(5)

  ! Boundary conditions at Su = 0
  ! Dirichlet boundary condition: theta(0) = 0
  FB(1)= thetaAt0
  ! Natural boundary conditions at Su = 0
  FB(2) = -d3TauBylAt0 - 5.677648734539267*dTauAt0 + 0.0016373261000040978*dTauAt0**3 - &
       0.0017178714720893825*d2TauAt0*dTauAt0*tauAt0 + &
       (-0.02385932600124142*d3TauBylAt0 + 0.0026100848980612*dTauAt0 + &
       1.9167302306221828e-7*dTauAt0**3)*tauAt0**2 - &
       1.5717438642139108e-7*d2TauAt0*dTauAt0*tauAt0**3 + &
       (-1.0914887945929934e-6*d3TauBylAt0 + 2.929698895585903e-7*dTauAt0 + &
       2.8629747849271183e-11*dTauAt0**3)*tauAt0**4 - &
       3.2888105052318203e-12*d2TauAt0*dTauAt0*tauAt0**5 + &
       (-1.5225974561258426e-11*d3TauBylAt0 + 2.296931032224594e-15*dTauAt0**3)*tauAt0**6 + &
       7.754652972284787e-20*dTauAt0**3*tauAt0**8

  FB(3) = d2TauAt0 + (-153.38910437962454 + 0.028713496087108503*dTauAt0**2)*tauAt0 + &
       0.02385932600124142*d2TauAt0*tauAt0**2 + &
       (0.08744332312100156 + 3.280560501121435e-6*dTauAt0**2)*tauAt0**3 + &
       1.0914887945929934e-6*d2TauAt0*tauAt0**4 + &
       (1.8084561083863596e-6 + 9.684280029977227e-11*dTauAt0**2)*tauAt0**5 + &
       1.5225974561258426e-11*d2TauAt0*tauAt0**6 + 6.112008991868504e-16*dTauAt0**2*tauAt0**7

  ! Boundary conditions at Su = 1, i.e. S = l
  ! Dirichlet boundary condition: theta(l) = 0
  FB(4) = thetaAtl 

  ! Natural boundary conditions at Su = 0
  FB(5) = -d3TauBylAtl - 5.677648734539267*dTauAtl + 0.0016373261000040978*dTauAtl**3 - &
       0.0017178714720893825*d2TauAtl*dTauAtl*tauAtl + &
       (-0.02385932600124142*d3TauBylAtl + 0.0026100848980612*dTauAtl + &
       1.9167302306221828e-7*dTauAtl**3)*tauAtl**2 - &
       1.5717438642139108e-7*d2TauAtl*dTauAtl*tauAtl**3 + &
       (-1.0914887945929934e-6*d3TauBylAtl + 2.929698895585903e-7*dTauAtl + &
       2.8629747849271183e-11*dTauAtl**3)*tauAtl**4 - &
       3.2888105052318203e-12*d2TauAtl*dTauAtl*tauAtl**5 + &
       (-1.5225974561258426e-11*d3TauBylAtl + 2.296931032224594e-15*dTauAtl**3)*tauAtl**6 + &
       7.754652972284787e-20*dTauAtl**3*tauAtl**8

  FB(6) = d2TauAtl + (-153.38910437962454 + 0.028713496087108503*dTauAtl**2)*tauAtl + &
       0.02385932600124142*d2TauAtl*tauAtl**2 + &
       (0.08744332312100156 + 3.280560501121435e-6*dTauAtl**2)*tauAtl**3 + &
       1.0914887945929934e-6*d2TauAtl*tauAtl**4 + &
       (1.8084561083863596e-6 + 9.684280029977227e-11*dTauAtl**2)*tauAtl**5 + &
       1.5225974561258426e-11*d2TauAtl*tauAtl**6 + 6.112008991868504e-16*dTauAtl**2*tauAtl**7

  FB(7) = SuAt0

END SUBROUTINE BCND
!----------------------------------------------------------------------
!----------------------------------------------------------------------

SUBROUTINE ICND(NDIM,PAR,ICP,NINT,U,UOLD,UDOT,UPOLD,FI,IJAC,DINT)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NDIM, ICP(*), NINT, IJAC
  DOUBLE PRECISION, INTENT(IN) :: PAR(*)
  DOUBLE PRECISION, INTENT(IN) :: U(NDIM), UOLD(NDIM), UDOT(NDIM), UPOLD(NDIM)
  DOUBLE PRECISION, INTENT(OUT) :: FI(NINT)
  DOUBLE PRECISION, INTENT(INOUT) :: DINT(NINT,*)
END SUBROUTINE ICND

SUBROUTINE FOPT(NDIM,U,ICP,PAR,IJAC,FS,DFDU,DFDP)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NDIM, ICP(*), IJAC
  DOUBLE PRECISION, INTENT(IN) :: U(NDIM), PAR(*)
  DOUBLE PRECISION, INTENT(OUT) :: FS
  DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDIM),DFDP(*)

END SUBROUTINE FOPT

!----------------------------------------------------------------------
!----------------------------------------------------------------------

SUBROUTINE PVLS(NDIM,U,PAR)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NDIM
  DOUBLE PRECISION, INTENT(IN) :: U(NDIM)
  DOUBLE PRECISION, INTENT(INOUT) :: PAR(*)
  DOUBLE PRECISION, EXTERNAL :: GETP

  ! We store L2 norm of twist tau in the third parameter. 
  PAR(3) = GETP('NRM', 2, U)
END SUBROUTINE PVLS
