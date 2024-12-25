! JB
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!
! Setting up bvp for 12-perversions solutions for the case of preferred
! long wavelength deformation  
! 
!----------------------------------------------------------------------
!----------------------------------------------------------------------

SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NDIM, IJAC, ICP(*)
  DOUBLE PRECISION, INTENT(IN) :: U(NDIM), PAR(*)
  DOUBLE PRECISION, INTENT(OUT) :: F(NDIM)
  DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDIM,NDIM),DFDP(NDIM,*)

  ! tau = nondimensional twisting strain
  ! dTau = d tau / dS
  ! d2Tau = d^2 tau / dS^2
  ! d3TauByl = (1/l) d^3 tau / dS^3
  ! MByF0 = M / F0, where M is moment 
  ! epsilon = nondimensional applied strain
  ! l = nondimensional length 
  ! Su = S/l, we use Su as the independent variable in our Auto-07p equations
  ! mExt  = the amplitude of applied external moment (used in guess step only)

  DOUBLE PRECISION theta, tau, dTau, d2Tau, d3TauByl, MByF0, epsilon, l, Su, mExt


  ! Constant and paramters
  l = 300./18. 
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
  ! is calculated from the constitutive relationship in section 3.4 in the paper

  F(5) = (10.900251965096425*d2Tau)/&
       (1 + 0.03772253492939769*tau**2 + 1.7723107713770452e-6*tau**4 + &
       2.554412451854936e-11*tau**6) + &
       MByF0/(1 + 0.03772253492939769*tau**2 + 1.7723107713770452e-6*tau**4 + &
       2.554412451854936e-11*tau**6) - &
       (50.46129268901878*(29.009999999999998 + epsilon)*tau)/&
       (1 + 0.03772253492939769*tau**2 + 1.7723107713770452e-6*tau**4 + &
       2.554412451854936e-11*tau**6) + &
       (-3.784596951676409*tau**3 + &
       dTau**2*(-0.3145465076109939*tau - 3.033988256915511e-6*tau**3) + &
       d3TauByl*dTau*(-2.5148356619598458*tau - 0.00023630810285027267*tau**3 - &
       5.108824903709871e-9*tau**5) + &
       d2Tau**2*(-0.11316760478819309*tau - 0.000010633864628262272*tau**3 - &
       2.2989712066694422e-10*tau**5) + &
       dTau**4*(3.349584502423124e-7*tau + 1.8036593059887357e-9*tau**3 + &
       3.663838204609135e-13*tau**5 + 2.1682883671936824e-17*tau**7) + &
       d2Tau*(-0.3145465076109939*tau**2 - 1.5169941284577556e-6*tau**4 + &
       dTau**2*(0.05660538894360891 - 0.00002059781235603992*tau**2 + &
       1.037335570432255e-9*tau**4 + 2.442558803072757e-13*tau**6 + &
       1.0841441835968412e-17*tau**8)))/&
       (1 + 0.03772253492939769*tau**2 + 1.7723107713770452e-6*tau**4 + &
       2.554412451854936e-11*tau**6)

  ! F(6) = M'= 0; the equilibrium equation
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

  ! We use suffix 'At0' with the variable name to imply the value of variable at Su = 0, i.e. S = 0. 
  ! For instance, thetaAt0  means theta at Su = 0
  DOUBLE PRECISION thetaAt0, tauAt0, dTauAt0, d2TauAt0, d3TauBylAt0, SuAt0

  ! We use suffix 'Atl' with the variable name to imply the value of variable at Su = 1, i.e. S = l.
  ! For instance, thetaAtl  means theta at Su = 1
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
  ! Natural boundary conditions 
  FB(2) =  -d3TauBylAt0 - 43.10163930355711*dTauAt0 + 0.003503060571691119*dTauAt0**3 - &
       0.004526704191527724*d2TauAt0*dTauAt0*tauAt0 + &
       (-0.03772253492939769*d3TauBylAt0 + 0.0038924613021862715*dTauAt0 + &
       3.156208685615939e-7*dTauAt0**3)*tauAt0**2 - &
       4.253545851304909e-7*d2TauAt0*dTauAt0*tauAt0**3 + &
       (-1.7723107713770452e-6*d3TauBylAt0 + 8.191768293671882e-7*dTauAt0 + &
       5.164204897326474e-11*dTauAt0**3)*tauAt0**4 - &
       9.195884826677771e-12*d2TauAt0*dTauAt0*tauAt0**5 + &
       (-2.554412451854936e-11*d3TauBylAt0 + 5.028672480328511e-15*dTauAt0**3)*tauAt0**6 + &
       2.168288367193683e-19*dTauAt0**3*tauAt0**8

  FB(3) =  d2TauAt0 + (-729.2609070243815 + 0.04310256978215163*dTauAt0**2)*tauAt0 + &
       0.03772253492939769*d2TauAt0*tauAt0**2 + &
       (0.12647362088247724 + 5.037042175865023e-6*dTauAt0**2)*tauAt0**3 + &
       1.7723107713770452e-6*d2TauAt0*tauAt0**4 + &
       (3.033988256915511e-6 + 1.556886285349002e-10*dTauAt0**2)*tauAt0**5 + &
       2.554412451854936e-11*d2TauAt0*tauAt0**6 + 1.0253919584499727e-15*dTauAt0**2*tauAt0**7

  ! Boundary conditions at Su = 1, i.e. Su = l
  ! Dirichlet boundary condition: theta(l) = 0
  FB(4) = thetaAtl 

  !Natural boundary conditions
  FB(5) =  -d3TauBylAtl - 43.10163930355711*dTauAtl + 0.003503060571691119*dTauAtl**3 - &
       0.004526704191527724*d2TauAtl*dTauAtl*tauAtl + &
       (-0.03772253492939769*d3TauBylAtl + 0.0038924613021862715*dTauAtl + &
       3.156208685615939e-7*dTauAtl**3)*tauAtl**2 - &
       4.253545851304909e-7*d2TauAtl*dTauAtl*tauAtl**3 + &
       (-1.7723107713770452e-6*d3TauBylAtl + 8.191768293671882e-7*dTauAtl + &
       5.164204897326474e-11*dTauAtl**3)*tauAtl**4 - &
       9.195884826677771e-12*d2TauAtl*dTauAtl*tauAtl**5 + &
       (-2.554412451854936e-11*d3TauBylAtl + 5.028672480328511e-15*dTauAtl**3)*tauAtl**6 + &
       2.168288367193683e-19*dTauAtl**3*tauAtl**8

  FB(6) = d2TauAtl + (-729.2609070243815 + 0.04310256978215163*dTauAtl**2)*tauAtl + &
       0.03772253492939769*d2TauAtl*tauAtl**2 + &
       (0.12647362088247724 + 5.037042175865023e-6*dTauAtl**2)*tauAtl**3 + &
       1.7723107713770452e-6*d2TauAtl*tauAtl**4 + &
       (3.033988256915511e-6 + 1.556886285349002e-10*dTauAtl**2)*tauAtl**5 + &
       2.554412451854936e-11*d2TauAtl*tauAtl**6 + 1.0253919584499727e-15*dTauAtl**2*tauAtl**7

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
  PAR(3)=GETP('NRM',2,U)
END SUBROUTINE PVLS
