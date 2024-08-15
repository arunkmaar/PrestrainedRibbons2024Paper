! add variables in the boundary conditions. 

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!   Prestrained ribbon
! To do:
! add comments
! tell what are different primary variables.

! remove rod's code
! note down the replacement rules, make a list and follow it for future purposes.
!
!----------------------------------------------------------------------
!----------------------------------------------------------------------

SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP)
  !--------- ----

  ! Evaluates the algebraic equations or ODE right hand side

  ! Input arguments :
  !      NDIM   :   Dimension of the algebraic or ODE system 
  !      U      :   State variables
  !      ICP    :   Array indicating the free parameter(s)
  !      PAR    :   Equation parameters

  ! Values to be returned :
  !      F      :   Equation or ODE right hand side values

  ! Normally unused Jacobian arguments : IJAC, DFDU, DFDP (see manual)

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NDIM, IJAC, ICP(*)
  DOUBLE PRECISION, INTENT(IN) :: U(NDIM), PAR(*)
  DOUBLE PRECISION, INTENT(OUT) :: F(NDIM)
  DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDIM,NDIM),DFDP(NDIM,*)


  ! Editing here onwards for prestrained ribbons  
  DOUBLE PRECISION theta, tau, dTau, d2Tau, d3TauBylVal, MByP0, ea, lVal, Su, mExt



  !--------------
  lVal = 300./18.
  !------  
  ea = PAR(1)
  mExt = PAR(2)


  ! Primary variables
  ! centerline displacement and stretch
  theta = U(1)
  tau = U(2)
  dTau = U(3)
  d2Tau = U(4)
  d3TauBylVal = U(5)
  MByP0 = U(6)
  Su = U(7)

  ! Equilibrium equations
  F(1) = tau* lVal
  F(2) = dTau* lVal
  F(3) = d2Tau* lVal
  F(4) = d3TauBylVal * lVal * lVal

  F(5) = (10.900251965096425*d2Tau)/&
       (1 + 0.03772253492939769*tau**2 + 1.7723107713770452e-6*tau**4 + &
       2.554412451854936e-11*tau**6) + &
       MByP0/(1 + 0.03772253492939769*tau**2 + 1.7723107713770452e-6*tau**4 + &
       2.554412451854936e-11*tau**6) - &
       (50.46129268901878*(29.009999999999998 + ea)*tau)/&
       (1 + 0.03772253492939769*tau**2 + 1.7723107713770452e-6*tau**4 + &
       2.554412451854936e-11*tau**6) + &
       (-3.784596951676409*tau**3 + &
       dTau**2*(-0.3145465076109939*tau - 3.033988256915511e-6*tau**3) + &
       d3TauBylVal*dTau*(-2.5148356619598458*tau - 0.00023630810285027267*tau**3 - &
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

  F(6) = -mExt*cos(40.8407*Su)
  F(7) = 1.

END SUBROUTINE FUNC

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

SUBROUTINE STPNT(NDIM,U,PAR,T)
  !--------- -----

  ! Input arguments :
  !      NDIM   :   Dimension of the algebraic or ODE system 

  ! Values to be returned :
  !      U      :   A starting solution vector
  !      PAR    :   The corresponding equation-parameter values

  ! Note : For time- or space-dependent solutions this subroutine has
  !        the scalar input parameter T contains the varying time or space
  !        variable value.

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NDIM
  DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM),PAR(*)
  DOUBLE PRECISION, INTENT(IN) :: T


  ! Initialize the equation parameters
  PAR(1) = 0.0
  PAR(2) = 0.

  ! Initialize the solution
  U(1)  = 0.0d0
  U(2)  = 0.0d0
  U(3)  = 0.0d0
  U(4)  = 0.0d0
  U(5)  = 0.0d0
  U(6)  = 0.0d0
  U(7)  = T 


END SUBROUTINE STPNT

!----------------------------------------------------------------------
!----------------------------------------------------------------------

SUBROUTINE BCND(NDIM,PAR,ICP,NBC,U0,U1,FB,IJAC,DBC)
  !--------- ----

  ! Boundary Conditions

  ! Input arguments :
  !      NDIM   :   Dimension of the ODE system 
  !      PAR    :   Equation parameters
  !      ICP    :   Array indicating the free parameter(s)
  !      NBC    :   Number of boundary conditions
  !      U0     :   State variable values at the left boundary
  !      U1     :   State variable values at the right boundary

  ! Values to be returned :
  !      FB     :   The values of the boundary condition functions 

  ! Normally unused Jacobian arguments : IJAC, DBC (see manual)

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NDIM, ICP(*), NBC, IJAC
  DOUBLE PRECISION, INTENT(IN) :: PAR(*), U0(NDIM), U1(NDIM)
  DOUBLE PRECISION, INTENT(OUT) :: FB(NBC)
  DOUBLE PRECISION, INTENT(INOUT) :: DBC(NBC,*)

  DOUBLE PRECISION thetaC, tauC, dTauC, d2TauC, d3TauBylValC, Suc

  DOUBLE PRECISION thetaM, tauM, dTauM, d2TauM, d3TauBylValM


  ! Dirichlet boundary conditions
  ! fixed end

  thetaC   = U0(1)
  tauC  = U0(2)
  dTauC = U0(3)
  d2TauC = U0(4)
  d3TauBylValC = U0(5)
  SuC= U0(7)

  thetaM = U1(1)
  tauM  = U1(2)
  dTauM = U1(3)
  d2TauM = U1(4)
  d3TauBylValM = U1(5)


  FB(1)= thetaC

  ! fixed end
  FB(2) =  -d3TauBylValC - 43.10163930355711*dTauC + 0.003503060571691119*dTauC**3 - &
       0.004526704191527724*d2TauC*dTauC*tauC + &
       (-0.03772253492939769*d3TauBylValC + 0.0038924613021862715*dTauC + &
       3.156208685615939e-7*dTauC**3)*tauC**2 - &
       4.253545851304909e-7*d2TauC*dTauC*tauC**3 + &
       (-1.7723107713770452e-6*d3TauBylValC + 8.191768293671882e-7*dTauC + &
       5.164204897326474e-11*dTauC**3)*tauC**4 - &
       9.195884826677771e-12*d2TauC*dTauC*tauC**5 + &
       (-2.554412451854936e-11*d3TauBylValC + 5.028672480328511e-15*dTauC**3)*tauC**6 + &
       2.168288367193683e-19*dTauC**3*tauC**8

  FB(3) =  d2TauC + (-729.2609070243815 + 0.04310256978215163*dTauC**2)*tauC + &
       0.03772253492939769*d2TauC*tauC**2 + &
       (0.12647362088247724 + 5.037042175865023e-6*dTauC**2)*tauC**3 + &
       1.7723107713770452e-6*d2TauC*tauC**4 + &
       (3.033988256915511e-6 + 1.556886285349002e-10*dTauC**2)*tauC**5 + &
       2.554412451854936e-11*d2TauC*tauC**6 + 1.0253919584499727e-15*dTauC**2*tauC**7

  FB(4) = thetaM 

  FB(5) =  -d3TauBylValM - 43.10163930355711*dTauM + 0.003503060571691119*dTauM**3 - &
       0.004526704191527724*d2TauM*dTauM*tauM + &
       (-0.03772253492939769*d3TauBylValM + 0.0038924613021862715*dTauM + &
       3.156208685615939e-7*dTauM**3)*tauM**2 - &
       4.253545851304909e-7*d2TauM*dTauM*tauM**3 + &
       (-1.7723107713770452e-6*d3TauBylValM + 8.191768293671882e-7*dTauM + &
       5.164204897326474e-11*dTauM**3)*tauM**4 - &
       9.195884826677771e-12*d2TauM*dTauM*tauM**5 + &
       (-2.554412451854936e-11*d3TauBylValM + 5.028672480328511e-15*dTauM**3)*tauM**6 + &
       2.168288367193683e-19*dTauM**3*tauM**8

  FB(6) = d2TauM + (-729.2609070243815 + 0.04310256978215163*dTauM**2)*tauM + &
       0.03772253492939769*d2TauM*tauM**2 + &
       (0.12647362088247724 + 5.037042175865023e-6*dTauM**2)*tauM**3 + &
       1.7723107713770452e-6*d2TauM*tauM**4 + &
       (3.033988256915511e-6 + 1.556886285349002e-10*dTauM**2)*tauM**5 + &
       2.554412451854936e-11*d2TauM*tauM**6 + 1.0253919584499727e-15*dTauM**2*tauM**7

  FB(7) = SuC

END SUBROUTINE BCND
!----------------------------------------------------------------------
!----------------------------------------------------------------------

SUBROUTINE ICND(NDIM,PAR,ICP,NINT,U,UOLD,UDOT,UPOLD,FI,IJAC,DINT)
  !--------- ----

  ! Integral Conditions

  ! Input arguments :
  !      NDIM   :   Dimension of the ODE system 
  !      PAR    :   Equation parameters
  !      ICP    :   Array indicating the free parameter(s)
  !      NINT   :   Number of integral conditions
  !      U      :   Value of the vector function U at `time' t

  ! The following input arguments, which are normally not needed,
  ! correspond to the preceding point on the solution branch
  !      UOLD   :   The state vector at 'time' t
  !      UDOT   :   Derivative of UOLD with respect to arclength
  !      UPOLD  :   Derivative of UOLD with respect to `time'

  ! Normally unused Jacobian arguments : IJAC, DINT

  ! Values to be returned :
  !      FI     :   The value of the vector integrand 

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NDIM, ICP(*), NINT, IJAC
  DOUBLE PRECISION, INTENT(IN) :: PAR(*)
  DOUBLE PRECISION, INTENT(IN) :: U(NDIM), UOLD(NDIM), UDOT(NDIM), UPOLD(NDIM)
  DOUBLE PRECISION, INTENT(OUT) :: FI(NINT)
  DOUBLE PRECISION, INTENT(INOUT) :: DINT(NINT,*)

  !X FI(1)=

END SUBROUTINE ICND

!----------------------------------------------------------------------
!----------------------------------------------------------------------

SUBROUTINE FOPT(NDIM,U,ICP,PAR,IJAC,FS,DFDU,DFDP)
  !--------- ----
  !
  ! Defines the objective function for algebraic optimization problems
  !
  ! Supplied variables :
  !      NDIM   :   Dimension of the state equation
  !      U      :   The state vector
  !      ICP    :   Indices of the control parameters
  !      PAR    :   The vector of control parameters
  !
  ! Values to be returned :
  !      FS      :   The value of the objective function
  !
  ! Normally unused Jacobian argument : IJAC, DFDP

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NDIM, ICP(*), IJAC
  DOUBLE PRECISION, INTENT(IN) :: U(NDIM), PAR(*)
  DOUBLE PRECISION, INTENT(OUT) :: FS
  DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDIM),DFDP(*)

  !X FS=

END SUBROUTINE FOPT

!----------------------------------------------------------------------
!----------------------------------------------------------------------

SUBROUTINE PVLS(NDIM,U,PAR)
  !--------- ----

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NDIM
  DOUBLE PRECISION, INTENT(IN) :: U(NDIM)
  DOUBLE PRECISION, INTENT(INOUT) :: PAR(*)
  DOUBLE PRECISION, EXTERNAL :: GETP

  !---------------------------------------------------------------------- 
  ! NOTE : 
  ! Parameters set in this subroutine should be considered as ``solution 
  ! measures'' and be used for output purposes only.
  ! 
  ! They should never be used as `true'' continuation parameters. 
  !
  ! They may, however, be added as ``over-specified parameters'' in the 
  ! parameter list associated with the AUTO-Constant NICP, in order to 
  ! print their values on the screen and in the ``p.xxx file.
  !
  ! They may also appear in the list associated with AUTO-Constant NUZR.
  !
  !---------------------------------------------------------------------- 
  ! For algebraic problems the argument U is, as usual, the state vector.
  ! For differential equations the argument U represents the approximate 
  ! solution on the entire interval [0,1]. In this case its values must 
  ! be accessed indirectly by calls to GETP, as illustrated below.
  !---------------------------------------------------------------------- 
  !
  ! Set PAR(2) equal to the L2-norm of U(1)
  !X PAR(2)=GETP('NRM',1,U)
  !
  ! Set PAR(3) equal to the minimum of U(2)
  !X PAR(3)=GETP('MIN',2,U)
  !
  ! Set PAR(4) equal to the value of U(2) at the left boundary.
  !X PAR(4)=GETP('BV0',2,U)
  !
  ! Set PAR(5) equal to the pseudo-arclength step size used.
  !X PAR(5)=GETP('STP',1,U)

  PAR(3)=GETP('NRM',2,U)
  ! PAR(3)=GETP('NRM',3,U)

  ! Set PAR(2) and PAR(3) equal to the L2-norm of ku and Dku
  !

  !---------------------------------------------------------------------- 
  ! The first argument of GETP may be one of the following:
  !        'NRM' (L2-norm),     'MAX' (maximum),
  !        'INT' (integral),    'BV0 (left boundary value),
  !        'MIN' (minimum),     'BV1' (right boundary value).
  !
  ! Also available are
  !   'STP' (Pseudo-arclength step size used).
  !   'FLD' (`Fold function', which vanishes at folds).
  !   'BIF' (`Bifurcation function', which vanishes at singular points).
  !   'HBF' (`Hopf function'; which vanishes at Hopf points).
  !   'SPB' ( Function which vanishes at secondary periodic bifurcations).
  !---------------------------------------------------------------------- 


END SUBROUTINE PVLS

!----------------------------------------------------------------------
!----------------------------------------------------------------------
