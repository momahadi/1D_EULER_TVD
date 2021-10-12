!...  
!...                        <><><><><><><><><><><><><>
!...  <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!...                        <><><><><><><><><><><><><>
!...  

Module Global

IMPLICIT NONE
INTEGER , PARAMETER ::       np = 8         ,    ni = 16001
INTEGER , PARAMETER ::     imax = 10001     , istep = 100
REAL(np), PARAMETER :: L_domain = 200.0_np  ,  tmax = 50.0_np , & 
                            cfl = 1.0_np    , &
                          alpha = 0.45      ,alphab = 0.4     , kapa1 = 0.3
REAL(np), PARAMETER ::    R_gas = 287.058_np, gamma = 1.4_np
REAL(np)            :: x(ni), ucol(ni,3), dx, dt
REAL(np)            :: p(ni), rho(ni), u(ni), T(ni)

END MODULE Global

!...  
!...                        <><><><><><><><><><><><><>
!...  <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!...                        <><><><><><><><><><><><><>
!...  
!...     _____________________________________________________________
!...
!...        NUMERICAL SOLUTION OF THE ONE DIMENSIONAL EULER'S
!...        EQUATIONS USING: 
!...        1- 4th-order 6-stage RK temporal scheme 
!...        2- 4th-order pentadiagonal compact Scheme
!...        3- Linear-Compact + Non-Linear-Characteristic 
!...     _____________________________________________________________

PROGRAM EULER_1D

USE Global
IMPLICIT NONE
INTEGER :: i 

!...
!...    Open Solution Output File
!...

OPEN (UNIT=8,FILE='SOLUTION.dat' ,FORM='FORMATTED',STATUS='UNKNOWN')

!...
!...  * * * * * * * * GRID * * * * * * * *
!...

dx = L_domain / (imax-1._np)

DO i=1,imax
  x(i)=- L_domain/2._np+dx*(i-1)
END DO

!...
!...  * * * * * * * * Solver * * * * * * * *
!...

CALL slvfld 

DO i=1,imax
  WRITE(8,10) x(i),ucol(i,1),ucol(i,2),ucol(i,3)
END DO

WRITE (*,*) 'FLOW FIELD SOLVED.'
10 FORMAT(f19.14,2X,f19.14,2X,f19.14,2X,f19.14)


END PROGRAM EULER_1D

!...  
!...                        <><><><><><><><><><><><><>
!...  <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!...                        <><><><><><><><><><><><><>
!...  

SUBROUTINE slvfld 
USE Global , ONLY : np, tmax, cfl, dt, istep, gamma, &
                     p , rho, u , t, ucol
IMPLICIT NONE
REAL(np) :: telaps, err 
INTEGER :: it

!...

telaps = 0.0_np

CALL inifld

!...
!...  TIME STEPPING
!...

it = 0

DO WHILE (telaps < tmax)
  ERR = 0.0_np
  it=it+1
  CALL tstep(dt,cfl)
  !dt = .01
  IF (telaps+dt > tmax) THEN
    dt=tmax-telaps
  END IF
  CALL ber46(it,ERR,istep,telaps,tmax)
  telaps=telaps+dt
  IF (MOD(it,istep) == 0) THEN
    WRITE (*,500) telaps
  END IF
END DO

WRITE (*,*) ' '
WRITE (*,*) '_______________________________'
WRITE (*,*) '-------------------------------'
WRITE (*,*) 'END'

430 CLOSE (UNIT=99)
500 FORMAT(1X,'ELAPSED TIME = ',f9.5,' Seconds')

END SUBROUTINE slvfld

!...  
!...                        <><><><><><><><><><><><><>
!...  <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!...                        <><><><><><><><><><><><><>
!...  

SUBROUTINE  inifld
USE Global , ONLY : np, tmax, cfl, dt, istep, gamma, R_gas, &
                     x, dx, p , rho, u , t, ucol, &
                     imax
IMPLICIT NONE
REAL(np) :: pi
INTEGER :: i  , k

!...
!...  INITIALIZE FLOW FIELD 
!...

pi=4.*ATAN(1.)
DO i=1,imax
  DO k =1,3
    ucol(i,k)= 0.0_np
  END DO
    p(i) = 0.0_np
  rho(i) = 0.0_np
    t(i) = 0.0_np
    u(i) = 0.0_np
  IF (x(i) <= 0.0_np) THEN
      p(i) = 2.0_np * 1.0_np
    rho(i) = 2.0_np
  ELSE
      p(i) = 1.0_np * 1.0_np
    rho(i) = 1.0_np
  END IF
  
  t(i) = p(i) / rho(i) / R_gas
  ucol(i,1)= rho(i)
  ucol(i,2)= rho(i) * u(i)
  ucol(i,3)= p(i) /(gamma-1.0_np) + rho(i) / 2.0_np * u(i)**2

  IF (ucol(i,3) < 0.0_np) THEN
    WRITE (*,*) 'ERROR -- NONPHYSICAL TERM - INIFLD'
    WRITE (*,*) 'GRID LOCATION:  (',i,')'
    WRITE (*,*) 'PHYSICAL LOCATION:  (',x(i),')'
    WRITE (*,*) 'DENSITY: ',ucol(i,1)
    WRITE (*,*) 'MOMENTUM:',ucol(i,2),'I '
    WRITE (*,*) 'ENERGY:  ',ucol(i,3)
    STOP
  END IF
END DO

END SUBROUTINE  inifld

!...  
!...                        <><><><><><><><><><><><><>
!...  <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!...                        <><><><><><><><><><><><><>
!...  

SUBROUTINE pv(ucol,p,rho,u,t)

USE Global , ONLY : np, ni, gamma, R_gas, imax
IMPLICIT NONE
INTEGER :: i
REAL(np) , INTENT(IN)  :: ucol(ni,3)
REAL(np) , INTENT(OUT) :: p(ni), rho(ni), u(ni), t(ni)

DO i=1,imax
  rho(i)= ucol(i,1)
    u(i)= ucol(i,2)/ucol(i,1)
    p(i)= (gamma-1.0_np)*(ucol(i,3) -0.5_np/ucol(i,1)*ucol(i,2)**2)
    t(i)= p(i)/ucol(i,1)/R_gas
END DO

END SUBROUTINE pv

!...  
!...                        <><><><><><><><><><><><><>
!...  <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!...                        <><><><><><><><><><><><><>
!...  

SUBROUTINE tstep (dt,cfl)

USE Global , ONLY : np, gamma, R_gas, &
                     p , rho, u , t, ucol, &
                     imax , dx !, dt, cfl
IMPLICIT NONE
INTEGER :: i
REAL(np) :: dtcfl, dt, cfl

dt = 1._np

DO i=1,imax
  dtcfl=cfl/(ABS(u(i))/dx+SQRT(gamma*p(i)/rho(i))/dx)
  dt=DMIN1(dtcfl,dt)
END DO

END SUBROUTINE tstep

!...  
!...                        <><><><><><><><><><><><><>
!...  <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!...                        <><><><><><><><><><><><><>
!...  
!... 
!... Ref: Gaitonde, D., Visbal, M.R., “High-Order Schemes for Navier–Stokes Equations: 
!... Algorithm and Implementation into FDL3DI,” Technical Report # AFRL-VA-WP-TR-1998-3060, 
!... Air Force Research Laboratory, Wright-Patterson Air Force Base,OH, 1998.
!... 

SUBROUTINE compactfilter(u,it,imax)

USE Global , ONLY : np, ni, alpha, alphab
REAL(np), INTENT(IN OUT)                 :: u(ni,3)
INTEGER, INTENT(IN OUT)                  :: it
INTEGER, INTENT(IN)                      :: imax
REAL(np) :: a(ni), b(ni), c(ni), d(ni)
REAL(np) :: a1(0:5), a2(0:5), a3(0:5), a4(0:5),  a5(0:5)
REAL(np) :: a6(0:5), a7(0:5), a8(0:5), a9(0:5), a10(0:5), a11(0:5)
INTEGER  :: i, ist, ind, k, ii

!...
!...  COEF.
!...

a1(2) =  1.0_np/64.0_np + 31.0_np*alphab/32.0_np 
a2(2) = 29.0_np/32.0_np +  3.0_np*alphab/16.0_np 
a3(2) = 15.0_np/64.0_np + 17.0_np*alphab/32.0_np 
a4(2) = -5.0_np/16.0_np +  5.0_np*alphab/ 8.0_np 
a5(2) = 15.0_np/64.0_np - 15.0_np*alphab/32.0_np 
a6(2) = -3.0_np/32.0_np +  3.0_np*alphab/16.0_np 
a7(2) =  1.0_np/64.0_np -  1.0_np*alphab/32.0_np 

a1(2) =  1.0_np/16.0_np + 7.0_np*alphab/8.0_np
a2(2) =  3.0_np/4.0_np  +        alphab/2.0_np
a3(2) =  3.0_np/8.0_np  +        alphab/4.0_np
a4(2) = -1.0_np/4.0_np  +        alphab/2.0_np
a5(2) =  1.0_np/16.0_np -        alphab/8.0_np
a6(2) =  0.0_np
a7(2) =  0.0_np

a1(3) = -1.0_np/64.0_np+         alphab/32.0_np
a2(3) =  3.0_np/32.0_np+ 13.0_np*alphab/16.0_np
a3(3) = 49.0_np/64.0_np+ 15.0_np*alphab/32.0_np
a4(3) =  5.0_np/16.0_np+  3.0_np*alphab/8.0_np
a5(3) =-15.0_np/64.0_np+ 15.0_np*alphab/32.0_np
a6(3) =  3.0_np/32.0_np-  3.0_np*alphab/16.0_np
a7(3) = -1.0_np/64.0_np+         alphab/32.0_np

a1(4) = 11.0_np/16.0_np+  5.0_np*alpha/8.0_np
a2(4) = 15.0_np/32.0_np+ 17.0_np*alpha/16.0_np
a3(4) = -3.0_np/16.0_np+  3.0_np*alpha/8.0_np
a4(4) =  1.0_np/32.0_np-         alpha/16.0_np

a1(5) =(93.0_np+ 70.0_np * alpha)/128.0_np
a2(5) =( 7.0_np+ 18.0_np * alpha)/ 16.0_np
a3(5) =(-7.0_np+ 14.0_np * alpha)/ 32.0_np
a4(5) =  1.0_np/ 16.0_np - alpha /  8.0_np
a5(5) = -1.0_np/128.0_np + alpha / 64.0_np

a1(0)=        (193.0_np+ 126.0_np*alpha)/256.0_np
a2(0)=        (105.0_np+ 302.0_np*alpha)/256.0_np
a3(0)=15.0_np*(-1.0_np +   2.0_np*alpha)/ 64.0_np
a4(0)=45.0_np*( 1.0_np -   2.0_np*alpha)/512.0_np
a5(0)= 5.0_np*(-1.0_np +   2.0_np*alpha)/256.0_np
a6(0)=        ( 1.0_np -   2.0_np*alpha)/512.0_np

!...
!...  I-DIRECTION
!...

ind=imax
ist=1

DO k=1,3
  DO i=ist+5,ind-5
    b(i)=alpha
    d(i)=1.0_np
    a(i)=alpha
    c(i)=a1(0)   *u(i,k) +a2(0)/2.*(u(i+1,k)+u(i-1,k))  &
        +a3(0)/2.*(u(i+2,k)+u(i-2,k)) +a4(0)/2.*(u(i+3,k)+u(i-3,k))  &
        +a5(0)/2.*(u(i+4,k)+u(i-4,k)) +a6(0)/2.*(u(i+5,k)+u(i-5,k))
  END DO
  DO i=2,3
    b(ist+i-1)=alphab
    d(ist+i-1)=1.0_np
    a(ist+i-1)=alphab
    c(ist+i-1)=a1(i) *u(ist  ,k) +a2(i) *u(ist+1,k)+a3(i) *u(ist+2 ,k)  &
        +a4(i) *u(ist+3,k)+a5(i) *u(ist+4 ,k)  &
        +a6(i) *u(ist+5,k)+a7(i) *u(ist+6 ,k)
  END DO
!...
  ii=ist+4-1
  b(ii)=alpha
  d(ii)=1.0_np
  a(ii)=alpha
  c(ii)=a1(4)   *u(ii,k) +a2(4)/2.*(u(ii+1,k)+u(ii-1,k))  &
      +a3(4)/2.*(u(ii+2,k)+u(ii-2,k)) +a4(4)/2.*(u(ii+3,k)+u(ii-3,k))
  
  ii=ist+5-1
  b(ii)=alpha
  d(ii)=1.0_np
  a(ii)=alpha
  c(ii)=a1(5)   *u(ii,k) +a2(5)/2.*(u(ii+1,k)+u(ii-1,k))  &
      +a3(5)/2.*(u(ii+2,k)+u(ii-2,k)) +a4(5)/2.*(u(ii+3,k)+u(ii-3,k))  &
      +a5(5)/2.*(u(ii+4,k)+u(ii-4,k))

!...

  b(ist)=0.0
  d(ist)=1.
  a(ist)=0.0
  c(ist)=u(ist,k)
  DO i=2,3
    b(ind-i+1)=alphab
    d(ind-i+1)=1.0_np
    a(ind-i+1)=alphab
    c(ind-i+1)=a1(i) *u(ind  ,k) +a2(i) *u(ind-1,k)+a3(i) *u(ind-2 ,k)  &
        +a4(i) *u(ind-3,k)+a5(i) *u(ind-4 ,k)  &
        +a6(i) *u(ind-5,k)+a7(i) *u(ind-6 ,k)
  END DO

!...

  ii=ind-4+1
  b(ii)=alpha
  d(ii)=1.0_np
  a(ii)=alpha
  c(ii)=a1(4)   *u(ii,k) +a2(4)/2.*(u(ii+1,k)+u(ii-1,k))  &
      +a3(4)/2.*(u(ii+2,k)+u(ii-2,k)) +a4(4)/2.*(u(ii+3,k)+u(ii-3,k))
  
  ii=ind-5+1
  b(ii)=alpha
  d(ii)=1.0_np
  a(ii)=alpha
  c(ii)=a1(5)   *u(ii,k) +a2(5)/2.*(u(ii+1,k)+u(ii-1,k))  &
      +a3(5)/2.*(u(ii+2,k)+u(ii-2,k)) +a4(5)/2.*(u(ii+3,k)+u(ii-3,k))  &
      +a5(5)/2.*(u(ii+4,k)+u(ii-4,k))

!...

  b(ind)=0.0
  d(ind)=1.0_np
  a(ind)=0.0
  c(ind)=u(ind,k)
!...
  CALL trid(ist,ind,b,d,a,c)
  DO i=ist,ind
    u(i,k)=c(i)
  END DO
END DO
!...

END SUBROUTINE compactfilter

!...  
!...                        <><><><><><><><><><><><><>
!...  <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!...                        <><><><><><><><><><><><><>
!...  

SUBROUTINE trid(il,iu,bb,dd,aa,cc)
USE Global , ONLY : np , ni
IMPLICIT NONE
INTEGER, INTENT(IN)       :: il
INTEGER, INTENT(IN)       :: iu
REAL(np), INTENT(IN)      :: bb(ni)
REAL(np), INTENT(IN OUT)      :: dd(ni)
REAL(np), INTENT(IN)      :: aa(ni)
REAL(np), INTENT(IN OUT)  :: cc(ni)
REAL(np) :: r
INTEGER :: i, j, lp

!...
!...  SUBROUTINE TRID SOLVE TRIDIAGONAL SYSTEM BY ELIMENATION
!...  IL = LOWER VALUE OF FIRST EQUATION
!...  IU = LOWER VALUE OF LAST  EQUATION
!...  BB = SUB DIAGONAL COEFF.
!...  DD =     DIAGONAL COEFF.
!...  AA = SUP DIAGONAL COEFF.
!...  CC = ELEMENT OF RIGHT-HAND SIDE
!...
!...  ESTABLISH UPPER TRANGULAR MATRIX
!...

 lp = il+1
DO  i=lp,iu
  r = bb(i)/dd(i-1)
  dd(i) = dd(i)-r*aa(i-1)
  cc(i) = cc(i)-r*cc(i-1)
END DO

!...
!...  BACK SUBSTITUTION
!...

 cc(iu) = cc(iu)/dd(iu)
DO  i = lp,iu
  j = iu-i+il
  cc(j) = (cc(j)-aa(j)*cc(j+1))/dd(j)
END DO

!...
!...  SOLUTION STORED IN CC
!...

END SUBROUTINE trid

!...  
!...                        <><><><><><><><><><><><><>
!...  <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!...                        <><><><><><><><><><><><><>
!...  

SUBROUTINE comp_invflux(ucol,e)
USE Global , ONLY : np, ni, imax, &
                    p , rho, u , t
IMPLICIT NONE
REAL(np), INTENT(IN)                         :: ucol(ni,3)
REAL(np), INTENT(OUT)                        :: e(ni,3)
INTEGER :: i

CALL pv(ucol,p,rho,u,t)

DO i=1,imax
  e(i,1)=ucol(i,2)
  e(i,2)=ucol(i,2)*u(i)+p(i)
  e(i,3)=(ucol(i,3)+p(i))*u(i)
END DO

END SUBROUTINE comp_invflux

!...  
!...                        <><><><><><><><><><><><><>
!...  <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!...                        <><><><><><><><><><><><><>
!...  
!... 
!... Yee , H. C., Sandham, N. D., Djomehri, M. J., “Low Dissipative High Order 
!... Shock-Capturing Methods Using Characteristic-Based Filters,” 
!... Journal of Computational Physics, Vol. 150, pp. 199-238, 1999.
!... 

SUBROUTINE cfilter(ucol)
USE Global, ONLY : np, ni, imax, R_gas , gamma, &
                   dx, dt, p, rho , u , t, kapa1
IMPLICIT NONE
REAL(np) ::  ucol(ni,3)
REAL(np) ::  alpha(ni,3),g(ni,3),phi(3)
REAL(np) ::  kapa(3)
REAL(np) ::  s1(ni,3)
REAL(np) ::  h(ni)
REAL(np) ::  rhob(ni),ub(ni),hb(ni),ab(ni)
REAL(np) ::  thetat(ni,3)
REAL(np) ::  thetai(ni,3)
REAL(np) ::  lambda(ni,3)
REAL(np) ::  dtdxi, delta, epsilon , delta2, beta, &
             rr, drho, dm, de, omega1, omega2, gama
INTEGER :: i, ip1, l, im1

dtdxi = dt/dx
CALL pv(ucol,p,rho,u,ab,h)

  delta = 1.0_np/16.0_np
epsilon = 1.0_np/10000000.0_np
 delta2 = 1.0_np/10000000.0_np
   beta = gamma-1.0_np
kapa(1) = kapa1 
kapa(2) = kapa1
kapa(3) = kapa1

DO i=1,imax
  h(i)=gamma*p(i)/rho(i)/beta+0.5_np*u(i)**2
END DO

!...
!...  I - DIRECTION
!...
!...  Roe's Averages , I=I+1/2
!...
!...  ALPHA=R**(-1)*DELTA(UCOL) , I=I+1/2
!...  R**-1 is The Inverse of The Eigenvectors Matrix
!...

DO i=1,imax-1
  ip1=i+1
  rr= DSQRT(rho(ip1)/rho(i))
  rhob(i)= rho(i)*rr
  ub(i)=(rr*u(ip1)+u(i))/(1.0_np+rr)
  hb(i)=(rr*h(ip1)+h(i))/(1.0_np+rr)
  ab(i)= DSQRT(beta*(hb(i)-0.5_np*ub(i)**2))
  lambda(i,1)=ub(i)
  lambda(i,2)=ub(i)+ab(i)
  lambda(i,3)=ub(i)-ab(i)
!...
  drho=ucol(ip1,1)-ucol(i,1)
  dm=ucol(ip1,2)-ucol(i,2)
  de=ucol(ip1,3)-ucol(i,3)
  alpha(i,1)=(1.0_np-0.5_np*beta*ub(i)**2/ab(i)**2)*drho +beta*ub(i)/ab(i)**2*dm  &
      -beta/ab(i)**2*de
  alpha(i,2)=1./rhob(i)/ab(i)*( ( 0.5_np*beta*ub(i)**2-ub(i)*ab(i))*drho  &
      +(ab(i)-beta*ub(i))*dm +beta*de)
  alpha(i,3)=1./rhob(i)/ab(i)*( (-0.5_np*beta*ub(i)**2-ub(i)*ab(i))*drho  &
      +(ab(i)+beta*ub(i))*dm -beta*de)
END DO

!...
!... I=I
!... Flux Limiters G , Theta(TELDA)
!...

DO i=2,imax-1
  im1=i-1
  DO l=1,3
    thetat(i,l)=  ABS(ABS(alpha(i,l))-ABS(alpha(im1,l)))/  &
        MAX(ABS(alpha(i,l))+ABS(alpha(im1,l)),epsilon)
    omega1=alpha(i,l)
    omega2=alpha(im1,l)
    IF((ABS(omega1) < ABS(omega2)).AND. (omega1*omega2 > 0.0_np)) THEN
      g(i,l)=omega1 !* sign(1.,OMEGA1)
    ELSE IF((ABS(omega1) > ABS(omega2)).AND.  &
          (omega1*omega2 > 0.0_np)) THEN
      g(i,l)=omega2 !* sign(1.,OMEGA2)
    ELSE
      g(i,l)=0.0
    END IF
  END DO
END DO

DO l=1,3
  g(1,l)=g(2,l)
  g(imax,l)=g(imax-1,l)
  thetat(1,l)=thetat(2,l)
  thetat(imax,l)=thetat(imax-1,l)
END DO

!...
!... I=I+1/2
!... Theta , GAMA , Phi , Filter Numerical Flux S1
!...

DO i=1,imax-1
  ip1=i+1
  DO l=1,3
    thetai(i,l)= DMAX1(thetat(i,l),thetat(ip1,l))
    gama=0.5*DSQRT(lambda(i,l)**2+delta) *(g(ip1,l)-g(i,l))  &
        *alpha(i,l)/(alpha(i,l)**2+epsilon)
    phi(l)=kapa(l)*thetai(i,l) *(0.5*DSQRT(lambda(i,l)**2+delta)  &
        *(g(ip1,l)+g(i,l))+ DSQRT((lambda(i,l)+gama)**2+delta)  &
        *alpha(i,l))
  END DO
  s1(i,1)=0.5*(phi(1) +0.5*rhob(i)/ab(i)*(phi(2)-phi(3)))
  s1(i,2)=0.5*(ub(i)*phi(1) +0.5*rhob(i)/ab(i)*(  &
      phi(2)*(ub(i)+ab(i)) -phi(3)*(ub(i)-ab(i))))
  s1(i,3)=0.5*(0.5*ub(i)**2*phi(1) +0.5*rhob(i)/ab(i)*(  &
      phi(2)*(hb(i)+ab(i)*ub(i)) -phi(3)*(hb(i)-ab(i)*ub(i))))
END DO

DO i=2,imax-1
  DO l=1,3
    ucol(i,l)=ucol(i,l) +dtdxi*(s1(i,l)-s1(i-1,l))
  END DO
END DO

END SUBROUTINE cfilter

!...  
!...                        <><><><><><><><><><><><><>
!...  <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!...                        <><><><><><><><><><><><><>
!...  
!... 
!... Ref: Berland, J., Bogey, C., and Bailly, C., “Low-dissipation and low-dispersion 
!... fourth-order Runge–Kutta algorithm,” Computers & Fluids, short communication, 2005.
!... 

SUBROUTINE ber46(it,ERR,istep,telaps,tmax)

USE Global, ONLY : np, ni, imax, R_gas , gamma, &
                   ucol, x, dx, dt, p, rho , u , T
IMPLICIT NONE
INTEGER, INTENT(IN OUT)             :: it
REAL(np), INTENT(OUT)               :: ERR
INTEGER, INTENT(IN)                 :: istep
REAL(np), INTENT(IN OUT)            :: telaps
REAL(np), INTENT(IN)                :: tmax
INTEGER :: s , i , k, ist, ind, im
REAL(np) :: ucoln(ni,3), ws(ni,3), us(ni,3)
REAL(np) :: e(ni,3), dedxi(ni,3), alpha(6), beta(6)
REAL(np) :: dtdxi , ERRP, PO, RHOTEMP, UTEMP, PTEMP, TTEMP


dtdxi = dt 
alpha(1) =  0.0_np
alpha(2) = -0.737101392796_np
alpha(3) = -1.634740794341_np
alpha(4) = -0.744739003780_np
alpha(5) = -1.469897351522_np
alpha(6) = -2.813971388035_np

beta(1) = 0.032918605146_np
beta(2) = 0.823256998199_np
beta(3) = 0.381530948900_np
beta(4) = 0.200092213184_np
beta(5) = 1.718581042715_np
beta(6) = 0.27_np

DO k=1,3
  ucoln(1   ,k)=ucol(1   ,k)
  ucoln(imax,k)=ucol(imax,k)
END DO

!...
!...  Initialization step
!...

DO i=1, imax
  DO k=1,3
    us(i,k)=ucol(i,k)
  END DO
END DO

!...
!...  Runge-Kutta Stages
!...

DO s=1, 6
  
  CALL pv(us,p,rho,u,t)
  CALL comp_invflux(us,e)
  
  ist=1
  ind=imax
  CALL cp4iv(e,dedxi,ist,ind,dx)
  
!...  <><><><><><><><><><>
!...  BOUNDARY CONDITIONS
!...  <><><><><><><><><><>

!     CALL C_BOUNDARY_COND(DEDXI,IT)
  
  DO i=2,imax-1
    DO k=1,3
      ws(i,k)= alpha(s)*ws(i,k) - dtdxi*dedxi(i,k)
      us(i,k)= us(i,k) + beta(s)*ws(i,k)
    END DO
  END DO
  
END DO

!...
!...  Update step
!...

DO i=2,imax-1
  DO k=1,3
    ucoln(i,k)=us(i,k)
  END DO
END DO

!...
!...  FILTERING THE SOLUTION
!...

CALL cfilter(ucoln,dt,it,dx,telaps,tmax)
CALL compactfilter(ucoln,it,imax)

!...
!...  CONVERGENCE HISTORY AND UPDATING FLOW VARIABLES
!...

ERR=0._np
errp=0._np

DO  i=1,imax
  po=(gamma-1._np)*( ucol(i,3)-0.5_np/ucol(i,1) *(ucol(i,2)**2))
  rhotemp=ucoln(i,1)
  utemp=ucoln(i,2)/ucoln(i,1)
  ptemp=(gamma-1._np)* (ucoln(i,3)-0.5_np*rhotemp*(utemp**2))
  ttemp=ptemp/rhotemp/R_gas
  errp=DMAX1(errp,ABS(po-ptemp))
  IF (errp > ERR) THEN
    im=i
  END IF
  ERR=DMAX1(ERR,errp)
  rho(i)=rhotemp
  u(i)=utemp
  p(i)=ptemp
  t(i)=ttemp
  ucol(i,1)=ucoln(i,1)
  ucol(i,2)=ucoln(i,2)
  ucol(i,3)=ucoln(i,3)
END DO

DO i=1,imax
  IF (p(i) <= 1.e-6 .OR. t(i) <= 1.e-6 .OR. rho(i) <= 1.e-6 .OR. ucol(i,3) <= 1.e-6) THEN
    WRITE (*,*) 'ERROR -- NONPHYSICAL TERM - PRECOR:'
    WRITE (*,*) 'Pressure =', p(i)
    WRITE (*,*) 'Temp. =', t(i)
    WRITE (*,*) 'DENSITY: ',ucol(i,1)
    WRITE (*,*) 'MOMENTUM:',ucol(i,2),'I '
    WRITE (*,*) 'ENERGY:  ',ucol(i,3)
    WRITE (*,*) 'GRID LOCATION:  (',i,')'
    WRITE (*,*) 'PHYSICAL LOCATION:  (',x(i),')'
    STOP
  END IF
END DO

IF (MOD(it,istep) == 0) THEN
  WRITE (*,*) '_________________________________________________'
  WRITE (*,*) ' '
  WRITE (*,500) it
  WRITE (*,600) ERR
  WRITE (*,20)
  WRITE (*,30) im, x(im)
END IF

20 FORMAT(1X,'Location   : ','__I__',2X,'___X(I)_____')
30 FORMAT(14X,i5,2X,f12.8)
500 FORMAT(1X,'FLOW Time-Step # : ',i6)
600 FORMAT(1X,'Max Variation  : ',f9.7)

END SUBROUTINE ber46

!...  
!...                        <><><><><><><><><><><><><>
!...  <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!...                        <><><><><><><><><><><><><>
!...  
!... 
!... Ref-1: Kim, J. W., and Lee, D. J., “Optimized Compact Finite Difference Schemes 
!... with Maximum Resolution,” AIAA Journal. Vol.34, pp. 887-893, 1996.
!... 
!... Ref-2: Kim, J. W., and Lee, D. J., “Implementation of Boundary Conditions for
!... the Optimized High-Order Compact Schemes,” 2nd AIAA/CEAS Aeroacoustics Conference, 
!... State College, PA, May 6-8, 1996.
!... 

SUBROUTINE cp4iv(u,du,ist,ind,dx)
USE Global, ONLY : np, ni, imax, R_gas , gamma
IMPLICIT NONE
REAL(np), INTENT(IN)                    :: u(ni,3)
REAL(np), INTENT(OUT)                   :: du(ni,3)
INTEGER, INTENT(IN)                     :: ist
INTEGER, INTENT(IN)                     :: ind
REAL(np), INTENT(IN)                    :: dx
REAL(np) :: a_(ni,3),b_(ni,3),c_(ni,3)
REAL(np) :: d_(ni,3),e_(ni,3),f_(ni,3),x(ni,3)
REAL(np) :: a    ,b    ,c    ,alpha    ,beta 
REAL(np) :: a_01 ,b_02 ,c_03 ,alpha_01 ,beta_02 
REAL(np) :: a_10 ,b_12 ,c_13 ,d_14     ,alpha_10 ,alpha_12 ,beta_13
REAL(np) :: a_20 ,b_21 ,c_23 ,d_24     ,e_25     ,beta_20  ,alpha_21 ,alpha_23 ,beta_24
INTEGER :: i, k

!...
!...  COEF.
!...

      a = 0.6511278808920836_np
      b = 0.2487500014377899_np
      c = 0.006144796612699781_np
  alpha = 0.5775233202590945_np
   beta = 0.08953895334666784_np

    a_01 = -3.061503488555582_np
    b_02 =  5.917946021057852_np
    c_03 =  0.4176795271056629_np
alpha_01 =  5.870156099940824_np
 beta_02 =  3.157271034936285_np

    a_10 = -0.5401943305881343_np
    b_12 =  0.8952361063034303_np
    c_13 =  0.2553815577627246_np
    d_14 =  0.007549029394582539_np
alpha_10 =  0.1663921564068434_np
alpha_12 =  0.7162501763222718_np
 beta_13 =  0.08619830787164529_np

    a_20 = -0.1327404414078232_np
    b_21 = -0.6819452549637237_np
    c_23 =  0.7109139355526556_np
    d_24 =  0.2459462758541114_np
    e_25 =  0.003965415751510620_np
 beta_20 =  0.03447751898726934_np
alpha_21 =  0.4406854601950040_np
alpha_23 =  0.6055509079866320_np
 beta_24 =  0.08141498512587530_np

!...
!...
!...

DO k=1,3
  e_(ist,k)= 0.0_np
  b_(ist,k)= 0.0_np
  d_(ist,k)= 1.0_np
  a_(ist,k)= alpha_01
  f_(ist,k)= beta_02
  c_(ist,k)= a_01 / dx * ( u(ist+1,k) - u(ist,k)) +  &
      b_02 / dx * ( u(ist+2,k) - u(ist,k)) + c_03 / dx * ( u(ist+3,k) - u(ist,k))
  
  e_(ist+1,k)= 0.0_np
  b_(ist+1,k)= alpha_10
  d_(ist+1,k)= 1.0_np
  a_(ist+1,k)= alpha_12
  f_(ist+1,k)= beta_13
  c_(ist+1,k)= a_10 / dx * ( u(ist  ,k) - u(ist+1,k)) +  &
      b_12 / dx * ( u(ist+2,k) - u(ist+1,k)) +  &
      c_13 / dx * ( u(ist+3,k) - u(ist+1,k)) +  &
      d_14 / dx * ( u(ist+4,k) - u(ist+1,k))
  
  e_(ist+2,k)= beta_20
  b_(ist+2,k)= alpha_21
  d_(ist+2,k)= 1.0_np
  a_(ist+2,k)= alpha_23
  f_(ist+2,k)= beta_24
  c_(ist+2,k)= a_20 / dx * ( u(ist  ,k) - u(ist+2,k)) +  &
      b_21 / dx * ( u(ist+1,k) - u(ist+2,k)) +  &
      c_23 / dx * ( u(ist+3,k) - u(ist+2,k)) +  &
      d_24 / dx * ( u(ist+4,k) - u(ist+2,k)) +  &
      e_25 / dx * ( u(ist+5,k) - u(ist+2,k))
  
!...
  
  DO i=ist+3,ind-3
    
    e_(i,k)= beta
    b_(i,k)= alpha
    d_(i,k)= 1.0_np
    a_(i,k)= alpha
    f_(i,k)= beta
    c_(i,k)= a / dx * (u(i+1,k)-u(i-1,k)) + b / dx * (u(i+2,k)-u(i-2,k)) +  &
        c / dx * (u(i+3,k)-u(i-3,k))
    
  END DO

!...
  
  e_(ind,k)= beta_02
  b_(ind,k)= alpha_01
  d_(ind,k)= 1.0_np
  a_(ind,k)= 0.0_np
  f_(ind,k)= 0.0_np
  c_(ind,k)= - a_01 / dx * ( u(ind-1,k) - u(ind,k))  &
      - b_02 / dx * ( u(ind-2,k) - u(ind,k))  &
      - c_03 / dx * ( u(ind-3,k) - u(ind,k))
  
  e_(ind-1,k)= beta_13
  b_(ind-1,k)= alpha_12
  d_(ind-1,k)= 1.0_np
  a_(ind-1,k)= alpha_10
  f_(ind-1,k)= 0.0_np
  c_(ind-1,k)= - a_10 / dx * ( u(ind  ,k) - u(ind-1,k))  &
      - b_12 / dx * ( u(ind-2,k) - u(ind-1,k))  &
      - c_13 / dx * ( u(ind-3,k) - u(ind-1,k))  &
      - d_14 / dx * ( u(ind-4,k) - u(ind-1,k))
  
  e_(ind-2,k)= beta_24
  b_(ind-2,k)= alpha_23
  d_(ind-2,k)= 1.0_np
  a_(ind-2,k)= alpha_21
  f_(ind-2,k)= beta_20
  c_(ind-2,k)= - a_20 / dx * ( u(ind  ,k) - u(ind-2,k))  &
      - b_21 / dx * ( u(ind-1,k) - u(ind-2,k))  &
      - c_23 / dx * ( u(ind-3,k) - u(ind-2,k))  &
      - d_24 / dx * ( u(ind-4,k) - u(ind-2,k))  &
      - e_25 / dx * ( u(ind-5,k) - u(ind-2,k))
  
END DO

!...

CALL pentadv(ist,ind,e_,b_,d_,a_,f_,c_,x)

DO k=1,3
  DO i=ist,ind
    du(i,k) = x(i,k)
  END DO
END DO

!...

END SUBROUTINE cp4iv

!...  
!...                        <><><><><><><><><><><><><>
!...  <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!...                        <><><><><><><><><><><><><>
!...  

SUBROUTINE pentadv(il,iu,e,b,d,a,f,c,x)

USE Global, ONLY : np, ni, imax, R_gas , gamma
IMPLICIT NONE
INTEGER, INTENT(IN)                      :: il
INTEGER, INTENT(IN)                      :: iu
REAL(np), INTENT(IN)                     :: e(ni,3)
REAL(np), INTENT(IN OUT)                 :: b(ni,3)
REAL(np), INTENT(IN OUT)                 :: d(ni,3)
REAL(np), INTENT(IN OUT)                 :: a(ni,3)
REAL(np), INTENT(IN)                     :: f(ni,3)
REAL(np), INTENT(OUT)                    :: c(ni,3)
REAL(np), INTENT(OUT)                    :: x(ni,3)
INTEGER :: i, k
REAL(np) :: r

!...
!...  SUBROUTINE PENTD SOLVE PENTADIAGONAL SYSTEM BY ELIMENATION
!...  IL = LOWER VALUE OF FIRST EQUATION
!...  IU = LOWER VALUE OF LAST  EQUATION
!...   E = SUB_SUB DIAGONAL COEFF.
!...   B = SUB  DIAGONAL COEFF.
!...   D =      DIAGONAL COEFF.
!...   A = SUP  DIAGONAL COEFF.
!...   F = SUP_SUP DIAGONAL COEFF.
!...   C = ELEMENT OF RIGHT-HAND SIDE
!...   X = SOLUTION VECTOR
!...
!...  ESTABLISH UPPER TRANGULAR MATRIX
!...

DO k=1,3
  DO i = il+1,iu-1
    r = b(i-1,k)/d(i-1,k)
    d(i,k) = d(i,k) - r * a(i-1,k)
    a(i,k) = a(i,k) - r * f(i-1,k)
    c(i,k) = c(i,k) - r * c(i-1,k)
    r = e(i-1,k)/d(i-1,k)
    b(i,k)   = b(i,k)   - r * a(i-1,k)
    d(i+1,k) = d(i+1,k) - r * f(i-1,k)
    c(i+1,k) = c(i+1,k) - r * c(i-1,k)
  END DO

!...
!...  BACK SUBSTITUTION
!...

  i=iu
  r = b(i-1,k) / d(i-1,k)
  d(i,k) =   d(i,k) - r * a(i-1,k)
  x(i,k) = ( c(i,k) - r * c(i-1,k) ) / d(i,k)
  x(i-1,k) = ( c(i-1,k) - a(i-1,k) * x(i,k) ) / d(i-1,k)
  DO i = iu-2,il,-1
    x(i,k) = ( c(i,k) - f(i,k) * x(i+2,k) - a(i,k) * x(i+1,k) ) / d(i,k)
  END DO
END DO

!...
!...  SOLUTION STORED IN X
!...

END SUBROUTINE pentadv

!...  
!...                      <><><><><><><><><><><><><>
!...  <><><><><><><><><><>  Charac. Boundry Cond's  <><><><><><><><><>
!...                      <><><><><><><><><><><><><>
!...  
!... 
!... Ref-1: non-reflecting boundary conditions (pages from my master thesis)
!... 
!... Ref-2: Poinsot, T. J., and Lele, S. K., “Boundary Conditions for Direct Simulations
!... of Compressible Viscous Flow,” Journal of Computational Physics, Vol. 101, 
!... No. 1, 1992, pp. 104–129.
!... 

SUBROUTINE c_boundary_cond(dedxi,it)

USE Global , ONLY : np, ni, imax, gamma, R_gas, x, dx, &
                     p , rho, u , t
IMPLICIT NONE
REAL(np), INTENT(OUT)            :: dedxi(ni,3)
INTEGER, INTENT(IN)              :: it
INTEGER :: i
REAL(np) l1 , l2 , l3 , d_1 , d_2 , d_3 , kout
REAL(np) uc , sound_speed , m_max , pinf
REAL(np) u_x , p_x , rho_x


i=1  !... Solid Wall
sound_speed = SQRT( gamma * p(i) / rho(i) )
uc = u(i)
u_x = 1./60.*( - 110. * u(i  ) / dx  &
    + 180. * u(i+1) / dx -  90. * u(i+2) / dx  &
    +  20. * u(i+3) / dx )
p_x = 1./60.*( - 110. * p(i  ) / dx  &
    + 180. * p(i+1) / dx -  90. * p(i+2) / dx  &
    +  20. * p(i+3) / dx )
rho_x = 1./60.*( - 110. * rho(i  ) / dx  &
    + 180. * rho(i+1) / dx -  90. * rho(i+2) / dx  &
    +  20. * rho(i+3) / dx )

!     L1 = ( Uc - sound_speed ) * ( p_x - rho(i,j) * sound_speed * u_x)
!     L2 = Uc * ( p_x - sound_speed**2 * rho_x )
!     L3 = ( Uc + sound_speed ) * ( p_x + rho(i,j) * sound_speed * u_x)

l1 = ( uc - sound_speed ) * ( p_x - rho(i) * sound_speed * u_x)
l2 = uc * ( p_x - sound_speed**2 * rho_x )
l3 = 0.0 ! ( Uc + sound_speed ) * ( p_x + rho(i) * sound_speed * u_x)

d_1 = 0.5 / sound_speed**2 * ( l1 - 2. * l2 + l3 )
d_2 = 0.5 / sound_speed / rho(i) * ( l3 - l1 )
d_3 = 0.5 * ( l1 + l3 )

dedxi(i,1) = d_1
dedxi(i,2) = u(i) * d_1 + rho(i) * d_2
dedxi(i,3) = 0.5 * u(i)**2   * d_1 + rho(i) * u(i)   * d_2 +  &
    1. / (gamma-1.) * d_3

!...  >>>>>>>>>>>>>>>>>>>>>>>>>>>>

i=imax !... Discharge Exit
sound_speed = SQRT( gamma * p(i) / rho(i) )
uc =  u(i)
u_x = -1./60.*( - 110. * u(i  ) / dx  &
    + 180. * u(i-1) / dx -  90. * u(i-2) / dx  &
    +  20. * u(i-3) / dx )
p_x = -1./60.*( - 110. * p(i  ) / dx  &
    + 180. * p(i-1) / dx -  90. * p(i-2) / dx  &
    +  20. * p(i-3) / dx )
rho_x = -1./60.*( - 110. * rho(i  ) / dx  &
    + 180. * rho(i-1) / dx -  90. * rho(i-2) / dx  &
    +  20. * rho(i-3) / dx )

l1 = 0.0 !( Uc - sound_speed ) * ( p_x - rho(i) * sound_speed * u_x)
l2 = uc * ( p_x - sound_speed**2 * rho_x )
l3 = ( uc + sound_speed ) * ( p_x + rho(i) * sound_speed * u_x)
IF (it > 1) THEN
  l1 = -l3
END IF

d_1 = 0.5 / sound_speed**2 * ( l1 - 2. * l2 + l3 )
d_2 = 0.5 / sound_speed / rho(i) * ( l3 - l1 )
d_3 = 0.5 * ( l1 + l3 )

dedxi(i,1) = d_1
dedxi(i,2) = u(i) * d_1 + rho(i) * d_2
dedxi(i,3) = 0.5 * u(i)**2   * d_1 + rho(i) * u(i)   * d_2  &
    + 1. / (gamma-1.) * d_3

!...  <><><><><><><><><><>

END SUBROUTINE c_boundary_cond