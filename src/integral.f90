!=================================================================
! INTEGRAL subroutines
!
! Subroutines for numerical integration. The functions to 
! integrate should be declared as double precision functions 
! in a module 'functions'. Optional arguments can be passed 
! to these functions using a derived data type variable 
! 'arguments', also defined in the module 'functions'. Uses 
! the module constants to control the tolerance.
!
! 2005 Pablo D. Mininni.
!      National Center for Atmospheric Research.
!      e-mail: mininni@ucar.edu 
!=================================================================

!*****************************************************************
      SUBROUTINE qromo(func,a,b,ss,arg)
!-----------------------------------------------------------------
!
! Romberg integration on an open interval. Returns as ss the 
! integral of the function func from a to b. The function is 
! not evaluated at the end points. This subroutine and the 
! following ones are based on the algorithms in Numerical 
! Recipes. 
!
! Parameters
!    func: function name
!    a   : lower limit of the integral
!    b   : upper limit of the integral
!    ss  : at the output contains the integral
!    arg : derived type variable containing arguments for func
!
      USE functions
      USE constants
      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(IN)  :: a,b
      DOUBLE PRECISION, INTENT(OUT) :: ss
      DOUBLE PRECISION, EXTERNAL    :: func
      DOUBLE PRECISION              :: dss
      TYPE(ARGUMENTS), INTENT(IN)   :: arg
      INTEGER, PARAMETER            :: JMAX = 16
      INTEGER, PARAMETER            :: JMAXP = JMAX+1
      INTEGER, PARAMETER            :: K = 5
      INTEGER, PARAMETER            :: KM = K-1
      DOUBLE PRECISION              :: h(JMAXP),s(JMAXP)
      INTEGER                       :: j

      h(1) = 1.d0
 Q1 : DO j = 1,JMAX
         CALL midpnt(func,a,b,s(j),j,arg)
         IF (j.ge.K) THEN
            CALL polint(h(j-KM),s(j-KM),K,0.d0,ss,dss)
            IF (ABS(dss).le.TOL*ABS(ss)) RETURN
         ENDIF
         s(j+1) = s(j)
         h(j+1) = h(j)/9.d0
      ENDDO Q1
      OPEN(1,file='error.log',position='append')
      WRITE(1,*) 'qromo: reached maximum number of iterations'
      CLOSE(1)

      RETURN
      END SUBROUTINE qromo

!*****************************************************************
      SUBROUTINE midpnt(func,a,b,s,n,arg)
!-----------------------------------------------------------------
!
! Computes the nth stage of refinement of an extended midpoint 
! rule, to integrate the function func between limits a and b 
! using Romberg integration. When called with n=1, the routine 
! returns as s the crudest estimate of the integral of f(x) 
! between a and b. Subsequent calls with n=2,3,... (in sequential 
! order) will improve the accuracy of s by adding (2/3)×3^(n-1) 
! additional interior points. The variable s should not be 
! modified between sequential calls.
!
! Parameters
!    func: function name
!    a   : lower limit of the integral
!    b   : upper limit of the integral
!    s   : at the output contains the integral
!    n   : order of refinement
!    arg : derived type variable containing arguments for func
!
      USE functions
      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(IN)  :: a,b
      DOUBLE PRECISION, INTENT(OUT) :: s
      DOUBLE PRECISION, EXTERNAL    :: func
      DOUBLE PRECISION              :: ddel,del,sum,tnm,x
      TYPE(ARGUMENTS), INTENT(IN)   :: arg
      INTEGER, INTENT(IN)           :: n
      INTEGER                       :: it,j

      IF (n.eq.1) THEN
         s = (b-a)*func(0.5d0*(a+b),arg)
      ELSE
         it = 3**(n-2)
         tnm = it
         del = (b-a)/(3.d0*tnm)
         ddel = del+del
         x = a+0.5d0*del
         sum=0.
 M1 : DO j = 1,it
         sum = sum+func(x,arg)
         x = x+ddel
         sum = sum+func(x,arg)
         x = x+del
      ENDDO M1
      s = (s+(b-a)*sum/tnm)/3.d0
      ENDIF

      RETURN
      END SUBROUTINE midpnt

!*****************************************************************
      SUBROUTINE polint(xa,ya,n,x,y,dy)
!-----------------------------------------------------------------
!
! Polynomial interpolation and extrapolation.
! Given arrays xa and ya, each of length n, and given 
! a value x, this routine returns a value y, and an error 
! estimate dy. If P(x) is the polynomial of degree n-1 
! such that P(xa_i) = ya_i, then the returned value is 
! y = P(x).
!
! Parameters
!    xa: x points for interpolation
!    ya: y points for interpolation
!    n : order of the interpolation
!    x : coordinate value for the extrapolation
!    y : extrapolation
!    dy: error in the extrapolation
!
      IMPLICIT NONE

      INTEGER, INTENT(IN)           :: n
      DOUBLE PRECISION, INTENT(IN)  :: xa(n),ya(n)
      DOUBLE PRECISION, INTENT(IN)  :: x
      DOUBLE PRECISION, INTENT(OUT) :: y,dy

      INTEGER, PARAMETER            :: NMAX = 10
      INTEGER                       :: i,m,ns
      DOUBLE PRECISION              :: c(NMAX),d(NMAX)
      DOUBLE PRECISION              :: den,dif,dift
      DOUBLE PRECISION              :: ho,hp,w

      ns = 1
      dif = ABS(x-xa(1))
 P1 : DO i = 1,n
         dift = ABS(x-xa(i))
         IF (dift.lt.dif) THEN
            ns = i
            dif = dift
         ENDIF
         c(i) = ya(i)
         d(i) = ya(i)
      ENDDO P1
      y = ya(ns)
      ns = ns-1
 P2 : DO m = 1,n-1
 P3 :    DO i = 1,n-m
            ho = xa(i)-x
            hp = xa(i+m)-x
            w = c(i+1)-d(i)
            den = ho-hp
            IF (den.eq.0.) THEN
               OPEN(1,file='error.log',position='append')
               WRITE(1,*) 'polint: inputs identical to roundoff error'
               CLOSE(1)
            ENDIF
            den = w/den
            d(i) = hp*den
            c(i) = ho*den
         ENDDO P3
         IF (2*ns.lt.n-m) THEN
            dy = c(ns+1)
         ELSE
            dy = d(ns)
            ns = ns-1
         ENDIF
         y = y+dy
      ENDDO P2

      RETURN
      END SUBROUTINE polint

!*****************************************************************
      SUBROUTINE gauleg(x1,x2,x,w,n)
!-----------------------------------------------------------------
!
! Computes the abscissas and weights for Gauss-Legendre n-point 
! quadrature. The weigh function for Gauss-Legendre quadrature 
! is W(x) = 1, and the limits of integration are x1 and x2.
!
! Parameters
!    x1 : lower limit of the integral
!    x2 : upper limit of the integral
!    x  : at the output contains the abscissas
!    w  : at the output contains the weights
!    n  : number of abscissas and weights to compute
!
      USE constants
      IMPLICIT NONE

      INTEGER, INTENT (IN)          :: n
      DOUBLE PRECISION, INTENT(OUT) :: x(n),w(n)
      DOUBLE PRECISION, INTENT(IN)  :: x1,x2
      DOUBLE PRECISION, PARAMETER   :: EPS = 3.d-14
      DOUBLE PRECISION              :: p1,p2,p3,pp
      DOUBLE PRECISION              :: xl,xm,z,z1
      INTEGER                       :: i,j,m

      m = (n+1)/2
      xm = 0.5d0*(x2+x1)
      xl = 0.5d0*(x2-x1)
 G1 : DO i = 1,m
         z = COS(PI*(i-.25d0)/(n+.5d0))
 1       CONTINUE
         p1 = 1.d0
         p2 = 0.d0
 G2 :    DO j = 1,n
            p3 = p2
            p2 = p1
            p1 = ((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
         ENDDO G2
         pp = n*(z*p1-p2)/(z*z-1.d0)
         z1 = z
         z = z1-p1/pp
         IF (ABS(z-z1).gt.EPS) GOTO 1
         x(i) = xm-xl*z
         x(n+1-i) = xm+xl*z
         w(i) = 2.d0*xl/((1.d0-z*z)*pp*pp)
         w(n+1-i) = w(i)
      ENDDO G1

      RETURN
      END SUBROUTINE gauleg

!*****************************************************************
      SUBROUTINE gaujac(x,w,n,alf,bet)
!-----------------------------------------------------------------
!
! Computes the abscissas and weights for Gauss-Jacobi n-point 
! quadrature. The weigh function for Gauss-Jacobi quadrature 
! is W(x) = (1-x)^alpha*(1-x)^beta, and the limits of integration 
! are -1 and 1. The abscissas x(1:n) are ordered from the largest 
! to the smallest.
!
! Parameters
!    x  : at the output contains the abscissas
!    w  : at the output contains the weights
!    n  : number of abscissas and weights to compute
!    alf: value of alpha in the weight function
!    bet: value of beta in the weight function
! 
      IMPLICIT NONE

      INTEGER, INTENT(IN)           :: n
      DOUBLE PRECISION, INTENT(OUT) :: x(n),w(n)
      DOUBLE PRECISION, INTENT(IN)  :: alf,bet
      DOUBLE PRECISION, PARAMETER   :: EPS = 3.d-14      
      DOUBLE PRECISION              :: alfbet,an,bn
      DOUBLE PRECISION              :: r1,r2,r3,gam1
      DOUBLE PRECISION              :: gam2,gam3,gam4
      DOUBLE PRECISION              :: a,b,c,p1,p2,p3
      DOUBLE PRECISION              :: pp,temp,z,z1
      INTEGER, PARAMETER            :: MAXIT = 10
      INTEGER                       :: i,its,j

 J1 : DO i = 1,n
         IF (i.eq.1) THEN
            an = alf/n
            bn = bet/n
            r1 = (1.+alf)*(2.78/(4.+n*n)+.768*an/n)
            r2 = 1.+1.48*an+.96*bn+.452*an*an+.83*an*bn
            z = 1.-r1/r2
         ELSE IF (i.eq.2) THEN
            r1 = (4.1+alf)/((1.+alf)*(1.+.156*alf))
            r2 = 1.+.06*(n-8.)*(1.+.12*alf)/n
            r3 = 1.+.012*bet*(1.+.25*abs(alf))/n
            z = z-(1.-z)*r1*r2*r3
         ELSE IF (i.eq.3) THEN
            r1 = (1.67+.28*alf)/(1.+.37*alf)
            r2 = 1.+.22*(n-8.)/n
            r3 = 1.+8.*bet/((6.28+bet)*n*n)
            z = z-(x(1)-z)*r1*r2*r3
         ELSE IF (i.eq.n-1) THEN
            r1 = (1.+.235*bet)/(.766+.119*bet)
            r2 = 1./(1.+.639*(n-4.)/(1.+.71*(n-4.)))
            r3 = 1./(1.+20.*alf/((7.5+alf)*n*n))
            z = z+(z-x(n-3))*r1*r2*r3
         ELSE IF (i.eq.n) THEN
            r1 = (1.+.37*bet)/(1.67+.28*bet)
            r2 = 1./(1.+.22*(n-8.)/n)
            r3 = 1./(1.+8.*alf/((6.28+alf)*n*n))
            z = z+(z-x(n-2))*r1*r2*r3
         ELSE
            z = 3.*x(i-1)-3.*x(i-2)+x(i-3)
         ENDIF
         alfbet = alf+bet
 J2 :    DO its = 1,MAXIT
            temp = 2.d0+alfbet
            p1 = (alf-bet+temp*z)/2.d0
            p2 = 1.d0
 J3 :       DO j = 2,n
               p3 = p2
               p2 = p1
               temp = 2*j+alfbet
               a = 2*j*(j+alfbet)*(temp-2.d0)
               b = (temp-1.d0)*(alf*alf-bet*bet+temp*(temp-2.d0)*z)
               c = 2.d0*(j-1+alf)*(j-1+bet)*temp
               p1 = (b*p2-c*p3)/a
            ENDDO J3
            pp = (n*(alf-bet-temp*z)*p1+2.d0*(n+alf)* &
                 (n+bet)*p2)/(temp*(1.d0-z*z))
            z1 = z
           z = z1-p1/pp
           IF (abs(z-z1).le.EPS) GOTO 2
           ENDDO J2
           OPEN(1,file='error.log',position='append')
           WRITE(1,*) 'gaujac: too many iterations'
           CLOSE(1)
 2         x(i) = z
           CALL gammln(alf+n,gam1)
           CALL gammln(bet+n,gam2)
           CALL gammln(n+1.d0,gam3)
           CALL gammln(n+alfbet+1.d0,gam4)
           w(i) = EXP(gam1+gam2-gam3-gam4)*temp*2.**alfbet/(pp*p2)
      ENDDO J1

      RETURN
      END SUBROUTINE gaujac

!*****************************************************************
      SUBROUTINE gauchv(x,w,n)
!-----------------------------------------------------------------
!
! Computes the abscissas and weights for Gauss-Chebyshev n-point 
! quadrature. The weigh function for Gauss-Chebyshev quadrature 
! is W(x) = 1/sqrt(1-x^2), and the limits of integration are -1 
! and 1.
!
! Parameters
!    x  : at the output contains the abscissas
!    w  : at the output contains the weights
!    n  : number of abscissas and weights to compute
! 
      USE constants
      IMPLICIT NONE

      INTEGER, INTENT(IN)           :: n
      DOUBLE PRECISION, INTENT(OUT) :: x(n),w(n)
      INTEGER                       :: i

      DO i = 1,n
         x(i) = cos(PI*(i-.5d0)/n)
         w(i) = PI/n
      END DO

      RETURN
      END SUBROUTINE gauchv
