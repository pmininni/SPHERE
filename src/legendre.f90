!=================================================================
! LEGENDRE subroutines
!
! Subroutines to compute the associated legendre polynomials 
! P_ml(x) and its derivatives. The polinomials are normalized 
! in the interval [-1,1] and the values for negative m are 
! computed according to the expressions in
!
! 2005 Pablo D. Mininni.
!      National Center for Atmospheric Research.
!      e-mail: mininni@ucar.edu 
!=================================================================

!*****************************************************************
      SUBROUTINE legder(l,m,x,plg,nor)
!-----------------------------------------------------------------
!
! Computes the derivative of the associated Legendre polynomial 
! sqrt(1-x^2)d P_ml(x)/d theta. To obtain the actual derivative 
! with respect to theta at the coordinate x, the result must be 
! divided by sqrt(1-x^2).
!
! Parameters
!    l  : integer >=0
!    m  : integer satisfying m<=l
!    x  : coordinate [=cos(theta)]
!    plg: at the output contains the derivative
!    nor: =0 unnormalized polynomial
!         =1 legendre normalization
!         =2 spherical harmonic normalization
!
      USE constants
      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(IN)  :: x
      DOUBLE PRECISION, INTENT(OUT) :: plg
      DOUBLE PRECISION              :: pl1,pl2
      DOUBLE PRECISION              :: norm,f1,f2
      INTEGER, INTENT(IN)           :: l,m
      INTEGER, INTENT(IN)           :: nor

      CALL plgndr(l,m,x,pl1,0)
      IF (l.gt.0.and.m.le.l-1) THEN
         CALL plgndr(l-1,m,x,pl2,0)
      ELSE
         pl2 = 0
      ENDIF
      plg = l*x*pl1-(l+m)*pl2
!
! Normalization
!
      IF (nor.gt.0) THEN
         IF (m.eq.0) THEN
            norm = SQRT(l+.5d0)
         ELSE IF (m.gt.0) THEN
            IF (l+m.le.32) THEN
               CALL factorial(l-m,f1)
               CALL factorial(l+m,f2)
               norm = SQRT((l+.5d0)*f1/f2)
            ELSE IF ((l+m.gt.32).and.(l-m.le.32)) THEN
               CALL factorial(l-m,f1)
               f1 = LOG(f1)
               CALL gammln(l+m+1d0,f2)
               norm = SQRT((l+.5d0)*EXP(f1-f2))
            ELSE
               CALL gammln(l-m+1d0,f1)
               CALL gammln(l+m+1d0,f2)
               norm = SQRT((l+.5d0)*EXP(f1-f2))
            ENDIF
         ENDIF
         IF (nor.eq.2) norm = norm/SQRT(2*PI)
         plg = plg*norm
      ENDIF

      RETURN
      END SUBROUTINE legder

!*****************************************************************
      SUBROUTINE plgndr(l,m,x,plg,nor)
!-----------------------------------------------------------------
!
! Computes the normalized associated Legendre polynomial 
! P_ml(x). This subroutine and the following ones are 
! based on the algorithms in Numerical Recipes. 
!
! Parameters
!    l  : integer >=0
!    m  : integer satisfying m<=l
!    x  : coordinate [=cos(theta)]
!    plg: at the output contains the value of P_ml(x)
!    nor: =0 unnormalized polynomial
!         =1 legendre normalization
!         =2 spherical harmonic normalization
!
      USE constants
      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(IN)  :: x
      DOUBLE PRECISION, INTENT(OUT) :: plg
      DOUBLE PRECISION              :: fact,pll,pmm
      DOUBLE PRECISION              :: pmmp1,somx2
      DOUBLE PRECISION              :: norm,f1,f2
      INTEGER, INTENT(IN)           :: l,m
      INTEGER, INTENT(IN)           :: nor
      INTEGER                       :: i,ll

      IF (m.lt.0.or.m.gt.l.or.ABS(x).gt.1.) THEN
         OPEN(1,file='error.log',position='append')
         WRITE(1,*) 'plgndr: bad arguments'
         CLOSE(1)
      ENDIF
      pmm = 1.d0
      IF (m.gt.0) THEN
         somx2 = SQRT((1.d0-x)*(1.d0+x))
         fact = 1.d0
 L1 :    DO i = 1,m
            pmm = -pmm*fact*somx2
            fact = fact+2.
         ENDDO L1
      ENDIF
      IF (l.eq.m) THEN
         plg = pmm
      ELSE
         pmmp1 = x*(2*m+1)*pmm
         IF (l.eq.m+1) THEN
            plg = pmmp1
         ELSE
 L2:        DO ll = m+2,l
               pll = (x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m)
               pmm = pmmp1
               pmmp1 = pll
            ENDDO L2
            plg = pll
         ENDIF
      ENDIF
!
! Normalization
!
      IF (nor.gt.0) THEN
         IF (m.eq.0) THEN
            norm = SQRT(l+.5d0)
         ELSE IF (m.gt.0) THEN
            IF (l+m.le.32) THEN
               CALL factorial(l-m,f1)
               CALL factorial(l+m,f2)
               norm = SQRT((l+.5d0)*f1/f2)
            ELSE IF ((l+m.gt.32).and.(l-m.le.32)) THEN
               CALL factorial(l-m,f1)
               f1 = LOG(f1)
               CALL gammln(l+m+1d0,f2)
               norm = SQRT((l+.5d0)*EXP(f1-f2))
            ELSE
               CALL gammln(l-m+1d0,f1)
               CALL gammln(l+m+1d0,f2)
               norm = SQRT((l+.5d0)*EXP(f1-f2))
            ENDIF
         ENDIF
         IF (nor.eq.2) norm = norm/SQRT(2*PI)
         plg = plg*norm
      ENDIF

      RETURN
      END SUBROUTINE plgndr

!*****************************************************************
      SUBROUTINE factorial(n,factrl)
!-----------------------------------------------------------------
!
! Computes the factorial of an integer in double precision. 
! For values of n smaller than 32, the factorial is computed 
! directly. In any other case, gamma functions are used.
!
! Parameters
!    n     : input variable
!    factrl: at the output contains the factorial
!
      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(OUT)   :: factrl
      DOUBLE PRECISION, SAVE          :: a(33)
      DOUBLE PRECISION                :: gamm
      INTEGER, INTENT(IN)             :: n
      INTEGER, SAVE                   :: ntop
      INTEGER                         :: j

      DATA ntop,a(1)/0,1./
      IF (n.lt.0) THEN
         OPEN(1,file='error.log',position='append')
         WRITE(1,*) 'factorial: negative input'
         CLOSE(1)
      ENDIF
      IF (n.le.ntop) THEN
         factrl = a(n+1)
      ELSE IF (n.le.32) THEN
 F1 :    DO j = ntop+1,n
            a(j+1) = j*a(j)
         ENDDO F1
         ntop = n
         factrl = a(n+1)
      ELSE
         CALL gammln(n+1d0,gamm)
         factrl = EXP(gamm)
      ENDIF

      RETURN
      END SUBROUTINE factorial

!*****************************************************************
      SUBROUTINE gammln(xx,gamm)
!-----------------------------------------------------------------
!
! Computes the natural logarithm of the Gamma function 
! Gamma(x), for x real larger than zero.
!
! Parameters
!    xx  : input variable
!    gamm: at the output contains ln[Gamma(xx)]
      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(IN)  :: xx
      DOUBLE PRECISION, INTENT(OUT) :: gamm
      DOUBLE PRECISION, SAVE        :: cof(6),stp
      DOUBLE PRECISION              :: ser,tmp,x,y
      INTEGER                       :: j

      DATA cof,stp/76.18009172947146d0,-86.50532032941677d0, &
           24.01409824083091d0,-1.231739572450155d0, &
           .1208650973866179d-2,-.5395239384953d-5, &
           2.5066282746310005d0/

      x = xx
      y = x
      tmp = x+5.5d0
      tmp = (x+0.5d0)*LOG(tmp)-tmp
      ser = 1.000000000190015d0
 G1 : DO j = 1,6
         y = y+1.d0
         ser = ser+cof(j)/y
      ENDDO G1
      gamm = tmp+LOG(stp*ser/x)

      RETURN
      END SUBROUTINE gammln
