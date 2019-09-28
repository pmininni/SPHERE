!=================================================================
! BESSEL subroutines
!
! Subroutines for computing BesselJ and Besselj functions and 
! its derivatives and zeros.
!
! 2005 Pablo D. Mininni.
!      National Center for Atmospheric Research.
!      e-mail: mininni@ucar.edu 
!=================================================================

!*****************************************************************
       SUBROUTINE zerobes(n,x1,x2,root)
!-----------------------------------------------------------------
!
! Computes the roots of the spherical Besselj functions between 
! x1 and x2. Uses a combination of bisection and Newton-Raphson 
! to find the zeros. The error in the root is controled by the 
! variable TOL in the module constants.
!
! Parameters
!    n   : order of the Bessel function (>=0)
!    x1  : beginning of the interval
!    x2  : end of the interval
!    root: at the output contains the root
!
      USE constants
      IMPLICIT NONE

      DOUBLE PRECISION, intent(IN)  :: x1,x2
      DOUBLE PRECISION, intent(OUT) :: root
      DOUBLE PRECISION              :: df,dx,dxold,f,fh
      DOUBLE PRECISION              :: fl,temp,xh,xl,order
      INTEGER, INTENT(IN)           :: n
      INTEGER, PARAMETER            :: MAXIT = 100 
      INTEGER                       :: j

      order = n+0.5d0
      CALL bessj(x1,order,fl,df)
      CALL bessj(x2,order,fh,df)
      IF ((fl.gt.0..and.fh.gt.0.).or.(fl.lt.0..and.fh.lt.0.)) THEN
         OPEN(1,file='error.log',position='append')
         WRITE(1,*) 'zerobes: interval contains no change of sign'
         CLOSE(1)
      ENDIF
      IF (ABS(fl).eq.0.) THEN
         root = x1
         RETURN
      ELSE IF (ABS(fh).eq.0.) THEN
         root = x2
         RETURN
      ELSE IF(fl.lt.0.) THEN
         xl = x1
         xh = x2
      ELSE
         xh = x1
         xl = x2
      ENDIF
      root = .5*(x1+x2)
      dxold = ABS(x2-x1)
      dx = dxold
      CALL bessj(root,order,f,df)
 Z1 : DO j = 1,MAXIT
         IF (((root-xh)*df-f)*((root-xl)*df-f).gt.0. &
            .or. abs(2.*f).gt.abs(dxold*df)) THEN
            dxold = dx
            dx = 0.5*(xh-xl)
            root = xl+dx
            IF (xl.eq.root) RETURN
         ELSE
            dxold = dx
            dx = f/df
            temp = root
            root = root-dx
            IF (temp.eq.root) RETURN
         ENDIF
         IF (abs(dx).lt.TOL) RETURN
         CALL bessj(root,order,f,df)
         IF (f.lt.0.) THEN
            xl = root
         ELSE
            xh = root
         ENDIF
      ENDDO Z1
      OPEN(1,file='error.log',position='append')
      WRITE(1,*) 'zerobes: root exceeding maximum iterations'
      CLOSE(1)

      RETURN
      END SUBROUTINE zerobes

!*****************************************************************
      SUBROUTINE sphbes(n,x,sj,sjp)
!-----------------------------------------------------------------
!
! Spherical Besselj function and its derivative. This subroutine 
! and the following ones were modified from Numerical Recipes.
!
! Parameters
!    x  : argument of the Bessel function (> = 0)
!    n  : order of the Bessel function (>=0)
!    sj : Bessel function j_n(x)
!    sjp: Derivative of the Bessel function
!
      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(IN)  :: x
      DOUBLE PRECISION, INTENT(OUT) :: sj,sjp
      DOUBLE PRECISION              :: factor,order,rj,rjp
      DOUBLE PRECISION, PARAMETER   :: RTPIO2 = 1.253314137315500d0
      INTEGER, INTENT(IN)           :: n

      order = n+0.5d0
      CALL bessj(x,order,rj,rjp)
      factor = RTPIO2/SQRT(x)
      sj = factor*rj
      sjp = factor*rjp-sj/(2.*x)

      RETURN
      END SUBROUTINE sphbes

!*****************************************************************
      SUBROUTINE bessj(x,xnu,rj,rjp)
!-----------------------------------------------------------------
!
! BesselJ function and its derivative. EPS controls the accuracy 
! and FPMIN is a number close to the smallest floating-point 
! number in the machine.
!
! Parameters
!    x  : argument of the Bessel function (>=0)
!    xnu: order of the Bessel function (>=0)
!    rj : Bessel function J_xnu(x)
!    rjp: Derivative of the Bessel function
!
      USE constants
      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(IN)  :: x,xnu
      DOUBLE PRECISION, INTENT(OUT) :: rj,rjp
      DOUBLE PRECISION, PARAMETER   :: EPS = 1.e-16
      DOUBLE PRECISION, PARAMETER   :: FPMIN = 1.e-30

      DOUBLE PRECISION   :: a,b,br,bi,c,cr,ci,d,del,den,di
      DOUBLE PRECISION   :: dlr,dli,dr,f,fact,gam,h,p,q,r
      DOUBLE PRECISION   :: rjl,rjl1,rjmu,rjp1,rjpl,rjtemp
      DOUBLE PRECISION   :: temp,w,xi,xi2,xmu,xmu2

      INTEGER, PARAMETER :: MAXIT = 10000
      INTEGER            :: i,isign,l,nl

      nl = MAX(0,INT(xnu-x+1.5d0))
      xmu = xnu-nl
      xmu2 = xmu*xmu
      xi = 1.d0/x
      xi2 = 2.d0*xi
      w = xi2/PI
      isign = 1
      h = xnu*xi
      IF (h.lt.FPMIN) h = FPMIN
      b = xi2*xnu
      d = 0.d0
      c = h
 I1 : DO i = 1,MAXIT
         b = b+xi2
         d = b-d
         IF (ABS(d).lt.FPMIN) d = FPMIN
         c = b-1.d0/c
         IF (ABS(c).lt.FPMIN) c=FPMIN
         d = 1.d0/d
         del = c*d
         h = del*h
         IF (d.lt.0.d0) isign = -isign
         IF (ABS(del-1.d0).lt.EPS) GOTO 10
      ENDDO I1
      OPEN(1,file='error.log',position='append')
      WRITE(1,*) 'bessj: x too large, try asymptotic expansion'
      CLOSE(1)
 10   CONTINUE
      rjl = isign*FPMIN
      rjpl = h*rjl
      rjl1 = rjl
      rjp1 = rjpl
      fact = xnu*xi
 I2 : DO l = nl,1,-1
         rjtemp = fact*rjl+rjpl
         fact = fact-xi
         rjpl = fact*rjtemp-rjl
         rjl = rjtemp
      ENDDO I2
      IF (rjl.eq.0.d0) rjl = EPS
      f = rjpl/rjl
      a = .25d0-xmu2
      p = -.5d0*xi
      q = 1.d0
      br = 2.d0*x
      bi = 2.d0
      fact = a*xi/(p*p+q*q)
      cr = br+q*fact
      ci = bi+p*fact
      den = br*br+bi*bi
      dr = br/den
      di = -bi/den
      dlr = cr*dr-ci*di
      dli = cr*di+ci*dr
      temp = p*dlr-q*dli
      q = p*dli+q*dlr
      p = temp
 I3 : DO i=2,MAXIT
         a = a+2*(i-1)
         bi = bi+2.d0
         dr = a*dr+br
         di = a*di+bi
         IF (ABS(dr)+ABS(di).lt.FPMIN) dr = FPMIN
         fact = a/(cr*cr+ci*ci)
         cr = br+cr*fact
         ci = bi-ci*fact
         IF (ABS(cr)+ABS(ci).lt.FPMIN) cr = FPMIN
         den = dr*dr+di*di
         dr = dr/den
         di = -di/den
         dlr = cr*dr-ci*di
         dli = cr*di+ci*dr
         temp = p*dlr-q*dli
         q = p*dli+q*dlr
         p = temp
         IF (ABS(dlr-1.d0)+ABS(dli).lt.EPS) GOTO 30
      ENDDO I3
      OPEN(1,file='error.log',position='append')
      WRITE(1,*) 'bessj: CF2 failed'
      CLOSE(1)
 30   CONTINUE
      gam = (p-f)/q
      rjmu = SQRT(w/((p-f)*gam+q))
      rjmu = SIGN(rjmu,rjl)
      fact = rjmu/rjl
      rj = rjl1*fact
      rjp = rjp1*fact

      RETURN
      END SUBROUTINE bessj
