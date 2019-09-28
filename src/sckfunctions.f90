!=================================================================
! SCKFUNCTIONS subroutines
!
! Subroutines to compute spherical Chandrasekhar-Kendall 
! functions for the SPHERE code.
!
! 2005 Pablo D. Mininni.
!      National Center for Atmospheric Research.
!      e-mail: mininni@ucar.edu 
!=================================================================

!*****************************************************************
       SUBROUTINE xsi2ck(xsi,lambda,normal,d,n,odir,ext,opt)
!-----------------------------------------------------------------
!
! Creates arrays with the values of the three cartesian 
! components of the field, in a cartesian grid. The box 
! has edges of size 2d, and the sphere inside the box 
! has radius 1. The arrays are exported to binary files 
! 
! Parameters
!    xsi   : amplitude of the Chandrasekhar-Kendall modes
!    lambda: table with the roots of the bessel functions
!    normal: normalization coefficients
!    d     : edge of the box (d>=2)
!    n     : box resolution in each direction
!    odir  : directory for binary output
!    ext   : extension used for the file
!    opt   : =0 velocity field
!            =1 magnetic field
!
      USE constants
      USE resolution
      IMPLICIT NONE

      DOUBLE COMPLEX, INTENT(IN)     :: xsi(l*(l+3)/2,2*q)
      DOUBLE COMPLEX                 :: eimp
      DOUBLE PRECISION, INTENT(IN)   :: lambda(l,q)
      DOUBLE PRECISION, INTENT(IN)   :: normal(l,q)
      DOUBLE PRECISION, INTENT(IN)   :: d
      DOUBLE PRECISION               :: lam1,norm
      DOUBLE PRECISION               :: plg,plgp
      DOUBLE PRECISION               :: sj,sjp
      DOUBLE PRECISION               :: fr,ft,ff
      DOUBLE PRECISION               :: ampr
      DOUBLE PRECISION               :: am1t,am2t,am1f,am2f
      DOUBLE PRECISION               :: sum1,sum2,sum3
      DOUBLE PRECISION               :: x,y,z
      DOUBLE PRECISION               :: ere,rho
      DOUBLE PRECISION               :: cth,sth
      DOUBLE PRECISION               :: cph,sph
      REAL                           :: fx(n,n,n)
      REAL                           :: fy(n,n,n)
      REAL                           :: fz(n,n,n)
      INTEGER, INTENT(IN)            :: n
      INTEGER, INTENT(IN)            :: opt
      INTEGER                        :: i,j,k
      INTEGER                        :: q1,l1,m1
      INTEGER                        :: ind1
      CHARACTER(len=100), INTENT(IN) :: odir
      CHARACTER(len=3), INTENT(IN)   :: ext
      CHARACTER                      :: field

      DO k = 1,n
         z = d*((DBLE(k)-1)/(DBLE(n)-1)-.5d0)
      DO j = 1,n
         y = d*((DBLE(j)-1)/(DBLE(n)-1)-.5d0)
      DO i = 1,n
         x = d*((DBLE(i)-1)/(DBLE(n)-1)-.5d0)
         ere = SQRT(x**2+y**2+z**2)
         rho = SQRT(x**2+y**2)
         IF (ere.lt.tiny) ere = tiny
         IF (rho.lt.tiny) rho = tiny
 SP :    IF (ere.le.1) THEN             ! points inside the sphere
         cth = z/ere                    ! cos(theta)
         sth = rho/ere                  ! sin(theta)
         cph = x/rho                    ! cos(phi)
         sph = y/rho                    ! sin(phi)
         fr = 0.d0
         ft = 0.d0
         ff = 0.d0
 SQ :    DO q1 = 1,2*q                  ! sum over q
 SL :    DO l1 = 1,l                    ! sum over l
            IF (q1.le.q) THEN
               lam1 = lambda(l1,q1)
               norm = normal(l1,q1)
            ELSE
               lam1 = -lambda(l1,2*q-q1+1)
               norm = normal(l1,2*q-q1+1)
            ENDIF
            CALL sphbes(l1,ABS(lam1)*ere,sj,sjp)
            sjp = sjp*ABS(lam1)*ere+sj
            ampr = l1*(l1+1)*norm*sj/ere
            am1t = -2*norm*lam1*sj/sth
            am2t = norm*sjp/ere
            am1f = -norm*lam1*sj
            am2f = -2*norm*sjp/(ere*sth)
            ind1 = l1*(l1+1)/2          ! m =0
            CALL plgndr(l1,0,cth,plg,2)
            CALL legder(l1,0,cth,plgp,2)
            IF (sth.ge.tiny) THEN
               plgp = plgp/sth
            ELSE
               plgp = 0.d0
            ENDIF
            sum1 = xsi(ind1,q1)*plg
            sum2 = xsi(ind1,q1)*plgp
            sum3 = 0.
 SM :       DO m1 = 1,l1                ! sum over 1<=m<=l
               ind1 = l1*(l1+1)/2+m1
               eimp = (cph+IM*sph)**m1  ! De Moivre's formula for e^(im.phi)
               CALL plgndr(l1,m1,cth,plg,2)
               CALL legder(l1,m1,cth,plgp,2)
               IF (sth.ge.tiny) THEN
                  plgp = plgp/sth
               ELSE
                  plgp = 0.d0
               ENDIF
               eimp = eimp*xsi(ind1,q1)
               sum1 = sum1+2*plg*DBLE(eimp)
               sum2 = sum2+2*plgp*DBLE(eimp)
               sum3 = sum3+m1*plg*AIMAG(eimp)
            END DO SM
            fr = fr+ampr*sum1
            ft = ft+am1t*sum3+am2t*sum2
            ff = ff+am1f*sum2+am2f*sum3
         END DO SL
         END DO SQ                      ! cartesian projection
            fx(i,j,k) = REAL(fr*sth*cph+ft*cth*cph-ff*sph)
            fy(i,j,k) = REAL(fr*sth*sph+ft*cth*sph+ff*cph)
            fz(i,j,k) = REAL(fr*cth-ft*sth)
         ELSE                           ! points outside the sphere
            fx(i,j,k) = 0.
            fy(i,j,k) = 0.
            fz(i,j,k) = 0.
         ENDIF SP
      END DO
      END DO
      END DO

      IF (opt.eq.0) THEN                ! exports the arrays
         field = 'v'
      ELSE
         field = 'b'
      ENDIF
      OPEN(1,file=trim(odir)//'/field_'//field//'x.'//ext//'.out', &
          form='unformatted')
      WRITE(1) fx
      CLOSE(1)
      OPEN(1,file=trim(odir)//'/field_'//field//'y.'//ext//'.out', &
          form='unformatted')
      WRITE(1) fy
      CLOSE(1)
      OPEN(1,file=trim(odir)//'/field_'//field//'z.'//ext//'.out', &
          form='unformatted')
      WRITE(1) fz
      CLOSE(1)

      RETURN
      END SUBROUTINE xsi2ck

!*****************************************************************
       SUBROUTINE xsi2ckp(xsi,lambda,normal,d,n,nx,ny,nz,opt)
!-----------------------------------------------------------------
!
! Computes the value of the three cartesian components of 
! a field in one point in real space, in a cartesian grid. 
! The box has edges of size 2d, and the sphere inside the 
! box has radius 1. The result is appended to a text file.
!
! Parameters
!    xsi   : amplitude of the Chandrasekhar-Kendall modes
!    lambda: table with the roots of the bessel functions
!    normal: normalization coefficients
!    d     : edge of the box (d>=2)
!    n     : box resolution in each direction
!    nx    : x coordinate where the field is computed
!    ny    : y coordinate where the field is computed
!    nz    : z coordinate where the field is computed
!    opt   : =0 velocity field
!            =1 magnetic field
!
      USE constants
      USE resolution
      IMPLICIT NONE

      DOUBLE COMPLEX, INTENT(IN)     :: xsi(l*(l+3)/2,2*q)
      DOUBLE COMPLEX                 :: eimp
      DOUBLE PRECISION, INTENT(IN)   :: lambda(l,q)
      DOUBLE PRECISION, INTENT(IN)   :: normal(l,q)
      DOUBLE PRECISION, INTENT(IN)   :: d
      DOUBLE PRECISION               :: lam1,norm
      DOUBLE PRECISION               :: plg,plgp
      DOUBLE PRECISION               :: sj,sjp
      DOUBLE PRECISION               :: fr,ft,ff
      DOUBLE PRECISION               :: ampr
      DOUBLE PRECISION               :: am1t,am2t,am1f,am2f
      DOUBLE PRECISION               :: sum1,sum2,sum3
      DOUBLE PRECISION               :: x,y,z
      DOUBLE PRECISION               :: ere,rho
      DOUBLE PRECISION               :: cth,sth
      DOUBLE PRECISION               :: cph,sph
      REAL                           :: fx
      REAL                           :: fy
      REAL                           :: fz
      INTEGER, INTENT(IN)            :: n
      INTEGER, INTENT(IN)            :: nx,ny,nz
      INTEGER, INTENT(IN)            :: opt
      INTEGER                        :: q1,l1,m1
      INTEGER                        :: ind1
      CHARACTER                      :: field

      z = d*((DBLE(nz)-1)/(DBLE(n)-1)-.5d0)
      y = d*((DBLE(ny)-1)/(DBLE(n)-1)-.5d0)
      x = d*((DBLE(nx)-1)/(DBLE(n)-1)-.5d0)
      ere = SQRT(x**2+y**2+z**2)
      rho = SQRT(x**2+y**2)
      IF (ere.lt.tiny) ere = tiny
      IF (rho.lt.tiny) rho = tiny
 SP : IF (ere.le.1) THEN                ! point inside the sphere
         cth = z/ere                    ! cos(theta)
         sth = rho/ere                  ! sin(theta)
         cph = x/rho                    ! cos(phi)
         sph = y/rho                    ! sin(phi)
         fr = 0.d0
         ft = 0.d0
         ff = 0.d0
 SQ :    DO q1 = 1,2*q                  ! sum over q
 SL :    DO l1 = 1,l                    ! sum over l
            IF (q1.le.q) THEN
               lam1 = lambda(l1,q1)
               norm = normal(l1,q1)
            ELSE
               lam1 = -lambda(l1,2*q-q1+1)
               norm = normal(l1,2*q-q1+1)
            ENDIF
            CALL sphbes(l1,ABS(lam1)*ere,sj,sjp)
            sjp = sjp*ABS(lam1)*ere+sj
            ampr = l1*(l1+1)*norm*sj/ere
            am1t = -2*norm*lam1*sj/sth
            am2t = norm*sjp/ere
            am1f = -norm*lam1*sj
            am2f = -2*norm*sjp/(ere*sth)
            ind1 = l1*(l1+1)/2          ! m =0
            CALL plgndr(l1,0,cth,plg,2)
            CALL legder(l1,0,cth,plgp,2)
            IF (sth.ge.tiny) THEN
               plgp = plgp/sth
            ELSE
               plgp = 0.d0
            ENDIF
            sum1 = xsi(ind1,q1)*plg
            sum2 = xsi(ind1,q1)*plgp
            sum3 = 0.
 SM :       DO m1 = 1,l1                ! sum over 1<=m<=l
               ind1 = l1*(l1+1)/2+m1
               eimp = (cph+IM*sph)**m1  ! De Moivre's formula for e^(im.phi)
               CALL plgndr(l1,m1,cth,plg,2)
               CALL legder(l1,m1,cth,plgp,2)
               IF (sth.ge.tiny) THEN
                  plgp = plgp/sth
               ELSE
                  plgp = 0.d0
               ENDIF
               eimp = eimp*xsi(ind1,q1)
               sum1 = sum1+2*plg*DBLE(eimp)
               sum2 = sum2+2*plgp*DBLE(eimp)
               sum3 = sum3+m1*plg*AIMAG(eimp)
            END DO SM
            fr = fr+ampr*sum1
            ft = ft+am1t*sum3+am2t*sum2
            ff = ff+am1f*sum2+am2f*sum3
         END DO SL
         END DO SQ                      ! cartesian projection
            fx = REAL(fr*sth*cph+ft*cth*cph-ff*sph)
            fy = REAL(fr*sth*sph+ft*cth*sph+ff*cph)
            fz = REAL(fr*cth-ft*sth)
      ELSE                              ! point outside the sphere
         fx = 0.
         fy = 0.
         fz = 0.
      ENDIF SP

      IF (opt.eq.0) THEN                ! writes the result
         field = 'v'
      ELSE
         field = 'b'
      ENDIF
      OPEN(1,file='field_'//field//'.txt',position='append')
      WRITE(1,10) fx,fy,fz
   10 FORMAT( E22.14,E22.14,E22.14 )
      CLOSE(1)

      RETURN
      END SUBROUTINE xsi2ckp
