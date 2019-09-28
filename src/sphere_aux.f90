!=================================================================
! SPHERE subroutines
!
! Auxiliary subroutines to compute energy, helicity, enstrophy, 
! spectra, and other global quantities in the hydrodynamic and 
! magnetohydrodynamic SPHERE spectral codes. 
!
! 2005 Pablo D. Mininni.
!      National Center for Atmospheric Research.
!      e-mail: mininni@ucar.edu 
!=================================================================

!*****************************************************************
      SUBROUTINE global(t,lambda,xsi,opt)
!-----------------------------------------------------------------
!
! Computes the total energy, helicity, and enstropy of the 
! field with expanding coefficients xsi. The output is 
! appended to the file '*global.txt' where * is 'k' or 'm' 
! according to the value of the input variable opt. In both 
! cases, the file contains four columns with time, energy, 
! helicity and enstrophy.
!
! Parameters
!    t     : time
!    lambda: table with the roots of the Besselj function
!    xsi   : expanding coefficients for the field
!    opt   : =0 kinetic energy
!            =1 magnetic energy
!
      USE resolution
      IMPLICIT NONE

      DOUBLE COMPLEX, INTENT(IN)     :: xsi(l*(l+3)/2,2*q)
      DOUBLE PRECISION, INTENT(IN)   :: lambda(l,q)
      DOUBLE PRECISION, INTENT(IN)   :: t
      DOUBLE PRECISION               :: ene,hel,ens
      DOUBLE PRECISION               :: lam1,lam2,tmp
      INTEGER, INTENT(IN)            :: opt
      INTEGER                        :: q1,l1,m1

      ene = 0.d0
      hel = 0.d0
      ens = 0.d0
      IF (opt.eq.0) THEN                ! computes kinetic quantities
         DO q1 = 1,2*q
         DO l1 = 1,l
            IF (q1.le.q) THEN
               lam1 = lambda(l1,q1)
            ELSE
               lam1 = -lambda(l1,2*q-q1+1)
            ENDIF
            lam2 = lam1**2
            tmp = ABS(xsi(l1*(l1+1)/2,q1))**2
            ene = ene+.5d0*tmp
            hel = hel+.5d0*tmp*lam1
            ens = ens+.5d0*tmp*lam2
            DO m1 = 1,l1
               tmp = ABS(xsi(l1*(l1+1)/2+m1,q1))**2
               ene = ene+tmp
               hel = hel+tmp*lam1
               ens = ens+tmp*lam2
            ENDDO
         ENDDO
         ENDDO
      ELSE                              ! computes magnetic quantities
         DO q1 = 1,2*q
         DO l1 = 1,l
            IF (q1.le.q) THEN
               lam1 = 1.d0/lambda(l1,q1)
               lam2 = lambda(l1,q1)**2
            ELSE
               lam1 = -1.d0/lambda(l1,2*q-q1+1)
               lam2 = lambda(l1,2*q-q1+1)**2
            ENDIF
            tmp = ABS(xsi(l1*(l1+1)/2,q1))**2
            ene = ene+.5d0*tmp
            hel = hel+.5d0*tmp*lam1
            ens = ens+.5d0*tmp*lam2
            DO m1 = 1,l1
               tmp = ABS(xsi(l1*(l1+1)/2+m1,q1))**2
               ene = ene+tmp
               hel = hel+tmp*lam1
               ens = ens+tmp*lam2
            ENDDO
         ENDDO
         ENDDO
      ENDIF

      IF (opt.eq.0) THEN                ! exports the results
         OPEN(1,file='kglobal.txt',position='append')
      ELSE
         OPEN(1,file='mglobal.txt',position='append')
      ENDIF
      WRITE(1,10) t,ene,hel,ens
   10 FORMAT( E13.6,E22.14,E22.14,E22.14 )
      CLOSE(1)

      RETURN
      END SUBROUTINE global

!*****************************************************************
      SUBROUTINE injections(t,lambda,xsi,normal,forv,wtx,wty,wtz)
!-----------------------------------------------------------------
!
! Computes the total energy, helicity, and enstropy of the 
! field with expanding coefficients xsi. The output is 
! appended to the file '*global.txt' where * is 'k' or 'm' 
! according to the value of the input variable opt. In both 
! cases, the file contains four columns with time, energy, 
! helicity and enstrophy.
!
! Parameters
!    t     : time
!    lambda: table with the roots of the Besselj function
!    xsi   : expanding coefficients for the field
!
      USE constants
      USE resolution
      IMPLICIT NONE

      DOUBLE COMPLEX, INTENT(IN)     :: xsi(l*(l+3)/2,2*q)
      DOUBLE COMPLEX, INTENT(IN)     :: forv(l*(l+3)/2,2*q)
      DOUBLE PRECISION, INTENT(IN)   :: normal(l,q)
      DOUBLE PRECISION, INTENT(IN)   :: lambda(l,q)
      DOUBLE PRECISION, INTENT(IN)   :: t
      DOUBLE PRECISION, INTENT(IN)   :: wtx,wty,wtz
      DOUBLE PRECISION               :: prec,forc
      DOUBLE PRECISION               :: lam1,lam2,tmp
      INTEGER                        :: q1,l1,m1
      DOUBLE PRECISION               :: lx,ly,lz
      DOUBLE PRECISION               :: norm
      DOUBLE PRECISION               :: ampl,stwo
      DOUBLE PRECISION               :: sj,sjp

      forc = 0.d0
      DO q1 = 1,2*q
      DO l1 = 1,l
         IF (q1.le.q) THEN
            lam1 = lambda(l1,q1)
         ELSE
            lam1 = -lambda(l1,2*q-q1+1)
         ENDIF
         tmp = CONJG(xsi(l1*(l1+1)/2,q1))*forv(l1*(l1+1)/2,q1)
         forc = forc+tmp
         DO m1 = 1,l1
            tmp = CONJG(xsi(l1*(l1+1)/2+m1,q1))*forv(l1*(l1+1)/2+m1,q1)
            forc = forc+2.d0*tmp
         ENDDO
      ENDDO
      ENDDO

      lx = 0.d0
      ly = 0.d0
      lz = 0.d0
      stwo = SQRT(2.d0)
      DO q1 = 1,2*q                     ! sum over q
         IF (q1.le.q) THEN
            lam1 = lambda(1,q1)
            norm = normal(1,q1)
         ELSE
            lam1 = -lambda(1,2*q-q1+1)
            norm = normal(1,2*q-q1+1)
         ENDIF
         CALL sphbes(1,ABS(lam1),sj,sjp)
         ampl = 4*SQRT(PI/3)*norm*sjp*ABS(lam1)/lam1
         lz = lz-ampl*xsi(1,q1)
         lx = lx+ampl*stwo*DBLE(xsi(2,q1))
         ly = ly-ampl*stwo*AIMAG(xsi(2,q1))
      ENDDO
      prec = wtx*lx+wty*ly+wtz*lz

      OPEN(1,file='injections.txt',position='append')
      WRITE(1,10) t,prec,forc
   10 FORMAT( E13.6,E22.14,E22.14 )
      CLOSE(1)

      RETURN
      END SUBROUTINE injections

!*****************************************************************
      SUBROUTINE momentum(lambda,normal,xsi,opt)
!-----------------------------------------------------------------
!
! Exports the angular momentum of the velocity field, or 
! the dipole moment of the magnetic field. The fields have 
! expanding coefficients xsi. The output is appended to the 
! file 'lmoment.txt' or 'dmoment.txt'. The files contain 
! three columns with the x, y, and z cartesian components 
! of the vector.
!
! Parameters
!    lambda: table with the roots of the Besselj function
!    lambda: table with the normalization coefficients
!    xsi   : expanding coefficients for the field
!    opt   : =0 angular momentum
!            =1 dipole moment
!
      USE constants
      USE resolution
      IMPLICIT NONE

      DOUBLE COMPLEX, INTENT(IN)     :: xsi(l*(l+3)/2,2*q)
      DOUBLE PRECISION, INTENT(IN)   :: lambda(l,q)
      DOUBLE PRECISION, INTENT(IN)   :: normal(l,q)
      DOUBLE PRECISION               :: norm,lam1
      DOUBLE PRECISION               :: ampl,stwo
      DOUBLE PRECISION               :: lx,ly,lz
      DOUBLE PRECISION               :: sj,sjp
      INTEGER, INTENT(IN)            :: opt
      INTEGER                        :: q1

      CALL momenaux(lambda,normal,xsi,lx,ly,lz,opt)
      IF (opt.eq.0) THEN                ! exports the results
         OPEN(1,file='lmoment.txt',position='append')
      ELSE
         OPEN(1,file='dmoment.txt',position='append')
      ENDIF
      WRITE(1,20) lx,ly,lz
   20 FORMAT( E22.14,E22.14,E22.14 )
      CLOSE(1)

      RETURN
      END SUBROUTINE momentum

!*****************************************************************
      SUBROUTINE momenaux(lambda,normal,xsi,lx,ly,lz,opt)
!-----------------------------------------------------------------
!
! Computes the three cartesian components of the angular 
! momentum of the velocity field, or the dipole moment of 
! the magnetic field. The fields have expanding coefficients 
! xsi.
!
! Parameters
!    lambda: table with the roots of the Besselj function
!    normal: table with the normalization coefficients
!    xsi   : expanding coefficients for the field
!    lx    : x component of the momentum
!    ly    : y component of the momentum
!    lz    : z component of the momentum
!    opt   : =0 angular momentum
!            =1 dipole moment
!
      USE constants
      USE resolution
      IMPLICIT NONE

      DOUBLE COMPLEX, INTENT(IN)     :: xsi(l*(l+3)/2,2*q)
      DOUBLE PRECISION, INTENT(IN)   :: lambda(l,q)
      DOUBLE PRECISION, INTENT(IN)   :: normal(l,q)
      DOUBLE PRECISION, INTENT(OUT)  :: lx,ly,lz
      DOUBLE PRECISION               :: norm,lam1
      DOUBLE PRECISION               :: ampl,stwo
      DOUBLE PRECISION               :: sj,sjp
      INTEGER, INTENT(IN)            :: opt
      INTEGER                        :: q1

      lx = 0.d0
      ly = 0.d0
      lz = 0.d0
      stwo = SQRT(2.d0)
      DO q1 = 1,2*q                     ! sum over q
         IF (q1.le.q) THEN
            lam1 = lambda(1,q1)
            norm = normal(1,q1)
         ELSE
            lam1 = -lambda(1,2*q-q1+1)
            norm = normal(1,2*q-q1+1)
         ENDIF
         CALL sphbes(1,ABS(lam1),sj,sjp)
         IF (opt.eq.0) THEN
            ampl = 4*SQRT(PI/3)*norm*sjp*ABS(lam1)/lam1
         ELSE
            ampl = 4*SQRT(PI/3)*norm*sjp*ABS(lam1)
         ENDIF
         lz = lz-ampl*xsi(1,q1)
         lx = lx+ampl*stwo*DBLE(xsi(2,q1))
         ly = ly-ampl*stwo*AIMAG(xsi(2,q1))
      ENDDO

      RETURN
      END SUBROUTINE momenaux

!*****************************************************************
      SUBROUTINE spectruml(lambda,xsi,ext,opt)
!*****************************************************************
!
! Computes the energy and helicity spectra in shells 
! of constant l. The spectra are written to the files 
! '*spectruml.ext.txt' and '*helicityl.ext.txt', where 
! * is 'k' or 'm' according to the value of the input 
! variable opt.
!
! Parameters
!    lambda: table with the roots of the Besselj function
!    xsi   : expanding coefficients for the field
!    ext   : extension used for the file
!    opt   : =0 kinetic spectrum
!            =1 magnetic spectrum
!
      USE resolution
      IMPLICIT NONE

      DOUBLE COMPLEX, INTENT(IN)     :: xsi(l*(l+3)/2,2*q)
      DOUBLE PRECISION, INTENT(IN)   :: lambda(l,q)
      DOUBLE PRECISION               :: enel(l),hell(l)
      DOUBLE PRECISION               :: lam1,tmp
      INTEGER, INTENT(IN)            :: opt
      INTEGER                        :: q1,l1,m1
      CHARACTER(len=3), INTENT(IN)   :: ext

      DO l1 = 1,l                       ! sets the spectra to zero
         enel(l1) = 0.d0
         hell(l1) = 0.d0
      END DO
      IF (opt.eq.0) THEN                ! computes the kinetic spectra
         DO q1 = 1,2*q
         DO l1 = 1,l
            IF (q1.le.q) THEN
               lam1 = lambda(l1,q1)
            ELSE
               lam1 = -lambda(l1,2*q-q1+1)
            ENDIF
            tmp = ABS(xsi(l1*(l1+1)/2,q1))**2
            enel(l1) = enel(l1)+.5d0*tmp
            hell(l1) = hell(l1)+.5d0*tmp*lam1
            DO m1 = 1,l1
               tmp = ABS(xsi(l1*(l1+1)/2+m1,q1))**2
               enel(l1) = enel(l1)+tmp
               hell(l1) = hell(l1)+tmp*lam1
            ENDDO
         ENDDO
         ENDDO
      ELSE                              ! computes the magnetic spectra
         DO q1 = 1,2*q
         DO l1 = 1,l
            IF (q1.le.q) THEN
               lam1 = 1.d0/lambda(l1,q1)
            ELSE
               lam1 = -1.d0/lambda(l1,2*q-q1+1)
            ENDIF
            tmp = ABS(xsi(l1*(l1+1)/2,q1))**2
            enel(l1) = enel(l1)+.5d0*tmp
            hell(l1) = hell(l1)+.5d0*tmp*lam1
            DO m1 = 1,l1
               tmp = ABS(xsi(l1*(l1+1)/2+m1,q1))**2
               enel(l1) = enel(l1)+tmp
               hell(l1) = hell(l1)+tmp*lam1
            ENDDO
         ENDDO
         ENDDO
      ENDIF

      IF (opt.eq.0) THEN                ! exports the energy
         OPEN(1,file='kspectruml.' // ext // '.txt')
      ELSE
         OPEN(1,file='mspectruml.' // ext // '.txt')
      ENDIF
      DO l1 = 1,l
         WRITE(1,30) l1,enel(l1)
   30    FORMAT( I4,E22.14 )
      ENDDO
      CLOSE(1)
      IF (opt.eq.0) THEN                ! exports the helicity
         OPEN(1,file='khelicityl.' // ext // '.txt')
      ELSE
         OPEN(1,file='mhelicityl.' // ext // '.txt')
      ENDIF
      DO l1 = 1,l
         WRITE(1,40) l1,hell(l1)
   40    FORMAT( I4,E22.14 )
      ENDDO
      CLOSE(1)

      RETURN
      END SUBROUTINE spectruml

!*****************************************************************
      SUBROUTINE spectrumq(lambda,xsi,ext,opt)
!*****************************************************************
!
! Computes the energy and helicity spectra in shells 
! of lambda. The number of bins in lambda are given by 
! q, and the width is max(lambda)/q. The spectra are 
! written to the files '*spectrumq.ext.txt' and 
! '*helicityq.ext.txt', where * is 'k' or 'm' according 
! to the value of the input variable opt.
!
! Parameters
!    lambda: table with the roots of the Besselj function
!    xsi   : expanding coefficients for the field
!    ext   : extension used for the file
!    opt   : =0 kinetic spectrum
!            =1 magnetic spectrum
!
      USE resolution
      IMPLICIT NONE

      DOUBLE COMPLEX, INTENT(IN)     :: xsi(l*(l+3)/2,2*q)
      DOUBLE PRECISION, INTENT(IN)   :: lambda(l,q)
      DOUBLE PRECISION               :: eneq(q),helq(q),lamq(q)
      DOUBLE PRECISION               :: lam1,tmp,bin
      INTEGER, INTENT(IN)            :: opt
      INTEGER                        :: q1,q2,l1,m1
      CHARACTER(len=3), INTENT(IN)   :: ext

      bin = (MAXVAL(lambda))/q
      DO q1 = 1,q                       ! sets the spectra to zero
         eneq(q1) = 0.d0
         helq(q1) = 0.d0
         lamq(q1) = q1*bin              ! generates the bins
      END DO
      IF (opt.eq.0) THEN                ! kinetic spectra
         DO q1 = 1,q                    ! computes the spectra (q>0)
         DO l1 = 1,l
            q2 = INT(lambda(l1,q1)/bin+.99)
            lam1 = lambda(l1,q1)
            tmp = ABS(xsi(l1*(l1+1)/2,q1))**2
            eneq(q2) = eneq(q2)+.5d0*tmp
            helq(q2) = helq(q2)+.5d0*tmp*lam1
            DO m1 = 1,l1
               tmp = ABS(xsi(l1*(l1+1)/2+m1,q1))**2
               eneq(q2) = eneq(q2)+tmp
               helq(q2) = helq(q2)+tmp*lam1
            ENDDO
         ENDDO
         ENDDO
         DO q1 = q+1,2*q                ! computes the spectra (q<0)
         DO l1 = 1,l
            q2 = INT(ABS(lambda(l1,2*q-q1+1))/bin+.99)
            lam1 = -lambda(l1,2*q-q1+1)
            tmp = ABS(xsi(l1*(l1+1)/2,q1))**2
            eneq(q2) = eneq(q2)+.5d0*tmp
            helq(q2) = helq(q2)+.5d0*tmp*lam1
            DO m1 = 1,l1
               tmp = ABS(xsi(l1*(l1+1)/2+m1,q1))**2
               eneq(q2) = eneq(q2)+tmp
               helq(q2) = helq(q2)+tmp*lam1
            ENDDO
         ENDDO
         ENDDO
      ELSE                              ! magnetic spectra
         DO q1 = 1,q                    ! computes the spectra (q>0)
         DO l1 = 1,l
            q2 = INT(lambda(l1,q1)/bin+.99)
            lam1 = 1.d0/lambda(l1,q1)
            tmp = ABS(xsi(l1*(l1+1)/2,q1))**2
            eneq(q2) = eneq(q2)+.5d0*tmp
            helq(q2) = helq(q2)+.5d0*tmp*lam1
            DO m1 = 1,l1
               tmp = ABS(xsi(l1*(l1+1)/2+m1,q1))**2
               eneq(q2) = eneq(q2)+tmp
               helq(q2) = helq(q2)+tmp*lam1
            ENDDO
         ENDDO
         ENDDO
         DO q1 = q+1,2*q                ! computes the spectra (q<0)
         DO l1 = 1,l
            q2 = INT(ABS(lambda(l1,2*q-q1+1))/bin+.99)
            lam1 = -1.d0/lambda(l1,2*q-q1+1)
            tmp = ABS(xsi(l1*(l1+1)/2,q1))**2
            eneq(q2) = eneq(q2)+.5d0*tmp
            helq(q2) = helq(q2)+.5d0*tmp*lam1
            DO m1 = 1,l1
               tmp = ABS(xsi(l1*(l1+1)/2+m1,q1))**2
               eneq(q2) = eneq(q2)+tmp
               helq(q2) = helq(q2)+tmp*lam1
            ENDDO
         ENDDO
         ENDDO
      ENDIF

      IF (opt.eq.0) THEN                ! exports the energy
         OPEN(1,file='kspectrumq.' // ext // '.txt')
      ELSE
         OPEN(1,file='mspectrumq.' // ext // '.txt')
      ENDIF
      DO q1 = 1,q
         WRITE(1,50) lamq(q1),eneq(q1)
   50    FORMAT( E22.14,E22.14 )
      ENDDO
      CLOSE(1)
      IF (opt.eq.0) THEN                ! exports the helicity
         OPEN(1,file='khelicityq.' // ext // '.txt')
      ELSE
         OPEN(1,file='mhelicityq.' // ext // '.txt')
      ENDIF
      DO q1 = 1,q
         WRITE(1,60) lamq(q1),helq(q1)
   60    FORMAT( E22.14,E22.14 )
      ENDDO
      CLOSE(1)

      RETURN
      END SUBROUTINE spectrumq

!*****************************************************************
      SUBROUTINE genext(ind,ext)
!*****************************************************************
!
! Converts the index ind into a character.
!
! Parameters
!    ind: index
!    ext: character
!
      IMPLICIT NONE

      INTEGER, INTENT(IN)           :: ind
      INTEGER                       :: ic,id,iu
      CHARACTER(len=3), INTENT(OUT) :: ext

      ic = 48+int(ind/100)
      id = 48+int(ind/10)-int(ind/100)*10
      iu = 48+int(ind)-int(ind/10)*10
      ext = char(ic) // char(id) // char(iu)

      RETURN
      END SUBROUTINE genext
