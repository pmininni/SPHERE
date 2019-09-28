!=================================================================
! TABLES subroutines
!
! Subroutines to compute tables for the SPHERE code.
!
! 2005 Pablo D. Mininni.
!      National Center for Atmospheric Research.
!      e-mail: mininni@ucar.edu 
!=================================================================

!*****************************************************************
       SUBROUTINE lambda_table(lambda)
!-----------------------------------------------------------------
!
! Creates a table lambda(l,q) with the roots of the 
! spherical Besselj functions for all values of l and q; 
! l is the order of the Besselj function, and q labels the 
! zero in sequential order.
!
! Parameters
!    lambda: at the output contains the table with the roots
!
      USE constants
      USE resolution
      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(OUT)  :: lambda(l,q)
      DOUBLE PRECISION               :: xini,xend,root
      INTEGER                        :: i,j
!
! q+1 roots of spherical Bessel j_1(x): (l,q) = (1,q)
!
      DO j = 1,q
         xini = PI*j
         xend = PI*(j+1)
         CALL zerobes(1,xini,xend,root)
         lambda(1,j) = root
      END DO
!
! q roots of spherical Bessel j_l(x): 2<=l<=q
!
      DO i = 2,l
         DO j = 1,q-1
            xini = lambda(i-1,j)
            xend = lambda(i-1,j+1)
            CALL zerobes(i,xini,xend,root)
            lambda(i,j) = root
         END DO
         xini = lambda(i-1,q)
         xend = 2*lambda(i-1,q)-lambda(i-1,q-1)
         CALL zerobes(i,xini,xend,root)
         lambda(i,q) = root
      END DO

      RETURN
      END SUBROUTINE lambda_table

!*****************************************************************
      SUBROUTINE normal_table(lambda,normal)
!-----------------------------------------------------------------
!
! Creates a table normal(l,q) with the normalization 
! coefficients of the Chandrasekhar-Kendall functions in 
! the sphere, as a function of l and q.
!
! Parameter
!    lambda: table with the roots of the Besselj function
!    normal: at the output contains the normalization factors
!
      USE resolution
      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(IN)   :: lambda(l,q)
      DOUBLE PRECISION, INTENT(OUT)  :: normal(l,q)
      DOUBLE PRECISION               :: tmp,tmq
      DOUBLE PRECISION               :: sj,sjp
      INTEGER                        :: i,j

      DO j = 1,q
         DO i = 1,l
            tmp = SQRT(i*(i+1.d0))
            CALL sphbes(i+1,lambda(i,j),sj,sjp)
            tmq = tmp*lambda(i,j)*ABS(sj)
            normal(i,j) = 1.d0/tmq
         END DO
      END DO

      RETURN
      END SUBROUTINE normal_table

!*****************************************************************
      SUBROUTINE cicoef_table(lambda,normal,c)
!-----------------------------------------------------------------
!
! Builds a table with coupling coefficients for the momentum 
! and induction equation C(1,2,3) = int[J(1).J(2)xJ(3)], where 
! J are the Chandrasekhar-Kendall functions. The coefficients 
! are given as a function of the indices (l3,ind2,q3,q2,ind1,q1), 
! where ind1 and ind2 are angular indices for the pair of values 
! (l1,m1) and (l2,m2) respectively. The indices q1, q2, and q3 
! run from 1 to 2*q, and the actual order of q is e.g. for q1, 
! q1 = 1,2,...,q,-q,-q+1,...,-1, where negative values of q1 
! denote negative values of the Chandrasekhar-Kendall eigenvalue 
! of the curl operator lambda(l,q). The angular indices are 
! defined as ind1=l1(l1+1)/2+m1 (0<=m1<=l1) and ind2=l2(l2+1)+m2 
! (-l2<=m2<=l2). Note that only the coefficients for positive 
! values of m1 are computed, since the evolution of the modes 
! with negative m follows from the relation 
! xsi(q,l,-m)=(-1)^m*CONJG(xsi(q,l,m)). The value of m3 is given 
! by the orthogonality condition m3=m1-m2.
!
! Parameter
!    lambda: table with the roots of the Besselj function
!    normal: table with the normalization factors
!    c     : at the output contains the coupling coefficients
!
      USE mpivars
      USE constants
      USE functions
      USE resolution
      IMPLICIT NONE

      DOUBLE COMPLEX, INTENT(OUT)  :: c(l,l*(l+2),2*q,2*q,ista:iend,qsta:qend)
      DOUBLE PRECISION, INTENT(IN) :: lambda(l,q)
      DOUBLE PRECISION, INTENT(IN) :: normal(l,q)
      DOUBLE PRECISION             :: xrad(nrad),wrad(nrad)
      DOUBLE PRECISION             :: xang(nang),wang(nang)
      DOUBLE PRECISION             :: lam1,lam2,lam3,norm
      DOUBLE PRECISION             :: R,RD1,RD2,RD3
      DOUBLE PRECISION             :: RDD1,RDD2,RDD3
      DOUBLE PRECISION             :: X,XD1,XD2,XD3
      DOUBLE PRECISION             :: XDD1,XDD2,XDD3
      TYPE(ARGUMENTS)              :: arg123,arg231,arg312
      INTEGER                      :: q1,q2,q3
      INTEGER                      :: l1,l2,l3
      INTEGER                      :: m1,m2,m3
      INTEGER                      :: ind1,ind2
      INTEGER                      :: odd1,odd2
      INTEGER                      :: i
!
! Computes weights and abscissas for quadratures
!
      CALL gauleg(0.d0,1.d0,xrad,wrad,nrad)
      CALL gauleg(-1.d0,1.d0,xang,wang,nang)
!
! Main loop: computes the coupling coefficients
!
 QL : DO q1 = qsta,qend
 IL : DO ind1 = ista,iend
      l1 = INT(.5*(SQRT(1+8*FLOAT(ind1))-1))
      IF (q1.le.q) THEN
         lam1 = lambda(l1,q1)
         norm = normal(l1,q1)
      ELSE
         lam1 = -lambda(l1,2*q-q1+1)
         norm = normal(l1,2*q-q1+1)
      ENDIF

 QP:  DO q2 = 1,q
      DO q3 = 1,q

 LP:  DO l2 = 1,l
      DO l3 = 1,l
!
! Computes the radial integrals
!
      R = 0
      RD1 = 0 ; RD2 = 0 ; RD3 = 0
      RDD1 = 0 ; RDD2 = 0 ; RDD3 = 0
      arg123%l1 = l1
      arg123%l2 = l2
      arg123%l3 = l3
      arg123%lam1 = ABS(lam1)
      arg123%lam2 = lambda(l2,q2)
      arg123%lam3 = lambda(l3,q3)
      arg231%l1 = l2
      arg231%l2 = l3
      arg231%l3 = l1
      arg231%lam1 = lambda(l2,q2)
      arg231%lam2 = lambda(l3,q3)
      arg231%lam3 = ABS(lam1)
      arg312%l1 = l3
      arg312%l2 = l1
      arg312%l3 = l2
      arg312%lam1 = lambda(l3,q3)
      arg312%lam2 = ABS(lam1)
      arg312%lam3 = lambda(l2,q2)
 IR : DO i = 1,nrad
         R = R+wrad(i)*TRIPLER(xrad(i),arg123)
         RD1 = RD1+wrad(i)*TRIPLERD(xrad(i),arg231)
         RD2 = RD2+wrad(i)*TRIPLERD(xrad(i),arg312)
         RD3 = RD3+wrad(i)*TRIPLERD(xrad(i),arg123)
         RDD1 = RDD1+wrad(i)*TRIPLERDD(xrad(i),arg123)
         RDD2 = RDD2+wrad(i)*TRIPLERDD(xrad(i),arg231)
         RDD3 = RDD3+wrad(i)*TRIPLERDD(xrad(i),arg312)
      ENDDO IR

      m1 = ind1-l1*(l1+1)/2          ! only positive values of m1
 MP : DO m2 = -l2,l2
      m3 = m1-m2                     ! orthogonality in m

      ind2 = l2*(l2+1)+m2            ! index ind2
      c(l3,ind2,q3,q2,ind1,q1) = 0   ! fill the array with zeros

 OM : IF (ABS(m3).le.l3) THEN        ! -l<=m<=l
!
! Computes the angular integrals
!
      X = 0
      XD1 = 0 ; XD2 = 0 ; XD3 = 0
      XDD1 = 0 ; XDD2 = 0 ; XDD3 = 0
      arg123%m1 = ABS(m1)
      arg123%m2 = ABS(m2)
      arg123%m3 = ABS(m3)
      arg231%m1 = ABS(m2)
      arg231%m2 = ABS(m3)
      arg231%m3 = ABS(m1)
      arg312%m1 = ABS(m3)
      arg312%m2 = ABS(m1)
      arg312%m3 = ABS(m2)
      odd1 = l1+l2+l3+m1+m2+m3
      odd2 = odd1/2
      odd2 = 2*odd2
      IF (odd1.eq.odd2) THEN         ! l1+l2+l3+m1+m2+m3 even
 XE :    DO i = 1,nang
            X = X+wang(i)*TRIPLEX(xang(i),arg123)
            XDD1 = XDD1+wang(i)*TRIPLEXDD(xang(i),arg123)
            XDD2 = XDD2+wang(i)*TRIPLEXDD(xang(i),arg231)
            XDD3 = XDD3+wang(i)*TRIPLEXDD(xang(i),arg312)
         ENDDO XE
      ELSE                           ! l1+l2+l3+m1+m2+m3 odd
 XO :    DO i = 1,nang
            XD1 = XD1+wang(i)*TRIPLEXD(xang(i),arg231)
            XD2 = XD2+wang(i)*TRIPLEXD(xang(i),arg312)
            XD3 = XD3+wang(i)*TRIPLEXD(xang(i),arg123)
         ENDDO XO
      ENDIF
      IF (m2.lt.0) THEN
         X = X*(-1)**ABS(m2)
         XD1 = XD1*(-1)**ABS(m2)
         XD2 = XD2*(-1)**ABS(m2)
         XD3 = XD3*(-1)**ABS(m2)
         XDD1 = XDD1*(-1)**ABS(m2)
         XDD2 = XDD2*(-1)**ABS(m2)
         XDD3 = XDD3*(-1)**ABS(m2)
      ENDIF
      IF (m3.lt.0) THEN
         X = X*(-1)**ABS(m3)
         XD1 = XD1*(-1)**ABS(m3)
         XD2 = XD2*(-1)**ABS(m3)
         XD3 = XD3*(-1)**ABS(m3)
         XDD1 = XDD1*(-1)**ABS(m3)
         XDD2 = XDD2*(-1)**ABS(m3)
         XDD3 = XDD3*(-1)**ABS(m3)
      ENDIF
!
! Coupling coefficients
!
      lam2 = lambda(l2,q2)           ! C(1,2,3) for q2,q3>0
      lam3 = lambda(l3,q3)
      c(l3,ind2,q3,q2,ind1,q1) = &
          2*PI*norm*normal(l2,q2)*normal(l3,q3)* &
          (l1*(l1+1)*((lam3*RD2-lam2*RD3)*(m2*m3*X-XDD1)+ &
          IM*(m3*XD2-m2*XD3)*(lam2*lam3*R+RDD1))+ &
          l2*(l2+1)*((lam3*RD1-lam1*RD3)*(m1*m3*X+XDD2)- &
          IM*(m3*XD1+m1*XD3)*(lam1*lam3*R+RDD2))+ &
          l3*(l3+1)*((lam1*RD2-lam2*RD1)*(m1*m2*X+XDD3)+ &
          IM*(m2*XD1+m1*XD2)*(lam1*lam2*R+RDD3)))
      lam2 = -lambda(l2,q2)         ! C(1,2,3) for q2<0
      c(l3,ind2,q3,2*q-q2+1,ind1,q1) = &
          2*PI*norm*normal(l2,q2)*normal(l3,q3)* &
          (l1*(l1+1)*((lam3*RD2-lam2*RD3)*(m2*m3*X-XDD1)+ &
          IM*(m3*XD2-m2*XD3)*(lam2*lam3*R+RDD1))+ &
          l2*(l2+1)*((lam3*RD1-lam1*RD3)*(m1*m3*X+XDD2)- &
          IM*(m3*XD1+m1*XD3)*(lam1*lam3*R+RDD2))+ &
          l3*(l3+1)*((lam1*RD2-lam2*RD1)*(m1*m2*X+XDD3)+ &
          IM*(m2*XD1+m1*XD2)*(lam1*lam2*R+RDD3)))
      lam2 = lambda(l2,q2)          ! C(1,2,3) for q3<0
      lam3 = -lambda(l3,q3)
      c(l3,ind2,2*q-q3+1,q2,ind1,q1) = &
          2*PI*norm*normal(l2,q2)*normal(l3,q3)* &
          (l1*(l1+1)*((lam3*RD2-lam2*RD3)*(m2*m3*X-XDD1)+ &
          IM*(m3*XD2-m2*XD3)*(lam2*lam3*R+RDD1))+ &
          l2*(l2+1)*((lam3*RD1-lam1*RD3)*(m1*m3*X+XDD2)- &
          IM*(m3*XD1+m1*XD3)*(lam1*lam3*R+RDD2))+ &
          l3*(l3+1)*((lam1*RD2-lam2*RD1)*(m1*m2*X+XDD3)+ &
          IM*(m2*XD1+m1*XD2)*(lam1*lam2*R+RDD3)))
      lam2 = -lambda(l2,q2)         ! C(1,2,3) for q2,q3<0
      c(l3,ind2,2*q-q3+1,2*q-q2+1,ind1,q1) = &
          2*PI*norm*normal(l2,q2)*normal(l3,q3)* &
          (l1*(l1+1)*((lam3*RD2-lam2*RD3)*(m2*m3*X-XDD1)+ &
          IM*(m3*XD2-m2*XD3)*(lam2*lam3*R+RDD1))+ &
          l2*(l2+1)*((lam3*RD1-lam1*RD3)*(m1*m3*X+XDD2)- &
          IM*(m3*XD1+m1*XD3)*(lam1*lam3*R+RDD2))+ &
          l3*(l3+1)*((lam1*RD2-lam2*RD1)*(m1*m2*X+XDD3)+ &
          IM*(m2*XD1+m1*XD2)*(lam1*lam2*R+RDD3)))

      ENDIF OM

      ENDDO MP

      ENDDO
      ENDDO LP

      ENDDO
      ENDDO QP

      ENDDO IL
      ENDDO QL

      RETURN
      END SUBROUTINE cicoef_table

!*****************************************************************
      SUBROUTINE coriol_table(lambda,normal,ixp,ixm,iz)
!-----------------------------------------------------------------
!
! Builds a table with coupling coefficients for the Coriolis 
! term in the momentum equation. The coefficients are given as 
! a function of the indices (l2,q2,ind1,q1) where ind1 is the 
! angular index for the pair of values (l1,m1). The indices q1 
! and q2 run from 1 to 2*q, and the actual order of q is e.g. 
! for q1, q1 = 1,2,...,q,-q,-q+1,...,-1, where negative values 
! of q1 denote negative values of the Chandrasekhar-Kendall 
! eigenvalue of the curl operator lambda(l,q). The angular 
! index is defined as ind1=l1(l1+1)/2+m1 (0<=m1<=l1). Note that 
! only the coefficients for positive values of m1 are computed. 
! The y-components of the integrals are given by Iyp = i*Ixp 
! and Iym = -i*Iym. The value of m2 is given by the 
! orthogonality conditions. 
!
! Parameter
!    lambda: table with the roots of the Besselj function
!    normal: table with the normalization factors
!    ixp   : coupling coefficients for (w_x+i w_y) and m2=m1+1
!    ixm   : coupling coefficients for (w_x-i w_y) and m2=m1-1
!    iz    : coupling coefficients for w_z and m2=m1
!
      USE mpivars
      USE constants
      USE functions
      USE resolution
      IMPLICIT NONE

      DOUBLE COMPLEX, INTENT(OUT)  :: ixp(l,2*q,ista:iend,qsta:qend)
      DOUBLE COMPLEX, INTENT(OUT)  :: ixm(l,2*q,ista:iend,qsta:qend)
      DOUBLE COMPLEX, INTENT(OUT)  :: iz(l,2*q,ista:iend,qsta:qend)
      DOUBLE PRECISION, INTENT(IN) :: lambda(l,q)
      DOUBLE PRECISION, INTENT(IN) :: normal(l,q)
      DOUBLE PRECISION             :: xrad(nrad),wrad(nrad)
      DOUBLE PRECISION             :: xang(nang),wang(nang)
      DOUBLE PRECISION             :: lam1,lam2,norm
      DOUBLE PRECISION             :: R,RD1,RDD
      DOUBLE PRECISION             :: QI,QD1,QD2
      DOUBLE PRECISION             :: W,WD1,WD2,WDD
      DOUBLE PRECISION             :: X,XD1,XD2,XDD
      DOUBLE PRECISION             :: Y,YD1,YD2,YDD
      DOUBLE PRECISION             :: Z,ZD1,ZD2,ZDD
      TYPE(ARGUMENTS)              :: arg12,arg21
      INTEGER                      :: q1,q2
      INTEGER                      :: l1,l2
      INTEGER                      :: m1,m2
      INTEGER                      :: ind1,ind2
      INTEGER                      :: odd1,odd2
      INTEGER                      :: i
!
! Computes weights and abscissas for quadratures
!
      CALL gauleg(0.d0,1.d0,xrad,wrad,nrad)
      CALL gauleg(-1.d0,1.d0,xang,wang,nang)
!
! Main loop: computes the coupling coefficients
!
 QL : DO q1 = qsta,qend
 IL : DO ind1 = ista,iend
      l1 = INT(.5*(SQRT(1+8*FLOAT(ind1))-1))
      IF (q1.le.q) THEN
         lam1 = lambda(l1,q1)
         norm = normal(l1,q1)
      ELSE
         lam1 = -lambda(l1,2*q-q1+1)
         norm = normal(l1,2*q-q1+1)
      ENDIF

 QP:  DO q2 = 1,q
 LP:  DO l2 = 1,l
!
! Computes the radial integrals
!
      R = 0 ; RD1 = 0 ; RDD = 0
      QI = 0 ; QD1 = 0 ; QD2 = 0
      arg12%l1 = l1
      arg12%l2 = l2
      arg12%l3 = 0
      arg12%lam1 = ABS(lam1)
      arg12%lam2 = lambda(l2,q2)
      arg12%lam3 = 0.d0
      arg21%l1 = l2
      arg21%l2 = l1
      arg21%l3 = 0
      arg21%lam1 = lambda(l2,q2)
      arg21%lam2 = ABS(lam1)
      arg21%lam3 = 0.d0
 IR : DO i = 1,nrad
         R = R+wrad(i)*DOUBLER(xrad(i),arg12)
         RD1 = RD1+wrad(i)*DOUBLERD(xrad(i),arg21)
         RDD = RDD+wrad(i)*DOUBLERDD(xrad(i),arg12)
         QI = QI+wrad(i)*DOUBLEQ(xrad(i),arg12)
         QD1 = QD1+wrad(i)*DOUBLEQD(xrad(i),arg21)
         QD2 = QD2+wrad(i)*DOUBLEQD(xrad(i),arg12)
      ENDDO IR

      m1 = ind1-l1*(l1+1)/2          ! only positive values of m1
      m2 = m1                        ! orthogonality in m for Iz
      iz(l2,q2,ind1,q1) = 0          ! fill the array with zeros

 LZ : IF (ABS(m2).le.l2) THEN        ! -l<=m<=l
!
! Computes the angular integrals for Iz
!
      X = 0 ; XD1 = 0 ; XD2 = 0 ; XDD = 0
      Z = 0 ; ZD1 = 0 ; ZD2 = 0
      arg12%m1 = ABS(m1)
      arg12%m2 = ABS(m2)
      arg12%m3 = 0
      arg21%m1 = ABS(m2)
      arg21%m2 = ABS(m1)
      arg21%m3 = 0
      odd1 = l1+l2+m1+m2
      odd2 = odd1/2
      odd2 = 2*odd2
      IF (odd1.eq.odd2) THEN         ! l1+l2+m1+m2 even
 ZE :    DO i = 1,nang
            XD1 = XD1+wang(i)*DOUBLEXD(xang(i),arg21)
            XD2 = XD2+wang(i)*DOUBLEXD(xang(i),arg12)
            Z = Z+wang(i)*DOUBLEZ(xang(i),arg12)
         ENDDO ZE
      ELSE                           ! l1+l2+m1+m2 odd
 ZO :    DO i = 1,nang
            X = X+wang(i)*DOUBLEX(xang(i),arg12)
            XDD = XDD+wang(i)*DOUBLEXDD(xang(i),arg12)
            ZD1 = ZD1+wang(i)*DOUBLEZD(xang(i),arg21)
            ZD2 = ZD2+wang(i)*DOUBLEZD(xang(i),arg12)
         ENDDO ZO
      ENDIF
      IF (m2.lt.0) THEN
         X = X*(-1)**ABS(m2)
         XD1 = XD1*(-1)**ABS(m2)
         XD2 = XD2*(-1)**ABS(m2)
         XDD = XDD*(-1)**ABS(m2)
         Z = Z*(-1)**ABS(m2)
         ZD1 = ZD1*(-1)**ABS(m2)
         ZD2 = ZD2*(-1)**ABS(m2)
      ENDIF
!
! Iz coupling coefficients
!
      lam2 = lambda(l2,q2)           ! Iz(1,2) for q2>0
      iz(l2,q2,ind1,q1) = 4*PI*norm*normal(l2,q2)* &
          (IM*(lam1*lam2*R+RDD)*(m1*XD2+m2*XD1)- &
          RD1*(lam1+lam2)*(m1*m2*X+XDD)+ &
          IM*Z*(m2*l1*(l1+1)*QD2+m1*l2*(l2+1)*QD1)- &
          QI*(lam2*l1*(l1+1)*ZD2-lam1*l2*(l2+1)*ZD1))
      lam2 = -lambda(l2,q2)          ! Iz(1,2) for q2<0
      iz(l2,2*q-q2+1,ind1,q1) = 4*PI*norm*normal(l2,q2)* &
          (IM*(lam1*lam2*R+RDD)*(m1*XD2+m2*XD1)- &
          RD1*(lam1+lam2)*(m1*m2*X+XDD)+ &
          IM*Z*(m2*l1*(l1+1)*QD2+m1*l2*(l2+1)*QD1)- &
          QI*(lam2*l1*(l1+1)*ZD2-lam1*l2*(l2+1)*ZD1))

      ENDIF LZ

      m2 = m1+1                      ! orthogonality in m for Ixp
      ixp(l2,q2,ind1,q1) = 0         ! fill the array with zeros

 IP : IF (ABS(m2).le.l2) THEN        ! -l<=m<=l
!
! Computes the angular integrals for Ixp
!
      W = 0 ; WD1 = 0 ; WD2 = 0 ; WDD = 0
      Y = 0 ; YD1 = 0 ; YD2 = 0
      arg12%m1 = ABS(m1)
      arg12%m2 = ABS(m2)
      arg12%m3 = 0
      arg21%m1 = ABS(m2)
      arg21%m2 = ABS(m1)
      arg21%m3 = 0
      odd1 = l1+l2+m1+m2
      odd2 = odd1/2
      odd2 = 2*odd2
      IF (odd1.eq.odd2) THEN         ! l1+l2+m1+m2 even
 YE :    DO i = 1,nang
            W = W+wang(i)*DOUBLEW(xang(i),arg12)
            WDD = WDD+wang(i)*DOUBLEWDD(xang(i),arg12)
            YD1 = YD1+wang(i)*DOUBLEYD(xang(i),arg21)
            YD2 = YD2+wang(i)*DOUBLEYD(xang(i),arg12)
         ENDDO YE
      ELSE                           ! l1+l2+m1+m2 odd
 YO :    DO i = 1,nang
            WD1 = WD1+wang(i)*DOUBLEWD(xang(i),arg21)
            WD2 = WD2+wang(i)*DOUBLEWD(xang(i),arg12)
            Y = Y+wang(i)*DOUBLEY(xang(i),arg12)
         ENDDO YO
      ENDIF
      IF (m2.lt.0) THEN
         W = W*(-1)**ABS(m2)
         WD1 = WD1*(-1)**ABS(m2)
         WD2 = WD2*(-1)**ABS(m2)
         WDD = WDD*(-1)**ABS(m2)
         Y = Y*(-1)**ABS(m2)
         YD1 = YD1*(-1)**ABS(m2)
         YD2 = YD2*(-1)**ABS(m2)
      ENDIF
!
! Ixp coupling coefficients
!
      lam2 = lambda(l2,q2)           ! Ixp(1,2) for q2>0
      ixp(l2,q2,ind1,q1) = 2*PI*norm*normal(l2,q2)* &
          (IM*(lam1*lam2*R+RDD)*(m1*WD2+m2*WD1)- &
          RD1*(lam1+lam2)*(m1*m2*W+WDD)- &
          IM*Y*(m2*l1*(l1+1)*QD2+m1*l2*(l2+1)*QD1)+ &
          QI*(lam2*l1*(l1+1)*YD2-lam1*l2*(l2+1)*YD1)+ &
          QI*W*(m2*lam2*l1*(l1+1)+m1*lam1*l2*(l2+1))- &
          IM*(l1*(l1+1)*QD2*WD2-l2*(l2+1)*QD1*WD1))
      lam2 = -lambda(l2,q2)          ! Ixp(1,2) for q2<0
      ixp(l2,2*q-q2+1,ind1,q1) = 2*PI*norm*normal(l2,q2)* &
          (IM*(lam1*lam2*R+RDD)*(m1*WD2+m2*WD1)- &
          RD1*(lam1+lam2)*(m1*m2*W+WDD)- &
          IM*Y*(m2*l1*(l1+1)*QD2+m1*l2*(l2+1)*QD1)+ &
          QI*(lam2*l1*(l1+1)*YD2-lam1*l2*(l2+1)*YD1)+ &
          QI*W*(m2*lam2*l1*(l1+1)+m1*lam1*l2*(l2+1))- &
          IM*(l1*(l1+1)*QD2*WD2-l2*(l2+1)*QD1*WD1))
      ENDIF IP

      m2 = m1-1                      ! orthogonality in m for Ixm
      ixm(l2,q2,ind1,q1) = 0         ! fill the array with zeros

 LM : IF (ABS(m2).le.l2) THEN        ! -l<=m<=l
!
! Computes the angular integrals for Ixm
!
      W = 0 ; WD1 = 0 ; WD2 = 0 ; WDD = 0
      Y = 0 ; YD1 = 0 ; YD2 = 0
      arg12%m1 = m1
      arg12%m2 = ABS(m2)
      arg12%m3 = 0
      arg21%m1 = ABS(m2)
      arg21%m2 = m1
      arg21%m3 = 0
      odd1 = l1+l2+m1+m2
      odd2 = odd1/2
      odd2 = 2*odd2
      IF (odd1.eq.odd2) THEN         ! l1+l2+m1+m2 even
 XE :    DO i = 1,nang
            W = W+wang(i)*DOUBLEW(xang(i),arg12)
            WDD = WDD+wang(i)*DOUBLEWDD(xang(i),arg12)
            YD1 = YD1+wang(i)*DOUBLEYD(xang(i),arg21)
            YD2 = YD2+wang(i)*DOUBLEYD(xang(i),arg12)
         ENDDO XE
      ELSE                           ! l1+l2+m1+m2 odd
 XO :    DO i = 1,nang
            WD1 = WD1+wang(i)*DOUBLEWD(xang(i),arg21)
            WD2 = WD2+wang(i)*DOUBLEWD(xang(i),arg12)
            Y = Y+wang(i)*DOUBLEY(xang(i),arg12)
         ENDDO XO
      ENDIF
      IF (m2.lt.0) THEN
         W = W*(-1)**ABS(m2)
         WD1 = WD1*(-1)**ABS(m2)
         WD2 = WD2*(-1)**ABS(m2)
         WDD = WDD*(-1)**ABS(m2)
         Y = Y*(-1)**ABS(m2)
         YD1 = YD1*(-1)**ABS(m2)
         YD2 = YD2*(-1)**ABS(m2)
      ENDIF
!
! Ixm coupling coefficients
!
      lam2 = lambda(l2,q2)           ! Ixm(1,2) for q2>0
      ixm(l2,q2,ind1,q1) = 2*PI*norm*normal(l2,q2)* &
          (IM*(lam1*lam2*R+RDD)*(m1*WD2+m2*WD1)- &
          RD1*(lam1+lam2)*(m1*m2*W+WDD)- &
          IM*Y*(m2*l1*(l1+1)*QD2+m1*l2*(l2+1)*QD1)+ &
          QI*(lam2*l1*(l1+1)*YD2-lam1*l2*(l2+1)*YD1)- &
          QI*W*(m2*lam2*l1*(l1+1)+m1*lam1*l2*(l2+1))+ &
          IM*(l1*(l1+1)*QD2*WD2-l2*(l2+1)*QD1*WD1))
      lam2 = -lambda(l2,q2)          ! Ixm(1,2) for q2<0
      ixm(l2,2*q-q2+1,ind1,q1) = 2*PI*norm*normal(l2,q2)* &
          (IM*(lam1*lam2*R+RDD)*(m1*WD2+m2*WD1)- &
          RD1*(lam1+lam2)*(m1*m2*W+WDD)- &
          IM*Y*(m2*l1*(l1+1)*QD2+m1*l2*(l2+1)*QD1)+ &
          QI*(lam2*l1*(l1+1)*YD2-lam1*l2*(l2+1)*YD1)- &
          QI*W*(m2*lam2*l1*(l1+1)+m1*lam1*l2*(l2+1))+ &
          IM*(l1*(l1+1)*QD2*WD2-l2*(l2+1)*QD1*WD1))
      ENDIF LM

      ENDDO LP
      ENDDO QP

      ENDDO IL
      ENDDO QL

      RETURN
      END SUBROUTINE coriol_table
