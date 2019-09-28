!=================================================================
      PROGRAM SPHERE_HD
!=================================================================
! SPHERE_HD code
!
! Spectral code to evolve the hydrodynamic equations in a 
! sphere using the Chandrasekhar-Kendall eigenfunctions of 
! the curl as an orthonormal basis. Runge-Kutta of adjustable 
! order is used to evolve the system in time (the order of the 
! integrator is set in the module 'rungekutta'). The code uses 
! external tables with the coupling coefficients created by 
! CREATE_TABLES to compute the nonlinear and Coriolis terms. 
! The tolerance used to compute integrals in the coupling 
! coefficients is set in the module 'constants'.
!
! The fluid fills a sphere of radius R=1, and the boundary 
! conditions in the rotating frame are u.n=w.n=0 (no velocity 
! and vorticity perpendicular to the surface of the sphere).
!
! Since the fields are real, only the coefficients of functions 
! with positive m are evolved in time (-l<=m<=l). For negative 
! values of m the coefficients are obtained as 
! xsi(q,l,-m)=(-1)^m*CONJG(xsi(q,l,m)).
!
! INPUT : parameter.txt
!         status.txt
!         omega.txt
!         lambda.in [binary file]
!         normal.in [binary file]
!         cicoef.in [MPI native file]
!         corixp.in [MPI native file]
!         corixm.in [MPI native file]
!         corioz.in [MPI native file]
! OUTPUT: kglobal.txt
!         lmoment.txt
!         kspectruml.XXX.txt
!         kspectrumq.XXX.txt
!         khelicityl.XXX.txt
!         khelicityq.XXX.txt
!         xsiv.XXX.out [binary file]
!
! 2005 Pablo D. Mininni.
!      National Center for Atmospheric Research.
!      e-mail: mininni@ucar.edu 
!=================================================================

      USE random
      USE mpivars
      USE constants
      USE resolution
      USE rungekutta
      IMPLICIT NONE

!
! Tables and expansion coefficients

      DOUBLE COMPLEX, ALLOCATABLE, DIMENSION(:,:,:,:,:,:) :: cicoef
      DOUBLE COMPLEX, ALLOCATABLE, DIMENSION(:,:,:,:)     :: corixp
      DOUBLE COMPLEX, ALLOCATABLE, DIMENSION(:,:,:,:)     :: corixm
      DOUBLE COMPLEX, ALLOCATABLE, DIMENSION(:,:,:,:)     :: corioz
      DOUBLE COMPLEX, ALLOCATABLE, DIMENSION(:,:)         :: forv
      DOUBLE COMPLEX, ALLOCATABLE, DIMENSION(:,:)         :: xsiv
      DOUBLE COMPLEX, ALLOCATABLE, DIMENSION(:,:)         :: xsi1,xss1
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:)       :: lambda
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:)       :: normal
      INTEGER, ALLOCATABLE, DIMENSION(:)                  :: xtype

!
! Some auxiliary matrices

      DOUBLE COMPLEX     :: nonlin
      DOUBLE COMPLEX     :: coriolis
      DOUBLE COMPLEX     :: precession
      DOUBLE COMPLEX     :: xv2,xv3
      DOUBLE COMPLEX     :: angle
      DOUBLE PRECISION   :: lam1,lam3
      DOUBLE PRECISION   :: wx,wy,wz
      DOUBLE PRECISION   :: wx0,wy0,wz0
      DOUBLE PRECISION   :: wtx,wty,wtz
      DOUBLE PRECISION   :: wparam0,wparam1,wparam2
      DOUBLE PRECISION   :: wparam3,wparam4,wparam5
      DOUBLE PRECISION   :: wparam6,wparam7,wparam8
      DOUBLE PRECISION   :: wparam9
      DOUBLE PRECISION   :: injt,tmp,tmq,tmr
      DOUBLE PRECISION   :: f0,p0
      DOUBLE PRECISION   :: nu,u0
      DOUBLE PRECISION   :: phase
      DOUBLE PRECISION   :: edge
      DOUBLE PRECISION   :: dt
      DOUBLE PRECISION   :: ampl,norm
      DOUBLE PRECISION   :: sj,sjp
      REAL               :: cputime1,cputime2
      INTEGER            :: stat,angu,outs
      INTEGER            :: rand,cort,seed
      INTEGER            :: fstep,gstep
      INTEGER            :: sstep,bstep
      INTEGER            :: timef,timeg
      INTEGER            :: times,timeb
      INTEGER            :: q1,q2,q3
      INTEGER            :: l1,l2,l3
      INTEGER            :: m1,m2,m3
      INTEGER            :: ind1,ind2
      INTEGER            :: ini,step
      INTEGER            :: t,o,k
      INTEGER            :: irank
      INTEGER            :: bench
      INTEGER            :: corio
      INTEGER            :: precs
      INTEGER            :: grid
      INTEGER            :: sindex,bindex
      INTEGER            :: length
      CHARACTER(len=110) :: cname
      CHARACTER(len=100) :: tdir,odir
      CHARACTER(len=3)   :: ext

!
! Namelists for the input files
      NAMELIST / status / stat,bench,corio,angu,outs,edge,grid,precs
      NAMELIST / parameter / dt,step,gstep,sstep,bstep
      NAMELIST / parameter / f0,p0,u0,nu,rand,cort,seed,tdir,odir
      NAMELIST / omega / wx,wy,wz
      NAMELIST / omegadot / injt
      NAMELIST / omegadot / wparam0,wparam1,wparam2,wparam3,wparam4
      NAMELIST / omegadot / wparam5,wparam6,wparam7,wparam8,wparam9

!
! Initializes the MPI library

      CALL MPI_INIT(ierr)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
!
! Reads from the external file 'status.txt'
! the status of a previous run (if any)
!     stat : last output of a previous run
!     bench: performs a benchmark run
!     corio: =0 no rotation
!            =1 constant rotation in the rotating frame
!     precs: =0 no precession
!            =1 precession
!     angu : =0 skips angular momentum and dipole moment computation
!            =1 computes angular momentum and dipole moment
!     outs : =0 writes files with the expansion coefficients
!            =1 writes also files with the cartesian fields 
!     edge : edge of the box where fields are projected (d>=2)
!     grid : box resolution in each direction

      IF (myrank.eq.0) THEN
         OPEN(1,file='parameter.inp',status='unknown',form="formatted")
         READ(1,NML=status)
         CLOSE(1)
      ENDIF
      CALL MPI_BCAST(stat,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(bench,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(corio,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(angu,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(outs,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(edge,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(grid,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(precs,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

!
! Reads from the external file 'parameter.txt' the 
! parameters that will be used during the integration
!     dt   : time step size
!     step : total number of time steps to compute
!     gstep: number of timesteps between global output
!     sstep: number of timesteps between spectrum output
!     bstep: number of timesteps between binary output
!     f0   : amplitude of the external mechanic force
!     p0   : relative change in the force phase
!     u0   : amplitude of the initial velocity field
!     nu   : kinematic viscosity
!     rand : =0 constant force
!            =1 random phases
!     cort : time correlation of the external force
!     seed : seed for the random number generator
!     tdir : directory with the tables
!     odir : directory for binary output

      IF (myrank.eq.0) THEN
         OPEN(1,file='parameter.inp',status='unknown',form="formatted")
         READ(1,NML=parameter)
         CLOSE(1)
      ENDIF
      CALL MPI_BCAST(dt,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(step,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(gstep,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(sstep,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(bstep,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(f0,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(p0,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(u0,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(nu,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(rand,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(cort,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(seed,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(tdir,100,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(odir,100,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)

!
! If the simulation is done in the rotating frame,
! reads the initial value of the rigid body rotation 
!     wx: initial value for the angular velocity in x
!     wy: initial value for the angular velocity in y
!     wz: initial value for the angular velocity in z

      IF (corio.eq.1) THEN
         IF (myrank.eq.0) THEN
            OPEN(1,file='parameter.inp',status='unknown',form="formatted")
            READ(1,NML=omega)
            CLOSE(1)
         ENDIF
         CALL MPI_BCAST(wx,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
         CALL MPI_BCAST(wy,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
         CALL MPI_BCAST(wz,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      ELSE
         wx = 0.d0 ; wy = 0.d0 ; wz = 0.d0
      ENDIF
!
! If the simulation is done in the rotating frame
! with precession reads the following parameters
!     injt: stat when precession starts working     
!     wparam0-9: parameters used for generating precession

      IF (precs.eq.1) THEN
         wx0 = wx
         wy0 = wy
         wz0 = wz
         IF (myrank.eq.0) THEN
            OPEN(1,file='parameter.inp',status='unknown',form="formatted")
            READ(1,NML=omegadot)
            CLOSE(1)
         ENDIF
         CALL MPI_BCAST(injt,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
         CALL MPI_BCAST(wparam0,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
         CALL MPI_BCAST(wparam1,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
         CALL MPI_BCAST(wparam2,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
         CALL MPI_BCAST(wparam3,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
         CALL MPI_BCAST(wparam4,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
         CALL MPI_BCAST(wparam5,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
         CALL MPI_BCAST(wparam5,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
         CALL MPI_BCAST(wparam6,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
         CALL MPI_BCAST(wparam7,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
         CALL MPI_BCAST(wparam8,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
         CALL MPI_BCAST(wparam9,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      ELSE
         wtx = 0.d0 ; wty = 0.d0 ; wtz = 0.d0
      ENDIF


!
! Allocates memory

      CALL distribute(l,q,nprocs,myrank,ista,iend,qsta,qend)
      CALL create_iotype(l,q,ista,iend,qsta,qend,itype)
      ALLOCATE( xtype(0:nprocs-1) )
      CALL create_block(l,q,nprocs,xtype)
      ALLOCATE( cicoef(l,l*(l+2),2*q,2*q,ista:iend,qsta:qend) )
      IF (corio.eq.1) THEN
         CALL create_iotypc(l,q,ista,iend,qsta,qend,ctype)
         ALLOCATE( corixp(l,2*q,ista:iend,qsta:qend) )
         ALLOCATE( corixm(l,2*q,ista:iend,qsta:qend) )
         ALLOCATE( corioz(l,2*q,ista:iend,qsta:qend) )
      ENDIF
      ALLOCATE( xsi1(l*(l+3)/2,2*q), xss1(l*(l+3)/2,2*q) )
      ALLOCATE( xsiv(l*(l+3)/2,2*q) )
      ALLOCATE( forv(l*(l+3)/2,2*q) )
      ALLOCATE( lambda(l,q) )
      ALLOCATE( normal(l,q) )

!
! Reads the tables from the external files

      IF (myrank.eq.0) THEN
         OPEN(1,file=trim(tdir) // '/lambda.in',form='unformatted')
         READ(1) lambda
         CLOSE(1)
         OPEN(1,file=trim(tdir) // '/normal.in',form='unformatted')
         READ(1) normal
         CLOSE(1)
      ENDIF
      CALL MPI_BCAST(lambda,l*q,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(normal,l*q,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      length = LEN(trim(tdir))
      cname(1:length+10) = trim(tdir) // '/cicoef.in'
      CALL MPI_FILE_OPEN(MPI_COMM_WORLD,trim(tdir) // '/cicoef.in', &
          MPI_MODE_RDONLY,MPI_INFO_NULL,fh,ierr)
      CALL MPI_FILE_SET_VIEW(fh,disp,MPI_DOUBLE_COMPLEX,itype,'native', &
          MPI_INFO_NULL,ierr)
      CALL MPI_FILE_READ_ALL(fh,cicoef, &
          4*(l+2)*l**2*q**2*(iend-ista+1)*(qend-qsta+1), &
          MPI_DOUBLE_COMPLEX,MPI_STATUS_IGNORE,ierr)
      CALL MPI_FILE_CLOSE(fh,ierr)
      IF (corio.eq.1) THEN
         cname(1:length+10) = trim(tdir) // '/corixp.in'
         CALL MPI_FILE_OPEN(MPI_COMM_WORLD,trim(tdir) // '/corixp.in', &
             MPI_MODE_RDONLY,MPI_INFO_NULL,fh,ierr)
         CALL MPI_FILE_SET_VIEW(fh,disp,MPI_DOUBLE_COMPLEX,ctype,'native', &
             MPI_INFO_NULL,ierr)
         CALL MPI_FILE_READ_ALL(fh,corixp, &
             2*l*q*(iend-ista+1)*(qend-qsta+1), &
             MPI_DOUBLE_COMPLEX,MPI_STATUS_IGNORE,ierr)
         CALL MPI_FILE_CLOSE(fh,ierr)
         cname(1:length+10) = trim(tdir) // '/corixm.in'
         CALL MPI_FILE_OPEN(MPI_COMM_WORLD,trim(tdir) // '/corixm.in', &
             MPI_MODE_RDONLY,MPI_INFO_NULL,fh,ierr)
         CALL MPI_FILE_SET_VIEW(fh,disp,MPI_DOUBLE_COMPLEX,ctype,'native', &
             MPI_INFO_NULL,ierr)
         CALL MPI_FILE_READ_ALL(fh,corixm, &
             2*l*q*(iend-ista+1)*(qend-qsta+1), &
             MPI_DOUBLE_COMPLEX,MPI_STATUS_IGNORE,ierr)
         CALL MPI_FILE_CLOSE(fh,ierr)
         cname(1:length+10) = trim(tdir) // '/corioz.in'
         CALL MPI_FILE_OPEN(MPI_COMM_WORLD,trim(tdir) // '/corioz.in', &
             MPI_MODE_RDONLY,MPI_INFO_NULL,fh,ierr)
         CALL MPI_FILE_SET_VIEW(fh,disp,MPI_DOUBLE_COMPLEX,ctype,'native', &
             MPI_INFO_NULL,ierr)
         CALL MPI_FILE_READ_ALL(fh,corioz, &
             2*l*q*(iend-ista+1)*(qend-qsta+1), &
             MPI_DOUBLE_COMPLEX,MPI_STATUS_IGNORE,ierr)
         CALL MPI_FILE_CLOSE(fh,ierr)
      ENDIF

!
! Sets the initial conditions. If stat is equal
! to zero a new run is started. In any other 
! case, a previous run is continued.

      timef = 0
      fstep = INT(cort/dt)
      INCLUDE 'initialf.f90'            ! external force

      IF (stat.eq.0) THEN

         ini = 1
         sindex = 0
         bindex = 0
         timeg = gstep
         times = sstep
         timeb = bstep
         INCLUDE 'initialv.f90'         ! initial conditions for v

      ELSE

         ini = INT(stat*bstep)+1
         bindex = stat
         sindex = INT(FLOAT(ini-1)/FLOAT(sstep))
         timeg = 0
         times = 0
         timeb = bstep
         IF (myrank.eq.0) THEN
            CALL genext(bindex,ext)
            OPEN(1,file=trim(odir) // '/xsiv.' // ext // '.out', &
                 form='unformatted')
            READ(1) xsiv
            CLOSE(1)
         ENDIF
         CALL MPI_BCAST(xsiv,(l+3)*l*q,MPI_DOUBLE_COMPLEX,0, &
              MPI_COMM_WORLD,ierr)
         bindex = bindex+1

      ENDIF

!
! Time integration scheme starts here

      IF (bench.eq.1) THEN
         CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
         CALL CPU_Time(cputime1)
      ENDIF

      xsi1 = xsiv
 RK : DO t = ini,step

! Every 'fstep' steps, updates the phase 
! of the external force

      IF ((timef.eq.fstep).and.(rand.eq.1)) THEN
         timef = 0
         IF (myrank.eq.0) phase = 2*pi*p0*randu(seed)
         CALL MPI_BCAST(phase,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
         angle = COS(phase)+IM*SIN(phase)
         DO q1 = qsta,qend
            DO ind1 = ista,iend
               forv(ind1,q1) = forv(ind1,q1)*angle
            END DO
         END DO
      ENDIF
      timef = timef+1

 BM : IF (bench.eq.0) THEN              ! skips output if doing benchmarks

! Every 'gstep' steps, generates external files 
! with global quantities

      IF (timeg.eq.gstep) THEN
         timeg = 0
         IF (myrank.eq.0) THEN
            CALL global((t-1)*dt,lambda,xsiv,0)
            CALL injections((t-1)*dt,lambda,xsiv,normal,forv,wtx,wty,wtz)
            IF (angu.eq.1) THEN
               CALL momentum(lambda,normal,xsiv,0)
            ENDIF
         ENDIF
      ENDIF
      timeg = timeg+1

! Every 'sstep' steps, generates external files 
! with the power spectrum

      IF (times.eq.sstep) THEN
         times = 0
         IF (myrank.eq.0) THEN
            CALL genext(sindex,ext)
            CALL spectruml(lambda,xsiv,ext,0)
            CALL spectrumq(lambda,xsiv,ext,0)
         ENDIF
         sindex = sindex+1
      ENDIF
      times = times+1

! Every 'bstep' steps, stores the results of the integration

      IF (timeb.eq.bstep) THEN
         timeb = 0
         IF (myrank.eq.0) THEN
            CALL genext(bindex,ext)
            OPEN(1,file=trim(odir) // '/xsiv.' // ext // '.out', &
                 form='unformatted')
            WRITE(1) xsiv
            CLOSE(1)
            IF (outs.eq.1) THEN
               CALL xsi2ck(xsiv,lambda,normal,edge,grid,odir,ext,0)
            ENDIF
         ENDIF
         bindex = bindex+1
      ENDIF
      timeb = timeb+1

      ENDIF BM

! Runge-Kutta of order 'ord'

 RD : DO o = ord,1,-1

      ! Update angular velocity if doing precession run
      IF (precs.eq.1) THEN
         INCLUDE 'initialwt.f90'
      ENDIF

      xss1 = xsi1
      DO q1 = qsta,qend
      DO ind1 = ista,iend
         l1 = INT(.5*(SQRT(1+8*FLOAT(ind1))-1))
         IF (q1.le.q) THEN
            lam1 = lambda(l1,q1)
         ELSE
            lam1 = -lambda(l1,2*q-q1+1)
         ENDIF
         m1 = ind1-l1*(l1+1)/2          ! only positive values of m1

         nonlin = 0.d0
 NL :    DO q2 = 1,2*q
         DO q3 = 1,2*q
         DO l2 = 1,l
         DO m2 = -l2,l2
         ind2 = l2*(l2+1)+m2
         IF (m2.ge.0) THEN              ! xsi(-m2) = (-1)^m2*CONJG(xsi(m2))
            xv2 = xss1(l2*(l2+1)/2+m2,q2)
         ELSE
            xv2 = (-1)**ABS(m2)*CONJG(xss1(l2*(l2+1)/2-m2,q2))
         ENDIF
         m3 = m1-m2                     ! orthogonality in m
         DO l3 = 1,l
 OM :    IF (ABS(m3).le.l3) THEN        ! -l<=m<=l
            IF (m3.ge.0) THEN           ! xsi(-m3) = (-1)^m3*CONJG(xsi(m3))
               xv3 = xss1(l3*(l3+1)/2+m3,q3)
            ELSE
               xv3 = (-1)**ABS(m3)*CONJG(xss1(l3*(l3+1)/2-m3,q3))
            ENDIF
            IF (q3.le.q) THEN
               lam3 = lambda(l3,q3)
            ELSE
               lam3 = -lambda(l3,2*q-q3+1)
            ENDIF
            nonlin = nonlin+lam3*cicoef(l3,ind2,q3,q2,ind1,q1)*xv2*xv3
         ENDIF OM
         END DO
         END DO
         END DO
         END DO
         END DO NL

         coriolis = 0.d0
 CO :    IF (corio.eq.1) THEN
         DO q2 = 1,2*q
         DO l2 = 1,l
            m2 = m1                     ! orthogonality in m for Iz
 IZ :       IF (ABS(m2).le.l2) THEN     ! -l<=m<=l
            IF (m2.ge.0) THEN           ! xsi(-m2) = (-1)^m2*CONJG(xsi(m2))
               xv2 = xss1(l2*(l2+1)/2+m2,q2)
            ELSE
               xv2 = (-1)**ABS(m2)*CONJG(xss1(l2*(l2+1)/2-m2,q2))
            ENDIF
            coriolis = coriolis+wz*corioz(l2,q2,ind1,q1)*xv2
            ENDIF IZ
            m2 = m1+1                   ! orthogonality in m for Ixp and Iyp
 IP :       IF (ABS(m2).le.l2) THEN     ! -l<=m<=l
            IF (m2.ge.0) THEN           ! xsi(-m2) = (-1)^m2*CONJG(xsi(m2))
               xv2 = xss1(l2*(l2+1)/2+m2,q2)
            ELSE
               xv2 = (-1)**ABS(m2)*CONJG(xss1(l2*(l2+1)/2-m2,q2))
            ENDIF
            coriolis = coriolis+(wx+IM*wy)*corixp(l2,q2,ind1,q1)*xv2
            ENDIF IP
            m2 = m1-1                   ! orthogonality in m for Ixm and Iym
 LM :       IF (ABS(m2).le.l2) THEN     ! -l<=m<=l
            IF (m2.ge.0) THEN           ! xsi(-m2) = (-1)^m2*CONJG(xsi(m2))
               xv2 = xss1(l2*(l2+1)/2+m2,q2)
            ELSE
               xv2 = (-1)**ABS(m2)*CONJG(xss1(l2*(l2+1)/2-m2,q2))
            ENDIF
            coriolis = coriolis+(wx-IM*wy)*corixm(l2,q2,ind1,q1)*xv2
            ENDIF LM
         ENDDO
         ENDDO
         ENDIF CO

         precession = 0.d0
 PR :    IF (precs.eq.1) THEN
         ! (l=1,m=0) or (l=1,m=1)
            IF ((ind1.eq.1).or.(ind1.eq.2)) THEN
               IF (q1.le.q) THEN
                  norm = normal(ind1,q1)
               ELSE
                  norm = normal(ind1,2*q-q1+1)
               ENDIF
               CALL sphbes(1,ABS(lam1),sj,sjp)
               ampl = 4*SQRT(PI/3)*norm*sjp*ABS(lam1)/lam1
               IF (ind1.eq.1) precession =  ampl*wtz
               IF (ind1.eq.2) precession = -ampl*(wtx-IM*wty)/SQRT(2.0)
            ENDIF
         ENDIF PR

         IF (m1.eq.0) THEN
            nonlin = DBLE(nonlin)
            coriolis = DBLE(coriolis)
            precession = DBLE(precession)
         ENDIF
         xsi1(ind1,q1) = xsiv(ind1,q1)+dt*(nonlin+coriolis+precession- &
                         nu*lam1**2*xss1(ind1,q1)+forv(ind1,q1))/o

      ENDDO
      ENDDO
 SY : DO irank = 0,nprocs-1             ! synchronization
         CALL MPI_BCAST(xsi1,1,xtype(irank),irank,MPI_COMM_WORLD,ierr)
      ENDDO SY
      ENDDO RD
      xsiv = xsi1

      ENDDO RK

!
! Computes the benchmark

      IF (bench.eq.1) THEN
         CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
         CALL CPU_Time(cputime2)
         IF (myrank.eq.0) THEN
            OPEN(1,file='benchmark.txt',position='append')
            WRITE(1,*) nprocs,(cputime2-cputime1)/(step-ini+1)
            CLOSE(1)
         ENDIF
      ENDIF

!
! End of SPHERE_HD

      IF (corio.eq.1) DEALLOCATE ( corixp,corixm,corioz )
      DEALLOCATE ( xsiv,xsi1,xss1 )
      DEALLOCATE ( cicoef )
      DEALLOCATE ( lambda )
      DEALLOCATE ( normal )
      DEALLOCATE ( xtype )
      CALL MPI_FINALIZE(ierr)

      END PROGRAM SPHERE_HD
