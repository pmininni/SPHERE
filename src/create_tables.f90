!=================================================================
      PROGRAM CREATE_TABLES
!=================================================================
! CREATE_TABLES code
!
! Generates binary files with tables of the roots of the 
! spherical bessel functions, normalization coefficients 
! for the spherical Chandrasekhar-Kendall functions, and 
! coupling coefficients for the Navier-Stokes and MHD 
! equations projected into this base. The output is used 
! by the codes SPHERE_HD and SPHERE_MHD.
!
! INPUT : create.inp
! OUTPUT: lambda.in [binary file]
!         normal.in [binary file]
!         cicoef.in [MPI native file]
!
! 2005 Pablo D. Mininni.
!      National Center for Atmospheric Research.
!      e-mail: mininni@ucar.edu 
!=================================================================

      USE mpivars
      USE functions
      USE resolution
      IMPLICIT NONE

!
! Table matrixes

      DOUBLE COMPLEX, ALLOCATABLE, DIMENSION(:,:,:,:,:,:) :: cicoef
      DOUBLE COMPLEX, ALLOCATABLE, DIMENSION(:,:,:,:)     :: corixp
      DOUBLE COMPLEX, ALLOCATABLE, DIMENSION(:,:,:,:)     :: corixm
      DOUBLE COMPLEX, ALLOCATABLE, DIMENSION(:,:,:,:)     :: corioz
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:)       :: lambda
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:)       :: normal

      INTEGER            :: corio,bound
      INTEGER            :: length
      CHARACTER(len=110) :: cname
      CHARACTER(len=100) :: odir

!
! Initializes the MPI library

      CALL MPI_INIT(ierr)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)

!
! Reads from the external file 'create.inp'
! information for the present computation
!     bound: if 1 computes coefficients for a potential field
!     corio: if 1 computes coefficients for the Corilis term
!     odir : directory for binary output

      IF (myrank.eq.0) THEN
         OPEN(1,file='create.inp',status='unknown')
         READ(1,*) bound
         READ(1,*) corio
         READ(1,'(a100)') odir
         CLOSE(1)
      ENDIF
      CALL MPI_BCAST(bound,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(corio,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(odir,100,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)

!
! Allocates memory

      CALL distribute(l,q,nprocs,myrank,ista,iend,qsta,qend)
      CALL create_iotype(l,q,ista,iend,qsta,qend,itype)
      ALLOCATE( cicoef(l,l*(l+2),2*q,2*q,ista:iend,qsta:qend) )
      IF (corio.eq.1) THEN
         CALL create_iotypc(l,q,ista,iend,qsta,qend,ctype)
         ALLOCATE( corixp(l,2*q,ista:iend,qsta:qend) )
         ALLOCATE( corixm(l,2*q,ista:iend,qsta:qend) )
         ALLOCATE( corioz(l,2*q,ista:iend,qsta:qend) )
      ENDIF
      ALLOCATE( lambda(l,q) )
      ALLOCATE( normal(l,q) )

!
! Creates the tables

      CALL lambda_table(lambda)
      CALL normal_table(lambda,normal)
      CALL cicoef_table(lambda,normal,cicoef)
      IF (corio.eq.1) CALL coriol_table(lambda,normal,corixp,corixm,corioz)

!
! Exports the results

      IF (myrank.eq.0) THEN
         OPEN(1,file=trim(odir) // '/lambda.in',form='unformatted')
         WRITE(1) lambda
         CLOSE(1)
         OPEN(1,file=trim(odir) // '/normal.in',form='unformatted')
         WRITE(1) normal
         CLOSE(1)
      ENDIF
      length = LEN(trim(odir))
      cname(1:length+10) = trim(odir) // '/cicoef.in'
      CALL MPI_FILE_OPEN(MPI_COMM_WORLD,trim(odir)//'/cicoef.in', &
          MPI_MODE_CREATE+MPI_MODE_WRONLY,MPI_INFO_NULL,fh,ierr)
      CALL MPI_FILE_SET_VIEW(fh,disp,MPI_DOUBLE_COMPLEX,itype,'native', &
          MPI_INFO_NULL,ierr)
      CALL MPI_FILE_WRITE_ALL(fh,cicoef, &
          4*(l+2)*l**2*q**2*(iend-ista+1)*(qend-qsta+1), &
          MPI_DOUBLE_COMPLEX,MPI_STATUS_IGNORE,ierr)
      CALL MPI_FILE_CLOSE(fh,ierr)
      IF (corio.eq.1) THEN
         cname(1:length+10) = trim(odir) // '/corixp.in'
         CALL MPI_FILE_OPEN(MPI_COMM_WORLD,trim(odir)//'/corixp.in', &
             MPI_MODE_CREATE+MPI_MODE_WRONLY,MPI_INFO_NULL,fh,ierr)
         CALL MPI_FILE_SET_VIEW(fh,disp,MPI_DOUBLE_COMPLEX,ctype,'native', &
             MPI_INFO_NULL,ierr)
         CALL MPI_FILE_WRITE_ALL(fh,corixp, &
             2*l*q*(iend-ista+1)*(qend-qsta+1), &
             MPI_DOUBLE_COMPLEX,MPI_STATUS_IGNORE,ierr)
         CALL MPI_FILE_CLOSE(fh,ierr)
         cname(1:length+10) = trim(odir) // '/corixm.in'
         CALL MPI_FILE_OPEN(MPI_COMM_WORLD,trim(odir)//'/corixm.in', &
             MPI_MODE_CREATE+MPI_MODE_WRONLY,MPI_INFO_NULL,fh,ierr)
         CALL MPI_FILE_SET_VIEW(fh,disp,MPI_DOUBLE_COMPLEX,ctype,'native', &
             MPI_INFO_NULL,ierr)
         CALL MPI_FILE_WRITE_ALL(fh,corixm, &
             2*l*q*(iend-ista+1)*(qend-qsta+1), &
             MPI_DOUBLE_COMPLEX,MPI_STATUS_IGNORE,ierr)
         CALL MPI_FILE_CLOSE(fh,ierr)
         cname(1:length+10) = trim(odir) // '/corioz.in'
         CALL MPI_FILE_OPEN(MPI_COMM_WORLD,trim(odir)//'/corioz.in', &
             MPI_MODE_CREATE+MPI_MODE_WRONLY,MPI_INFO_NULL,fh,ierr)
         CALL MPI_FILE_SET_VIEW(fh,disp,MPI_DOUBLE_COMPLEX,ctype,'native', &
             MPI_INFO_NULL,ierr)
         CALL MPI_FILE_WRITE_ALL(fh,corioz, &
             2*l*q*(iend-ista+1)*(qend-qsta+1), &
             MPI_DOUBLE_COMPLEX,MPI_STATUS_IGNORE,ierr)
         CALL MPI_FILE_CLOSE(fh,ierr)
      ENDIF

!
! End of CREATE_TABLES

      CALL MPI_FINALIZE(ierr)
      IF (corio.eq.1) DEALLOCATE ( corixp,corixm,corioz )
      DEALLOCATE ( cicoef )
      DEALLOCATE ( lambda )
      DEALLOCATE ( normal )

      END PROGRAM CREATE_TABLES
