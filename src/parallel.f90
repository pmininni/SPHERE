!=================================================================
! PARALLEL subroutines
!
! Subroutines to distribute the coupling coefficients between 
! different processors using MPI. Since the number of elements 
! in the array containing the coupling coefficients scales as 
! 4(L+2)(L+3)L^3*Q^3 (L and Q are respectively the maximum 
! values of l and q) and the number of expanding coefficients 
! for the fields scales as (L+3)L*Q, only the coupling 
! coefficients are distributed.
!
! 2005 Pablo D. Mininni.
!      National Center for Atmospheric Research.
!      e-mail: mininni@ucar.edu 
!=================================================================

!*****************************************************************
      SUBROUTINE distribute(l,q,nprocs,irank,ista,iend,qsta,qend)
!-----------------------------------------------------------------
!
! Distributes the last two indices of the array containing 
! the coupling coefficients in nproc processors. If nprocs 
! is smaller than 2*q, only the last index of the array is 
! distributed. In any other case, the distribution takes 
! place in both directions.
!
! Parameters
!     l     : the maximum value of l
!     q     : the maximum value of q
!     nprocs: the number of processors
!     irank : the rank of the processor
!     ista  : start value for the local l and m coordinates
!     iend  : end value for the local l and m coordinates
!     qsta  : start value for the local q-coordinate
!     qend  : end value for the local q-coordinate
!-----------------------------------------------------------------

      IMPLICIT NONE

      INTEGER, INTENT(IN)   :: l,q
      INTEGER, INTENT(IN)   :: nprocs
      INTEGER, INTENT(IN)   :: irank
      INTEGER, INTENT(OUT)  :: ista,iend
      INTEGER, INTENT(OUT)  :: qsta,qend
      INTEGER               :: lrank,qrank
      INTEGER               :: nl,nq

      IF (nprocs.le.2*q) THEN
         ista = 1
         iend = l*(l+3)/2
         CALL range(1,2*q,nprocs,irank,qsta,qend)
      ELSE
         CALL factorize(nprocs,q,nq,nl)
         lrank = INT(FLOAT(irank)/nq)
         qrank = irank-lrank*nq
         CALL range(1,l*(l+3)/2,nl,lrank,ista,iend)
         CALL range(1,2*q,nq,qrank,qsta,qend)
      ENDIF

      RETURN
      END SUBROUTINE distribute

!*****************************************************************
      SUBROUTINE create_iotype(l,q,ista,iend,qsta,qend,itype)
!-----------------------------------------------------------------
!
! Creates an MPI derived data type for MPI I/O of the 
! distributed coupling coefficients array. 
!
! Parameters
!     l     : the maximum value of l
!     q     : the maximum value of q
!     ista  : start value for the local l and m coordinates
!     iend  : end value for the local l and m coordinates
!     qsta  : start value for the local q-coordinate
!     qend  : end value for the local q-coordinate
!     itype : MPI data type for the block
!-----------------------------------------------------------------


      IMPLICIT NONE
      INCLUDE 'mpif.h'

      INTEGER, INTENT(IN)   :: l,q
      INTEGER, INTENT(IN)   :: ista,iend
      INTEGER, INTENT(IN)   :: qsta,qend
      INTEGER, INTENT(OUT)  :: itype
      INTEGER, DIMENSION(6) :: sizes,subsizes,starts
      INTEGER               :: ierr

      sizes(1) = l
      sizes(2) = l*(l+2)
      sizes(3) = 2*q
      sizes(4) = 2*q
      sizes(5) = l*(l+3)/2
      sizes(6) = 2*q
      subsizes(1) = l
      subsizes(2) = l*(l+2)
      subsizes(3) = 2*q
      subsizes(4) = 2*q
      subsizes(5) = iend-ista+1
      subsizes(6) = qend-qsta+1
      starts(1) = 0
      starts(2) = 0
      starts(3) = 0
      starts(4) = 0
      starts(5) = ista-1
      starts(6) = qsta-1
      CALL MPI_TYPE_CREATE_SUBARRAY(6,sizes,subsizes,starts, &
           MPI_ORDER_FORTRAN,MPI_DOUBLE_COMPLEX,itype,ierr)
      CALL MPI_TYPE_COMMIT(itype,ierr)

      RETURN
      END SUBROUTINE create_iotype

!*****************************************************************
      SUBROUTINE create_iotypc(l,q,ista,iend,qsta,qend,itype)
!-----------------------------------------------------------------
!
! Creates an MPI derived data type for MPI I/O of the 
! distributed Coriolis coefficients array. 
!
! Parameters
!     l     : the maximum value of l
!     q     : the maximum value of q
!     ista  : start value for the local l and m coordinates
!     iend  : end value for the local l and m coordinates
!     qsta  : start value for the local q-coordinate
!     qend  : end value for the local q-coordinate
!     itype : MPI data type for the block
!-----------------------------------------------------------------


      IMPLICIT NONE
      INCLUDE 'mpif.h'

      INTEGER, INTENT(IN)   :: l,q
      INTEGER, INTENT(IN)   :: ista,iend
      INTEGER, INTENT(IN)   :: qsta,qend
      INTEGER, INTENT(OUT)  :: itype
      INTEGER, DIMENSION(4) :: sizes,subsizes,starts
      INTEGER               :: ierr

      sizes(1) = l
      sizes(2) = 2*q
      sizes(3) = l*(l+3)/2
      sizes(4) = 2*q
      subsizes(1) = l
      subsizes(2) = 2*q
      subsizes(3) = iend-ista+1
      subsizes(4) = qend-qsta+1
      starts(1) = 0
      starts(2) = 0
      starts(3) = ista-1
      starts(4) = qsta-1
      CALL MPI_TYPE_CREATE_SUBARRAY(4,sizes,subsizes,starts, &
           MPI_ORDER_FORTRAN,MPI_DOUBLE_COMPLEX,itype,ierr)
      CALL MPI_TYPE_COMMIT(itype,ierr)

      RETURN
      END SUBROUTINE create_iotypc

!*****************************************************************
      SUBROUTINE create_block(l,q,nprocs,itype)
!-----------------------------------------------------------------
!
! Creates an array of MPI derived data types for communication 
! between processors of the expanding coefficients array. 
!
! Parameters
!     nprocs: the number of processors
!     itype : the MPI derived data type
!-----------------------------------------------------------------

      IMPLICIT NONE
      INCLUDE 'mpif.h'

      INTEGER, INTENT(IN)   :: l,q
      INTEGER, INTENT(IN)   :: nprocs
      INTEGER, INTENT(OUT)  :: itype(0:nprocs-1)
      INTEGER, DIMENSION(2) :: sizes,subsizes,starts
      INTEGER               :: ista,iend
      INTEGER               :: qsta,qend
      INTEGER               :: itemp
      INTEGER               :: irank
      INTEGER               :: ierr

      DO irank = 0,nprocs-1
         CALL distribute(l,q,nprocs,irank,ista,iend,qsta,qend)
         sizes(1) = l*(l+3)/2
         sizes(2) = 2*q
         subsizes(1) = iend-ista+1
         subsizes(2) = qend-qsta+1
         starts(1) = ista-1
         starts(2) = qsta-1
         CALL MPI_TYPE_CREATE_SUBARRAY(2,sizes,subsizes,starts, &
              MPI_ORDER_FORTRAN,MPI_DOUBLE_COMPLEX,itemp,ierr)
         CALL MPI_TYPE_COMMIT(itemp,ierr)
         itype(irank) = itemp
      ENDDO

      RETURN
      END SUBROUTINE create_block

!*****************************************************************
      SUBROUTINE range(n1,n2,nprocs,irank,ista,iend)
!-----------------------------------------------------------------
!
! Soubroutine to compute the local coordinate range 
! when splitting the original array into the nodes
!
! Parameters
!     n1    : the minimum value in the splitted dimension
!     n2    : the maximum value in the splitted dimension
!     nprocs: the number of processors
!     irank : the rank of the processor
!     ista  : start value for the local coordinate
!     iend  : end value for the local coordinate
!-----------------------------------------------------------------

      IMPLICIT NONE

      INTEGER, INTENT(OUT) :: ista,iend
      INTEGER, INTENT(IN)  :: n1,n2,irank
      INTEGER, INTENT(IN)  :: nprocs
      INTEGER              :: myrank
      INTEGER              :: iwork1,iwork2

      iwork1 = (n2-n1+1)/nprocs
      iwork2 = MOD(n2-n1+1,nprocs)
      ista = irank*iwork1+n1+MIN(irank,iwork2)
      iend = ista+iwork1-1
      IF (iwork2.gt.irank) iend = iend+1

      RETURN
      END SUBROUTINE range

!*****************************************************************
      SUBROUTINE factorize(nprocs,q,nq,nl)
!-----------------------------------------------------------------
!
! Writes the total number of processors as a product of two 
! integer numbers nq and nl, where nq is smaller or equal to 
! 2*q.
!
! Parameters
!     nprocs: number of processors
!     q     : maximum value of q
!     nq    : number of blocks in the q direction
!     nl    : number of blocks in the l,m direction
!-----------------------------------------------------------------

      IMPLICIT NONE

      REAL                 :: df,rn
      INTEGER, INTENT(OUT) :: nq,nl
      INTEGER, INTENT(IN)  :: nprocs,q
      INTEGER              :: i,in

      rn = FLOAT(nprocs)
      nq = 1
      nl = 1
      df = rn-2*INT(rn/2)
      DO WHILE (df.eq.0)                ! multiple of two
         IF (nq.le.q) THEN
            nq = nq*2
         ELSE
            nl = nl*2
         ENDIF
         rn=rn/2
         df = rn-2*INT(rn/2)
      ENDDO
      df = rn-3*INT(rn/3)
      DO WHILE (df.eq.0)                ! multiple of three
         IF (3*nq.le.2*q) THEN
            nq = nq*3
         ELSE
            nl = nl*3
         ENDIF
         rn = rn/3
         df = rn-3*INT(rn/3)
      END DO
      DO i = 6,INT(SQRT(rn))+1,6        ! 6i-1 and 6i+1 up to sqrt(n)
         df = rn-(i-1)*INT(rn/(i-1))
         DO WHILE (df.eq.0)
            IF ((i-1)*nq.le.2*q) THEN
               nq = nq*(i-1)
            ELSE
               nl = nl*(i-1)
            ENDIF
            rn = rn/(i-1)
            df = rn-(i-1)*INT(rn/(i-1))
         ENDDO
         df = rn-(i+1)*INT(rn/(i+1))
         DO WHILE (df.eq.0)
            IF ((i+1)*nq.le.2*q) THEN
               nq = nq*(i+1)
            ELSE
               nl = nl*(i+1)
            ENDIF
            rn = rn/(i+1)
            df = rn-(i+1)*INT(rn/(i+1))
         ENDDO
      ENDDO
      IF (rn.gt.1) THEN                 ! multiple of itself
         IF (INT(rn)*nq.le.2*q) THEN
            nq = nq*INT(rn)
         ELSE
            nl = nl*INT(rn)
         ENDIF
      ENDIF

      RETURN
      END SUBROUTINE factorize
