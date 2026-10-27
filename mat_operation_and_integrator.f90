MODULE MATRIX_OPERATOR
CONTAINS

FUNCTION CDET(MAT) RESULT(determinant)
IMPLICIT NONE

COMPLEX(KIND=8),INTENT(IN),DIMENSION(:,:)  :: MAT
COMPLEX(KIND=8) :: DETERMINANT

COMPLEX(KIND=8), ALLOCATABLE :: A(:,:)
INTEGER(KIND=4) :: N, LDA, INFO, COUNTER, I
INTEGER(KIND=4), ALLOCATABLE, DIMENSION(:) :: ipivot
N=SIZE(MAT,1)
lda=SIZE(MAT,1)
ALLOCATE(IPIVOT(N))
ALLOCATE(A(N,N))
A=MAT
IPIVOT=0
CALL zgetrf(n,n,a,lda,ipivot,info)
IF (info /= 0) THEN
  STOP 'NUMERICAL ERROR IN ZRETRF'
END IF
counter=0
determinant=DCMPLX(1.0D0,0.0D0)
DO i=1,n
  IF (ipivot(i).NE.i) counter=counter+1
  determinant=determinant*a(i,i)
ENDDO

IF(MOD(counter,2)==1) determinant=-1.0D0*determinant
DEALLOCATE(IPIVOT,A)
END FUNCTION CDET

FUNCTION cinv(A) RESULT(INVERSE_A)
implicit none

integer(KIND=4) :: N_dim
complex(KIND=8), intent(in),DIMENSION(:,:)  :: A
complex(KIND=8) :: inverse_A(SIZE(A,1),SIZE(A,1))
logical :: success

integer(kind=4) :: dim4
integer(kind=4) :: info
integer(kind=4), allocatable :: IPIV(:)
complex(KIND=8), allocatable :: copy_of_A(:,:), work(:)
N_DIM=SIZE(A,1)
dim4 = N_DIM
allocate(copy_of_A(N_DIM,N_DIM))
allocate(work(N_DIM))
allocate(ipiv(N_DIM))

copy_of_A=A

Call ZGETRF( dim4, dim4, copy_of_A, dim4, IPIV, INFO )

call ZGETRI( dim4, copy_of_A, dim4, IPIV, WORK, dim4, INFO )

inverse_A=copy_of_A
if (info ==0 ) then; success=.true. ; endif
if (info /=0 ) then; success=.false. ; endif
IF (.NOT. SUCCESS) THEN
  WRITE(*,*) "CALLING ZGETRI WRONG IN FUNCTION CINV."

END IF
deallocate(IPIV, COPY_OF_A, WORK)
end FUNCTION CINV

complex(KIND=8) function ctra(A) result(TraceA)
implicit none
complex(KIND=8), dimension(:,:), intent(in) :: A
integer :: i,n
n=size(A,1)
TraceA=0
do i=1,n
  TraceA=TraceA+A(i,i)
end do
end function ctra

FUNCTION CLOGDET(MAT) RESULT(LOGDET)
IMPLICIT NONE

COMPLEX(KIND=8),INTENT(IN),DIMENSION(:,:)  :: MAT
COMPLEX(KIND=8) :: LOGDET

COMPLEX(KIND=8), ALLOCATABLE :: A(:,:)
INTEGER(KIND=4) :: N, LDA, INFO, COUNTER, I
INTEGER(KIND=4), ALLOCATABLE, DIMENSION(:) :: ipivot
N=SIZE(MAT,1)
lda=SIZE(MAT,1)
ALLOCATE(IPIVOT(N))
ALLOCATE(A(N,N))
A=MAT
IPIVOT=0
CALL zgetrf(n,n,a,lda,ipivot,info)
IF (info /= 0) THEN
  STOP 'NUMERICAL ERROR IN ZRETRF'
END IF
counter=0
LOGDET=0
DO i=1,n
  IF (ipivot(i).NE.i) counter=counter+1
  LOGDET = LOGDET + LOG(A(I,I))
ENDDO

IF(MOD(counter,2)==1) LOGDET = LOGDET + DCMPLX(0.D0,3.1415926535897932384)
LOGDET = DCMPLX(REAL(LOGDET),MOD(AIMAG(LOGDET),2*3.1415926535897932384))
DEALLOCATE(IPIVOT,A)
END FUNCTION CLOGDET
END MODULE MATRIX_OPERATOR
