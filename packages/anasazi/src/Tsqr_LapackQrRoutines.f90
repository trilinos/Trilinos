module LapackQrRoutines
  implicit none

  interface 
     subroutine DLARFP( n, alpha, x, incx, tau )
       integer, intent(in)     :: n, incx
       real (8), intent(inout) :: alpha
       real (8), intent(out)   :: tau
       real (8), intent(inout) :: x(*)
     end subroutine DLARFP

     subroutine SLARFP( n, alpha, x, incx, tau )
       integer, intent(in)     :: n, incx
       real (4), intent(inout) :: alpha
       real (4), intent(out)   :: tau
       real (4), intent(inout) :: x(*)
     end subroutine SLARFP

     subroutine DLARFG( n, alpha, x, incx, tau )
       integer, intent(in)     :: n, incx
       real (8), intent(inout) :: alpha
       real (8), intent(out)   :: tau
       real (8), intent(inout) :: x(*)
     end subroutine DLARFG

     subroutine SLARFG( n, alpha, x, incx, tau )
       integer, intent(in)     :: n, incx
       real (4), intent(inout) :: alpha
       real (4), intent(out)   :: tau
       real (4), intent(inout) :: x(*)
     end subroutine SLARFG
     
     subroutine DGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
       real (8), intent(in)    :: ALPHA,BETA
       integer, intent(in)     :: INCX,INCY,LDA,M,N
       character, intent(in)   :: TRANS
       real (8), intent(in)    :: A(LDA,*), X(*)
       real (8), intent(inout) :: Y(*)
     end subroutine DGEMV

     subroutine SGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
       real (4), intent(in)    :: ALPHA,BETA
       integer, intent(in)     :: INCX,INCY,LDA,M,N
       character, intent(in)   :: TRANS
       real (4), intent(in)    :: A(LDA,*), X(*)
       real (4), intent(inout) :: Y(*)
     end subroutine SGEMV
     
     subroutine DGER(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)
       real (8), intent(in)    :: ALPHA
       integer, intent(in)     :: INCX,INCY,LDA,M,N
       real (8), intent(inout) :: A(LDA,*)
       real (8), intent(in)    :: X(*),Y(*)
     end subroutine DGER

     subroutine SGER(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)
       real (4), intent(in)    :: ALPHA
       integer, intent(in)     :: INCX,INCY,LDA,M,N
       real (4), intent(inout) :: A(LDA,*)
       real (4), intent(in)    :: X(*),Y(*)
     end subroutine SGER

     subroutine DGEQR2( M, N, A, LDA, TAU, WORK, INFO )
       integer, intent(in)     :: M, N, LDA
       real (8), intent(inout) :: A(LDA,*)
       real (8), intent(out)   :: TAU(*), WORK(*)
       integer, intent(out)    :: INFO
     end subroutine DGEQR2

     subroutine SGEQR2( M, N, A, LDA, TAU, WORK, INFO )
       integer, intent(in)     :: M, N, LDA
       real (4), intent(inout) :: A(LDA,*)
       real (4), intent(out)   :: TAU(*), WORK(*)
       integer, intent(out)    :: INFO
     end subroutine SGEQR2

     subroutine DGEQRF( M, N, A, LDA, TAU, WORK, LWORK, INFO )
       integer, intent(in)     :: M, N, LDA, LWORK
       real (8), intent(inout) :: A(LDA,*)
       real (8), intent(out)   :: TAU(*), WORK(*)
       integer, intent(out)    :: INFO
     end subroutine DGEQRF

     subroutine SGEQRF( M, N, A, LDA, TAU, WORK, LWORK, INFO )
       integer, intent(in)     :: M, N, LDA, LWORK
       real (4), intent(inout) :: A(LDA,*)
       real (4), intent(out)   :: TAU(*), WORK(*)
       integer, intent(out)    :: INFO
     end subroutine SGEQRF

     subroutine DORMQR( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK, LWORK, INFO )
       character, intent(in) :: SIDE, TRANS
       integer, intent(in)     :: M, N, K, LDA, LDC, LWORK
       real (8), intent(in)    :: A(LDA,*), TAU(*)
       real (8), intent(inout) :: C(LDC,*)
       real (8), intent(out)   :: WORK(*)
       integer, intent(out)    :: INFO
     end subroutine DORMQR

     subroutine SORMQR( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK, LWORK, INFO )
       character, intent(in) :: SIDE, TRANS
       integer, intent(in)     :: M, N, K, LDA, LDC, LWORK
       real (4), intent(in)    :: A(LDA,*), TAU(*)
       real (4), intent(inout) :: C(LDC,*)
       real (4), intent(out)   :: WORK(*)
       integer, intent(out)    :: INFO
     end subroutine SORMQR

     subroutine DORM2R( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK, INFO )
       character, intent(in) :: SIDE, TRANS
       integer, intent(in)     :: M, N, K, LDA, LDC
       real (8), intent(in)    :: A(LDA,*), TAU(*)
       real (8), intent(inout) :: C(LDC,*)
       real (8), intent(out)   :: WORK(*)
       integer, intent(out)    :: INFO
     end subroutine DORM2R

     subroutine SORM2R( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK, INFO )
       character, intent(in) :: SIDE, TRANS
       integer, intent(in)     :: M, N, K, LDA, LDC
       real (4), intent(in)    :: A(LDA,*), TAU(*)
       real (4), intent(inout) :: C(LDC,*)
       real (4), intent(out)   :: WORK(*)
       integer, intent(out)    :: INFO
     end subroutine SORM2R

     subroutine DORGQR( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
       integer, intent(in)     :: M, N, K, LDA, LWORK
       real (8), intent(inout) :: A(LDA,*)
       real (8), intent(in)    :: TAU(*)
       real (8), intent(out)   :: WORK(*)
       integer, intent(out)    :: INFO
     end subroutine DORGQR

     subroutine SORGQR( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
       integer, intent(in)     :: M, N, K, LDA, LWORK
       real (4), intent(inout) :: A(LDA,*)
       real (4), intent(in)    :: TAU(*)
       real (4), intent(out)   :: WORK(*)
       integer, intent(out)    :: INFO
     end subroutine SORGQR

     subroutine DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
       real (8), intent(in)    :: ALPHA,BETA
       integer, intent(in)     :: K,LDA,LDB,LDC,M,N
       character, intent(in)   :: TRANSA,TRANSB
       real (8), intent(in)    :: A(LDA,*), B(LDB,*)
       real (8), intent(inout) :: C(LDC,*)
     end subroutine DGEMM

     subroutine SGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
       real (4), intent(in)    :: ALPHA,BETA
       integer, intent(in)     :: K,LDA,LDB,LDC,M,N
       character, intent(in)   :: TRANSA,TRANSB
       real (4), intent(in)    :: A(LDA,*), B(LDB,*)
       real (4), intent(inout) :: C(LDC,*)
     end subroutine SGEMM

     subroutine DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )
       character, intent(in)   :: JOBZ, UPLO
       integer, intent(in)     :: LDA, LWORK, N
       real (8), intent(inout) :: A(LDA,*)
       real (8), intent(out)   :: W(*), WORK(*)
       integer, intent(out)    :: INFO
     end subroutine DSYEV

     subroutine SSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )
       character, intent(in)   :: JOBZ, UPLO
       integer, intent(in)     :: LDA, LWORK, N
       real (4), intent(inout) :: A(LDA,*)
       real (4), intent(out)   :: W(*), WORK(*)
       integer, intent(out)    :: INFO
     end subroutine SSYEV

     real(8) function DNRM2( N, X, INCX )
       integer, intent(in)  :: N, INCX
       real (8), intent(in) :: X(*)
     end function DNRM2

     real(4) function SNRM2( N, X, INCX )
       integer, intent(in)  :: N, INCX
       real (4), intent(in) :: X(*)
     end function SNRM2

     real(8) function DDOT( N, DX, INCX, DY, INCY )
       integer, intent(in)  :: N, INCX, INCY
       real (8), intent(in) :: DX(*), DY(*)
     end function DDOT

     real(4) function SDOT( N, DX, INCX, DY, INCY )
       integer, intent(in)  :: N, INCX, INCY
       real (4), intent(in) :: DX(*), DY(*)
     end function SDOT

     integer function IDAMAX( N, X, INCX )
       integer, intent(in)  :: N, INCX
       real (8), intent(in) :: X(*)
     end function IDAMAX

     integer function ISAMAX( N, X, INCX )
       integer, intent(in)  :: N, INCX
       real (4), intent(in) :: X(*)
     end function ISAMAX

     subroutine DGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, INFO )
       character, intent(in) :: JOBU, JOBVT
       integer, intent(in)   :: LDA, LDU, LDVT, LWORK, M, N
       integer, intent(out)  :: INFO
       real (8), intent(inout) :: A(LDA,*)
       real (8), intent(out) :: S(*), U(LDU,*), VT(LDVT,*), WORK(*)
     end subroutine DGESVD

     subroutine SGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, INFO )
       character, intent(in) :: JOBU, JOBVT
       integer, intent(in)   :: LDA, LDU, LDVT, LWORK, M, N
       integer, intent(out)  :: INFO
       real (4), intent(inout) :: A(LDA,*)
       real (4), intent(out) :: S(*), U(LDU,*), VT(LDVT,*), WORK(*)
     end subroutine SGESVD
  end interface

end module LapackQrRoutines
