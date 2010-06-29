module TsqrCombine
  use LapackQrRoutines
  implicit none

contains

  ! Apply the Q factor stored in [R; A] to [C_top; C_bot].  The C
  ! blocks are allowed, but not required, to have different leading
  ! dimensions (ldc_top resp. ldc_bottom).  R is upper triangular, so
  ! we do not need it; the Householder reflectors representing the Q
  ! factor are stored compactly in A (specifically, in all of A, not
  ! just the lower triangle).
  !
  ! In the "sequential under parallel" version of TSQR, this function
  ! belongs to the sequential part (i.e., operating on cache blocks on
  ! a single processor).
  !
  ! This routine used to be called "tsqr_apply_inner2" (the original
  ! tsqr_apply_inner assumed that C_top and C_bot had the same leading
  ! dimension, but was otherwise the same).  After that, it was called
  ! "tsqr_apply_inner".  Now we have different function names, based
  ! on input scalar datatypes: s_apply_inner, d_apply_inner, etc.  We
  ! also removed the ISO_C_BINDING dependency, since that doesn't work
  ! with older compilers.
  !
  ! trans [in]     'N' means apply Q, anything else (such as 'T') 
  !                means apply Q^T
  ! m [in]         number of rows of A
  ! ncols_C [in]   number of columns of [C_top; C_bot]
  ! ncols_Q [in]   number of columns of [R; A]
  ! A [in]         m by ncols_Q matrix, in which the Householder 
  !                reflectors representing the Q factor are stored
  ! lda [in]       leading dimension of A
  ! tau [in]       array of length ncols_Q, storing the scaling factors 
  !                for the Householder reflectors representing Q 
  ! C_top [inout]  ncols_Q by ncols_C matrix
  ! ldc_top [in]   leading dimension of C_top
  ! C_bot [inout]  m by ncols_C matrix
  ! ldc_bot [in]   leading dimension of C_bot
  ! work [out]     workspace array of length ncols_C
  ! 
  subroutine d_apply_inner( trans, m, ncols_C, ncols_Q, A, lda, tau, &
       C_top, ldc_top, C_bot, ldc_bot, work )
    implicit none
    
    character, intent(in)        :: trans
    integer, intent(in), value   :: m, ncols_Q, ncols_C, lda, ldc_top, ldc_bot
    real(8), intent(in)          :: A(lda,ncols_Q)
    real(8), intent(inout)       :: C_top(ldc_top,ncols_C), C_bot(ldc_bot,ncols_C)
    real(8), intent(in)          :: tau(ncols_Q)
    real(8), intent(out), target :: work(ncols_C)

    real (8), pointer :: y(:)
    integer           :: j, k

    y => work(1:ncols_C)
    y = 0
    if (trans == 'N' .or. trans == 'n') then
       do k = ncols_Q, 1, -1
          ! The m bottom elements of the length m+1 "sparse" Householder
          ! reflector are stored in A(1:m,k).  The topmost element is
          ! implicit and is 1.0.  If we want to exploit its sparsity, we
          ! can't use DLARF to apply it to [C_top; C_bot]; we have to
          ! roll our own routine, which you see here.
          
          ! $y^T := A(1:m,k)^T C_bot(1:m, 1:ncols_C)$, so $y := C_bot(1:m, 1:ncols_C)^T A(1:m, k)$.
          !call DGEMV( 'T', m, ncols_C, 1.0d0, C_bot(1:m, 1:ncols_C), ldc, A(1:m, k), 1, 0.0d0, y(1:ncols_C), 1 )
          call DGEMV( 'T', m, ncols_C, 1.0d0, C_bot, ldc_bot, A(1, k), 1, 0.0d0, y, 1 )
          
          ! $y^T := y^T + C_top(k, 1:ncols_C)$
          do j = 1, ncols_C
             y(j) = y(j) + C_top(k, j)
          end do
          
          ! Update C_top(k, 1:ncols_C)
          do j = 1, ncols_C
             C_top(k, j) = C_top(k, j) - tau(k) * y(j)
          end do
          
          ! Update C_bot(1:m, 1:ncols_C)
          !call DGER( m, ncols_C, -tau(k), A(1:m,k), 1, y(1:ncols_C), 1, C_bot(1:m, 1:ncols_C), ldc )
          call DGER( m, ncols_C, -tau(k), A(1,k), 1, y, 1, C_bot, ldc_bot )
       end do
    else
       do k = 1, ncols_Q
          ! The m bottom elements of the length m+1 "sparse" Householder
          ! reflector are stored in A(1:m,k).  The topmost element is
          ! implicit and is 1.0.  If we want to exploit its sparsity, we
          ! can't use DLARF to apply it to [C_top; C_bot]; we have to
          ! roll our own routine, which you see here.
          
          ! $y^T := A(1:m,k)^T C_bot(1:m, 1:ncols_C)$, so $y := C_bot(1:m, 1:ncols_C)^T A(1:m, k)$.
          !call DGEMV( 'T', m, ncols_C, 1.0d0, C_bot(1:m, 1:ncols_C), ldc_bot, A(1:m, k), 1, 0.0d0, y(1:ncols_C), 1 )
          call DGEMV( 'T', m, ncols_C, 1.0d0, C_bot, ldc_bot, A(1, k), 1, 0.0d0, y, 1 )
          
          ! $y^T := y^T + C_top(k, 1:ncols_C)$
          do j = 1, ncols_C
             y(j) = y(j) + C_top(k, j)
          end do
          
          ! Update C_top(k, 1:ncols_C)
          do j = 1, ncols_C
             C_top(k, j) = C_top(k, j) - tau(k) * y(j)
          end do
          
          ! Update C_bot(1:m, 1:ncols_C)
          !call DGER( m, ncols_C, -tau(k), A(1:m,k), 1, y(1:ncols_C), 1, C_bot(1:m, 1:ncols_C), ldc_bot )
          call DGER( m, ncols_C, -tau(k), A(1,k), 1, y, 1, C_bot(1, 1), ldc_bot )
       end do
    end if
  end subroutine d_apply_inner


  subroutine s_apply_inner( trans, m, ncols_C, ncols_Q, A, lda, tau, &
       C_top, ldc_top, C_bot, ldc_bot, work )
    implicit none
    
    character, intent(in)        :: trans
    integer, intent(in), value   :: m, ncols_Q, ncols_C, lda, ldc_top, ldc_bot
    real(4), intent(in)          :: A(lda,ncols_Q)
    real(4), intent(inout)       :: C_top(ldc_top,ncols_C), C_bot(ldc_bot,ncols_C)
    real(4), intent(in)          :: tau(ncols_Q)
    real(4), intent(out), target :: work(ncols_C)

    real(4), pointer :: y(:)
    integer          :: j, k

    y => work(1:ncols_C)
    y = 0
    if (trans == 'N' .or. trans == 'n') then
       do k = ncols_Q, 1, -1
          ! The m bottom elements of the length m+1 "sparse"
          ! Householder reflector are stored in A(1:m,k).  The topmost
          ! element is implicit and is 1.0.  If we want to exploit its
          ! sparsity, we can't use SLARF to apply it to [C_top;
          ! C_bot]; we have to roll our own routine, which you see
          ! here.
          !
          ! $y^T := A(1:m,k)^T C_bot(1:m, 1:ncols_C)$, so $y := C_bot(1:m, 1:ncols_C)^T A(1:m, k)$.
          call SGEMV( 'T', m, ncols_C, 1.0, C_bot, ldc_bot, A(1, k), 1, 0.0, y, 1 )
          
          ! $y^T := y^T + C_top(k, 1:ncols_C)$
          do j = 1, ncols_C
             y(j) = y(j) + C_top(k, j)
          end do
          
          ! Update C_top(k, 1:ncols_C)
          do j = 1, ncols_C
             C_top(k, j) = C_top(k, j) - tau(k) * y(j)
          end do
          
          ! Update C_bot(1:m, 1:ncols_C)
          call SGER( m, ncols_C, -tau(k), A(1,k), 1, y, 1, C_bot, ldc_bot )
       end do
    else
       do k = 1, ncols_Q
          ! The m bottom elements of the length m+1 "sparse"
          ! Householder reflector are stored in A(1:m,k).  The topmost
          ! element is implicit and is 1.0.  If we want to exploit its
          ! sparsity, we can't use SLARF to apply it to [C_top;
          ! C_bot]; we have to roll our own routine, which you see
          ! here.
          !
          ! $y^T := A(1:m,k)^T C_bot(1:m, 1:ncols_C)$, so $y := C_bot(1:m, 1:ncols_C)^T A(1:m, k)$.
          call SGEMV( 'T', m, ncols_C, 1.0, C_bot, ldc_bot, A(1, k), 1, 0.0, y, 1 )
          
          ! $y^T := y^T + C_top(k, 1:ncols_C)$
          do j = 1, ncols_C
             y(j) = y(j) + C_top(k, j)
          end do
          
          ! Update C_top(k, 1:ncols_C)
          do j = 1, ncols_C
             C_top(k, j) = C_top(k, j) - tau(k) * y(j)
          end do
          
          ! Update C_bot(1:m, 1:ncols_C)
          call SGER( m, ncols_C, -tau(k), A(1,k), 1, y, 1, C_bot(1, 1), ldc_bot )
       end do
    end if
  end subroutine s_apply_inner

  ! Perform one "inner" QR factorization step of sequential / parallel
  ! TSQR.  (In either case, only one processor calls this routine.)
  !
  ! In the "sequential under parallel" version of TSQR, this function
  ! belongs to the sequential part (i.e., operating on cache blocks on
  ! a single processor).  Only the first cache block $A_0$ is factored
  ! as $Q_0 R_0 = A_0$ (see tsqr_factor_first); subsequent cache blocks
  ! $A_k$ are factored using this routine, which combines them with 
  ! $R_{k-1}$.
  !
  ! Here is the matrix to factor:
  ! \[
  ! \begin{pmatrix}
  ! R_{k-1} \\      % $A_{k-1}$ is $m_{k-1} \times n$ with $m_{k-1} \geq n$
  ! A_k     \\      % $m_k \times n$ with $m_k \geq n$
  ! \end{pmatrix}
  ! \]
  !
  ! Since $R_{k-1}$ is n by n upper triangular, we can store the
  ! Householder reflectors (representing the Q factor of $[R_{k-1};
  ! A_k]$) entirely in $A_k$ (specifically, in all of $A_k$, not just
  ! below the diagonal).
  !
  ! m [in]     Number of rows in the "bottom" block to factor.  The number of
  !            rows in the top block doesn't matter, given the assumptions 
  !            above, as long as $m_{k-1} \geq n$.
  ! n [in]     Number of columns (same in both blocks)
  ! R [inout]  "Top" upper triangular n by n block $R_{k-1}$.
  !            Overwritten with the new R factor $R_k$ of $[R_{k-1}; A_k]$.
  ! ldr [in]   Leading dimension of R
  ! A [inout]  "Bottom" dense m by n block $A_k$.  Overwritten with the 
  !            Householder reflectors representing the Q factor of 
  !            $[R_{k-1}; A_k]$.
  ! tau [out]  Scaling factors of the Householder reflectors.  Corresponds 
  !            to the output of DGEQRF of the same name.
  ! work [out] Workspace (length >= n; don't need lwork or workspace query)
  !
  subroutine d_factor_inner( m, n, R, ldr, A, lda, tau, work )
    implicit none

    integer, value, intent(in)    :: m, n, ldr, lda
    real(8), intent(inout)        :: R(ldr,n)
    real(8), intent(inout)        :: A(lda,n)
    real(8), intent(out)          :: tau(n)
    ! "target" means it's legitimate to have pointers alias this
    ! array, internal to this function.
    real(8), intent(out), target  :: work(n)

    real(8), pointer  :: y(:)
    integer           :: j, k

    y => work(1:n)
    y = 0
    do k = 1, n-1
       ! Form the "sparse" Householder reflector, using DLARFP instead
       ! of DLARFG (need LAPACK 3.2, thanks Jason Riedy!) so that the
       ! diagonal elements of the R factor are nonnegative.
       !call DLARFP( m + 1, R(k,k), A(1:m,k), 1, tau(k) )
       call DLARFP( m + 1, R(k,k), A(1,k), 1, tau(k) )

       ! $y^T := A(1:m,k)^T A(1:m, k+1:n)$, so $y := A(1:m, k+1:n)^T A(1:m, k)$.
       ! BEGIN mfh 07 June 2009
       ! Another segfault here, which could mean that the use of subarrays 
       ! must somehow involve allocating and copying data -- not what we want.
       !call DGEMV( 'T', m, n-k, 1.0d0, A(1:m, k+1:n), lda, A(1:m, k), 1, 0.0d0, y(1:n-k), 1 )
       call DGEMV( 'T', m, n-k, 1.0d0, A(1, k+1), lda, A(1, k), 1, 0.0d0, y, 1 )

       ! $y^T := y^T + R(k, k+1:n)$
       do j = k+1, n
          y(j-k) = y(j-k) + R(k, j)
       end do

       ! Update R(k, k+1:n)
       do j = k+1, n
          R(k, j) = R(k, j) - tau(k) * y(j-k)
       end do

       ! Update A(1:m, k+1:n)
       ! BEGIN mfh 07 June 2009
       ! free() was triggering a segfault here, which means the use of subarrays 
       ! must somehow involve allocating and copying data -- not what we want.
       !call DGER( m, n-k, -tau(k), A(1:m,k), 1, y(1:n-k), 1, A(1:m, k+1:n), lda )
       call DGER( m, n-k, -tau(k), A(1,k), 1, y, 1, A(1, k+1), lda )
       ! END mfh 07 June 2009
    end do

    ! Compute the Householder reflector for the last column.  This
    ! last iteration doesn't require an update of the trailing matrix,
    ! because there is no trailing matrix left!
    !call DLARFP( m + 1, R(n,n), A(1:m,n), 1, tau(n) )
    call DLARFP( m + 1, R(n,n), A(1,n), 1, tau(n) )
  end subroutine d_factor_inner


  subroutine s_factor_inner( m, n, R, ldr, A, lda, tau, work )
    implicit none

    integer, value, intent(in)    :: m, n, ldr, lda
    real(4), intent(inout)        :: R(ldr,n)
    real(4), intent(inout)        :: A(lda,n)
    real(4), intent(out)          :: tau(n)
    ! "target" means it's legitimate to have pointers alias this
    ! array, internal to this function.
    real(4), intent(out), target  :: work(n)

    real(4), pointer  :: y(:)
    integer           :: j, k

    y => work(1:n)
    y = 0
    do k = 1, n-1
       ! Form the "sparse" Householder reflector, using SLARFP instead
       ! of SLARFG (need LAPACK >= 3.2, thanks Jason Riedy!) so that
       ! the diagonal elements of the R factor are nonnegative.
       call SLARFP( m + 1, R(k,k), A(1,k), 1, tau(k) )

       ! $y^T := A(1:m,k)^T A(1:m, k+1:n)$, so $y := A(1:m, k+1:n)^T A(1:m, k)$.
       call SGEMV( 'T', m, n-k, 1.0, A(1, k+1), lda, A(1, k), 1, 0.0, y, 1 )

       ! $y^T := y^T + R(k, k+1:n)$
       do j = k+1, n
          y(j-k) = y(j-k) + R(k, j)
       end do

       ! Update R(k, k+1:n)
       do j = k+1, n
          R(k, j) = R(k, j) - tau(k) * y(j-k)
       end do

       ! Update A(1:m, k+1:n)
       call SGER( m, n-k, -tau(k), A(1,k), 1, y, 1, A(1, k+1), lda )
    end do

    ! Compute the Householder reflector for the last column.  This
    ! last iteration doesn't require an update of the trailing matrix,
    ! because there is no trailing matrix left!
    call SLARFP( m + 1, R(n,n), A(1,n), 1, tau(n) )
  end subroutine s_factor_inner


  ! Compute QR factorization of [R_top; R_bot].  Store resulting R
  ! factor in R_top and Householder reflectors in R_bot.
  !
  ! n [in]         Number of rows and columns of each of R_top and R_bot
  ! R_top [inout]  n by n upper triangular matrix 
  ! ldr_top [in]   Leading dimension of R_top
  ! R_bot [inout]  n by n upper triangular matrix 
  ! ldr_bot [in]   Leading dimension of R_bot
  ! tau [out]      Scaling factors for Householder reflectors
  ! work [out]     Workspace array (of length >= n)
  !
  subroutine d_factor_pair( n, R_top, ldr_top, R_bot, ldr_bot, tau, work )
    implicit none

    integer, value, intent(in)    :: n, ldr_top, ldr_bot
    real(8), intent(inout)        :: R_top(ldr_top,n), R_bot(ldr_bot,n)
    real(8), intent(out)          :: tau(n)
    ! "target" means it's legitimate to have pointers alias this
    ! array, internal to this function.
    real(8), intent(out), target  :: work(n)

    real(8), pointer  :: y(:)
    integer           :: j, k

    y => work(1:n)
    y = 0
    do k = 1, n-1
       ! Form the "sparse" Householder reflector, using DLARFP instead
       ! of DLARFG (need LAPACK 3.2, thanks Jason Riedy!) so that the
       ! diagonal elements of the R factor are nonnegative.  Length of
       ! this reflector is k+1, including the one top element (on the
       ! diagonal of R_top) and k bottom elements (above and including
       ! the diagonal of R_bot).
       call DLARFP( k + 1, R_top(k,k), R_bot(1,k), 1, tau(k) )

       ! $y^T := R_bot(1:k,k)^T R_bot(1:k, k+1:n)$, so $y := R_bot(1:k, k+1:n)^T R_bot(1:k, k)$.
       call DGEMV( 'T', k, n-k, 1.0d0, R_bot(1, k+1), ldr_bot, R_bot(1, k), 1, 0.0d0, y, 1 )

       ! $y^T := y^T + R_top(k, k+1:n)$
       do j = k+1, n
          y(j-k) = y(j-k) + R_top(k, j)
       end do

       ! Update R_top(k, k+1:n)
       do j = k+1, n
          R_top(k, j) = R_top(k, j) - tau(k) * y(j-k)
       end do

       ! Update R_bot(1:k, k+1:n)
       call DGER( k, n-k, -tau(k), R_bot(1,k), 1, y, 1, R_bot(1, k+1), ldr_bot )
    end do

    ! Compute the Householder reflector for the last column.  This
    ! last iteration doesn't require an update of the trailing matrix,
    ! because there is no trailing matrix left!
    call DLARFP( n + 1, R_top(n,n), R_bot(1,n), 1, tau(n) )
  end subroutine d_factor_pair


  subroutine s_factor_pair( n, R_top, ldr_top, R_bot, ldr_bot, tau, work )
    implicit none

    integer, value, intent(in)    :: n, ldr_top, ldr_bot
    real(4), intent(inout)        :: R_top(ldr_top,n), R_bot(ldr_bot,n)
    real(4), intent(out)          :: tau(n)
    ! "target" means it's legitimate to have pointers alias this
    ! array, internal to this function.
    real(4), intent(out), target  :: work(n)

    real(4), pointer  :: y(:)
    integer           :: j, k

    y => work(1:n)
    y = 0
    do k = 1, n-1
       ! Form the "sparse" Householder reflector, using SLARFP instead
       ! of SLARFG (need LAPACK >= 3.2, thanks Jason Riedy!) so that
       ! the diagonal elements of the R factor are nonnegative.
       ! Length of this reflector is k+1, including the one top
       ! element (on the diagonal of R_top) and k bottom elements
       ! (above and including the diagonal of R_bot).
       call SLARFP( k + 1, R_top(k,k), R_bot(1,k), 1, tau(k) )

       ! $y^T := R_bot(1:k,k)^T R_bot(1:k, k+1:n)$, so $y := R_bot(1:k, k+1:n)^T R_bot(1:k, k)$.
       call SGEMV( 'T', k, n-k, 1.0, R_bot(1, k+1), ldr_bot, R_bot(1, k), 1, 0.0, y, 1 )

       ! $y^T := y^T + R_top(k, k+1:n)$
       do j = k+1, n
          y(j-k) = y(j-k) + R_top(k, j)
       end do

       ! Update R_top(k, k+1:n)
       do j = k+1, n
          R_top(k, j) = R_top(k, j) - tau(k) * y(j-k)
       end do

       ! Update R_bot(1:k, k+1:n)
       call SGER( k, n-k, -tau(k), R_bot(1,k), 1, y, 1, R_bot(1, k+1), ldr_bot )
    end do

    ! Compute the Householder reflector for the last column.  This
    ! last iteration doesn't require an update of the trailing matrix,
    ! because there is no trailing matrix left!
    call SLARFP( n + 1, R_top(n,n), R_bot(1,n), 1, tau(n) )
  end subroutine s_factor_pair


  ! Apply Q factor (or Q^T) of the 2*ncols_Q by ncols_Q matrix 
  ! [R_top; R_bot] (stored in R_bot and tau) to the 2*ncols_Q by ncols_C 
  ! matrix [C_top; C_bot].  The two blocks C_top and C_bot may have 
  ! different leading dimensions (ldc_top resp. ldc_bot).
  !
  subroutine d_apply_pair( trans, ncols_C, ncols_Q, R_bot, ldr_bot, &
       tau, C_top, ldc_top, C_bot, ldc_bot, work )
    implicit none
    
    character, intent(in)        :: trans
    integer, intent(in), value   :: ncols_Q, ncols_C, ldr_bot, ldc_top, ldc_bot
    real(8), intent(in)          :: R_bot(ldr_bot,ncols_Q)
    real(8), intent(inout)       :: C_top(ldc_top,ncols_C), C_bot(ldc_bot,ncols_C)
    real(8), intent(in)          :: tau(ncols_Q)
    real(8), intent(out), target :: work(ncols_C)

    real(8), pointer :: y(:)
    integer          :: j, k

    y => work(1:ncols_C)
    y = 0
    if (trans == 'N' .or. trans == 'n') then
       do k = ncols_Q, 1, -1
          ! The k bottom elements of the length k+1 "sparse" Householder
          ! reflector are stored in R_bot(1:k,k).  The topmost element is
          ! implicit and is 1.0.  If we want to exploit its sparsity, we
          ! can't use DLARF to apply it to [C_top; C_bot]; we have to
          ! roll our own routine, which you see here.
          
          ! $y^T := R_bot(1:k,k)^T C_bot(1:k, 1:ncols_C)$, so $y := C_bot(1:k, 1:ncols_C)^T R_bot(1:k, k)$.
          call DGEMV( 'T', k, ncols_C, 1.0d0, C_bot, ldc_bot, R_bot(1, k), 1, 0.0d0, y, 1 )
          
          ! $y^T := y^T + C_top(k, 1:ncols_C)$
          do j = 1, ncols_C
             y(j) = y(j) + C_top(k, j)
          end do
          
          ! Update C_top(k, 1:ncols_C)
          do j = 1, ncols_C
             C_top(k, j) = C_top(k, j) - tau(k) * y(j)
          end do
          
          ! Update C_bot(1:k, 1:ncols_C)
          call DGER( k, ncols_C, -tau(k), R_bot(1,k), 1, y, 1, C_bot, ldc_bot )
       end do
    else
       do k = 1, ncols_Q
          ! The k bottom elements of the length k+1 "sparse" Householder
          ! reflector are stored in R_bot(1:k,k).  The topmost element is
          ! implicit and is 1.0.  If we want to exploit its sparsity, we
          ! can't use DLARF to apply it to [C_top; C_bot]; we have to
          ! roll our own routine, which you see here.
          
          ! $y^T := R_bot(1:k,k)^T C_bot(1:k, 1:ncols_C)$, so $y := C_bot(1:k, 1:ncols_C)^T R_bot(1:k, k)$.
          call DGEMV( 'T', k, ncols_C, 1.0d0, C_bot, ldc_bot, R_bot(1, k), 1, 0.0d0, y, 1 )
          
          ! $y^T := y^T + C_top(k, 1:ncols_C)$
          do j = 1, ncols_C
             y(j) = y(j) + C_top(k, j)
          end do
          
          ! Update C_top(k, 1:ncols_C)
          do j = 1, ncols_C
             C_top(k, j) = C_top(k, j) - tau(k) * y(j)
          end do
          
          ! Update C_bot(1:k, 1:ncols_C)
          call DGER( k, ncols_C, -tau(k), R_bot(1,k), 1, y, 1, C_bot(1, 1), ldc_bot )
       end do
    end if
  end subroutine d_apply_pair


  subroutine s_apply_pair( trans, ncols_C, ncols_Q, R_bot, ldr_bot, &
       tau, C_top, ldc_top, C_bot, ldc_bot, work )
    implicit none
    
    character, intent(in)        :: trans
    integer, intent(in), value   :: ncols_Q, ncols_C, ldr_bot, ldc_top, ldc_bot
    real(4), intent(in)          :: R_bot(ldr_bot,ncols_Q)
    real(4), intent(inout)       :: C_top(ldc_top,ncols_C), C_bot(ldc_bot,ncols_C)
    real(4), intent(in)          :: tau(ncols_Q)
    real(4), intent(out), target :: work(ncols_C)

    real(4), pointer :: y(:)
    integer          :: j, k

    y => work(1:ncols_C)
    y = 0
    if (trans == 'N' .or. trans == 'n') then
       do k = ncols_Q, 1, -1
          ! The k bottom elements of the length k+1 "sparse" Householder
          ! reflector are stored in R_bot(1:k,k).  The topmost element is
          ! implicit and is 1.0.  If we want to exploit its sparsity, we
          ! can't use DLARF to apply it to [C_top; C_bot]; we have to
          ! roll our own routine, which you see here.
          
          ! $y^T := R_bot(1:k,k)^T C_bot(1:k, 1:ncols_C)$, so $y := C_bot(1:k, 1:ncols_C)^T R_bot(1:k, k)$.
          call SGEMV( 'T', k, ncols_C, 1.0, C_bot, ldc_bot, R_bot(1, k), 1, 0.0, y, 1 )
          
          ! $y^T := y^T + C_top(k, 1:ncols_C)$
          do j = 1, ncols_C
             y(j) = y(j) + C_top(k, j)
          end do
          
          ! Update C_top(k, 1:ncols_C)
          do j = 1, ncols_C
             C_top(k, j) = C_top(k, j) - tau(k) * y(j)
          end do
          
          ! Update C_bot(1:k, 1:ncols_C)
          call SGER( k, ncols_C, -tau(k), R_bot(1,k), 1, y, 1, C_bot, ldc_bot )
       end do
    else
       do k = 1, ncols_Q
          ! The k bottom elements of the length k+1 "sparse" Householder
          ! reflector are stored in R_bot(1:k,k).  The topmost element is
          ! implicit and is 1.0.  If we want to exploit its sparsity, we
          ! can't use DLARF to apply it to [C_top; C_bot]; we have to
          ! roll our own routine, which you see here.
          
          ! $y^T := R_bot(1:k,k)^T C_bot(1:k, 1:ncols_C)$, so $y := C_bot(1:k, 1:ncols_C)^T R_bot(1:k, k)$.
          call SGEMV( 'T', k, ncols_C, 1.0, C_bot, ldc_bot, R_bot(1, k), 1, 0.0, y, 1 )
          
          ! $y^T := y^T + C_top(k, 1:ncols_C)$
          do j = 1, ncols_C
             y(j) = y(j) + C_top(k, j)
          end do
          
          ! Update C_top(k, 1:ncols_C)
          do j = 1, ncols_C
             C_top(k, j) = C_top(k, j) - tau(k) * y(j)
          end do
          
          ! Update C_bot(1:k, 1:ncols_C)
          call SGER( k, ncols_C, -tau(k), R_bot(1,k), 1, y, 1, C_bot(1, 1), ldc_bot )
       end do
    end if
  end subroutine s_apply_pair


end module TsqrCombine
