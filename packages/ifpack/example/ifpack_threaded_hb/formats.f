      subroutine coocsr(nrow,nnz,a,ir,jc,ao,jao,iao)
c----------------------------------------------------------------------- 
      real*8 a(*),ao(*),x
      integer ir(*),jc(*),jao(*),iao(*)
c-----------------------------------------------------------------------
c  Coordinate     to   Compressed Sparse Row 
c----------------------------------------------------------------------- 
c converts a matrix that is stored in coordinate format
c  a, ir, jc into a row general sparse ao, jao, iao format.
c
c on entry:
c--------- 
c nrow	= dimension of the matrix 
c nnz	= number of nonzero elements in matrix
c a,
c ir, 
c jc    = matrix in coordinate format. a(k), ir(k), jc(k) store the nnz
c         nonzero elements of the matrix with a(k) = actual real value of
c 	  the elements, ir(k) = its row number and jc(k) = its column 
c	  number. The order of the elements is arbitrary. 
c
c on return:
c----------- 
c ir 	is destroyed
c
c ao, jao, iao = matrix in general sparse matrix format with ao 
c 	continung the real values, jao containing the column indices, 
c	and iao being the pointer to the beginning of the row, 
c	in arrays ao, jao.
c
c Notes:
c------ This routine is NOT in place.  See coicsr
c
c------------------------------------------------------------------------
      do 1 k=1,nrow+1
         iao(k) = 0
 1    continue
c determine row-lengths.
      do 2 k=1, nnz
         iao(ir(k)) = iao(ir(k))+1
 2    continue
c starting position of each row..
      k = 1
      do 3 j=1,nrow+1
         k0 = iao(j)
         iao(j) = k
         k = k+k0
 3    continue
c go through the structure  once more. Fill in output matrix.
      do 4 k=1, nnz
         i = ir(k)
         j = jc(k)
         x = a(k)
         iad = iao(i)
         ao(iad) =  x
         jao(iad) = j
         iao(i) = iad+1
 4    continue
c shift back iao
      do 5 j=nrow,1,-1
         iao(j+1) = iao(j)
 5    continue
      iao(1) = 1
      return
c------------- end of coocsr ------------------------------------------- 
      end
c----------------------------------------------------------------------- 
      subroutine ssrcsr (nrow,a,ja,ia,nzmax,ao,jao,iao,indu,ierr)
      integer ia(nrow+1),iao(nrow+1),ja(*),jao(nzmax),indu(nrow+1)
      real*8 a(*),ao(nzmax)
c-----------------------------------------------------------------------
c Symmetric Sparse Row   to    (regular) Compressed Sparse Row
c----------------------------------------------------------------------- 
c this subroutine converts  a symmetric  matrix in which only the lower 
c part is  stored in compressed sparse row format, i.e.,
c a matrix stored in symmetric sparse format, into a fully stored matrix
c i.e., a matrix where both the lower and upper parts are stored in 
c compressed sparse row format. the algorithm is in place (i.e. result 
c may be overwritten onto the input matrix a, ja, ia ----- ). 
c the output matrix delivered by ssrcsr is such that each row starts with
c the elements of the lower part followed by those of the upper part.
c----------------------------------------------------------------------- 
c on entry:
c--------- 
c	
c nrow  = row dimension of inout matrix
c a, 
c ia, 
c ja    = matrix in compressed sparse row format. This is assumed to be
c         a lower triangular matrix. 
c
c nzmax	= size of arrays ao and jao. ssrcsr will abort if the storage 
c	   provided in a, ja is not sufficient to store A. See ierr. 
c	
c on return:
c----------
c ao, iao, 
c   jao = output matrix in compressed sparse row format. The resulting 
c         matrix is symmetric and is equal to A+A**T - D, if
c         A is the original lower triangular matrix. ao, jao, iao,
c         can be the same as a, ja, ia in the calling sequence.
c      
c indu  = integer array of length nrow+1. If the input matrix is such 
c         that the last element in each row is its diagonal element then
c         on return, indu will contain the pointers to the diagonal 
c         element in each row of the output matrix. Otherwise used as
c         work array.
c ierr  = integer. Serving as error message. If the length of the arrays
c         ao, jao exceeds nzmax, ierr returns the minimum value
c         needed for nzmax. otherwise ierr=0 (normal return).
c 
c----------------------------------------------------------------------- 
      ierr = 0
      do 1 i=1,nrow+1
         indu(i) = 0     
 1    continue
c     
c     compute  number of elements in each row of strict upper part. 
c     put result in indu(i+1)  for row i. 
c     
      do 3 i=1, nrow
         do 2 k=ia(i),ia(i+1)-1 
            j = ja(k)
            if (j .lt. i) indu(j+1) = indu(j+1)+1
 2       continue 
 3    continue
c-----------
c     find addresses of first elements of ouput matrix. result in indu
c-----------
      indu(1) = 1 
      do 4 i=1,nrow
         lenrow = ia(i+1)-ia(i)
         indu(i+1) = indu(i) + indu(i+1) + lenrow
 4    continue
c--------------------- enough storage in a, ja ? --------
      nnz = indu(nrow+1)-1 
      if (nnz .gt. nzmax) then
         ierr = nnz
         return
      endif
c
c     now copy lower part (backwards).
c     
      kosav = indu(nrow+1)
      do 6 i=nrow,1,-1
         klast = ia(i+1)-1
         kfirst = ia(i)
         iao(i+1) = kosav
         ko = indu(i) 
         kosav = ko
         do 5 k = kfirst, klast
            ao(ko) = a(k)
            jao(ko) = ja(k)
	    ko = ko+1
 5       continue
         indu(i) = ko 
 6    continue
      iao(1) = 1
c
c     now copy upper part. Go through the structure of ao, jao, iao
c     that has already been copied (lower part). indu(i) is the address
c     of the next free location in row i for ao, jao.
c     
      do 8 i=1,nrow
c     i-th row is now in ao, jao, iao structure -- lower half part
         do 9 k=iao(i), iao(i+1)-1 
            j = jao(k)
            if (j .ge. i)  goto 8
            ipos = indu(j)
            ao(ipos) = ao(k)
            jao(ipos) = i
            indu(j) = indu(j) + 1 
 9       continue
 8    continue
      return
c----- end of ssrcsr --------------------------------------------------- 
c----------------------------------------------------------------------- 
      end
c----------------------------------------------------------------------- 
      subroutine csrmsr (n,a,ja,ia,ao,jao,wk,iwk)
      real*8 a(*),ao(*),wk(n)
      integer ia(n+1),ja(*),jao(*),iwk(n+1)
c----------------------------------------------------------------------- 
c Compressed Sparse Row   to      Modified - Sparse Row 
c                                 Sparse row with separate main diagonal
c----------------------------------------------------------------------- 
c converts a general sparse matrix a, ja, ia into 
c a compressed matrix using a separated diagonal (referred to as
c the bell-labs format as it is used by bell labs semi conductor
c group. We refer to it here as the modified sparse row format.
c Note: this has been coded in such a way that one can overwrite
c the output matrix onto the input matrix if desired by a call of
c the form 
c
c     call csrmsr (n, a, ja, ia, a, ja, wk,iwk)
c
c In case ao, jao, are different from a, ja, then one can
c use ao, jao as the work arrays in the calling sequence:
c
c     call csrmsr (n, a, ja, ia, ao, jao, ao,jao)
c
c----------------------------------------------------------------------- 
c
c on entry :
c---------
c a, ja, ia = matrix in csr format. note that the 
c	     algorithm is in place: ao, jao can be the same
c            as a, ja, in which case it will be overwritten on it
c            upon return.
c	 
c on return :
c-----------
c
c ao, jao  = sparse matrix in modified sparse row storage format:
c	   +  ao(1:n) contains the diagonal of the matrix. 
c	   +  ao(n+2:nnz) contains the nondiagonal elements of the
c             matrix, stored rowwise.
c	   +  jao(n+2:nnz) : their column indices
c	   +  jao(1:n+1) contains the pointer array for the nondiagonal
c             elements in ao(n+1:nnz) and jao(n+2:nnz).
c             i.e., for i .le. n+1 jao(i) points to beginning of row i 
c	      in arrays ao, jao.
c	       here nnz = number of nonzero elements+1 
c work arrays:
c------------
c wk	= real work array of length n
c iwk   = integer work array of length n+1
c
c notes: 
c------- 
c        Algorithm is in place.  i.e. both:
c
c          call csrmsr (n, a, ja, ia, ao, jao, ao,jao)
c          (in which  ao, jao, are different from a, ja)
c           and
c          call csrmsr (n, a, ja, ia, a, ja, wk,iwk) 
c          (in which  wk, jwk, are different from a, ja)
c        are OK.
c--------
c coded by Y. Saad Sep. 1989. Rechecked Feb 27, 1990.
c-----------------------------------------------------------------------
      icount = 0
c
c store away diagonal elements and count nonzero diagonal elements.
c
      do 1 i=1,n
         wk(i) = 0.0d0
         iwk(i+1) = ia(i+1)-ia(i)
         do 2 k=ia(i),ia(i+1)-1
            if (ja(k) .eq. i) then
               wk(i) = a(k)
               icount = icount + 1 
               iwk(i+1) = iwk(i+1)-1
            endif
 2       continue
 1    continue
c     
c compute total length
c     
      iptr = n + ia(n+1) - icount
c     
c     copy backwards (to avoid collisions)
c     
      do 500 ii=n,1,-1
         do 100 k=ia(ii+1)-1,ia(ii),-1
            j = ja(k)
            if (j .ne. ii) then
               ao(iptr) = a(k)
               jao(iptr) = j 
               iptr = iptr-1
            endif
 100     continue
 500  continue
c
c compute pointer values and copy wk(*)
c
      jao(1) = n+2
      do 600 i=1,n
         ao(i) = wk(i) 
         jao(i+1) = jao(i)+iwk(i+1)
 600  continue
      return	
c------------ end of subroutine csrmsr ---------------------------------
c----------------------------------------------------------------------- 
      end
c----------------------------------------------------------------------- 
      subroutine csrcsc (n,job,ipos,a,ja,ia,ao,jao,iao)
      integer ia(n+1),iao(n+1),ja(*),jao(*)
      real*8  a(*),ao(*)
c-----------------------------------------------------------------------
c Compressed Sparse Row     to      Compressed Sparse Column
c
c (transposition operation)   Not in place. 
c----------------------------------------------------------------------- 
c -- not in place --
c this subroutine transposes a matrix stored in a, ja, ia format.
c ---------------
c on entry:
c----------
c n	= dimension of A.
c job	= integer to indicate whether to fill the values (job.eq.1) of the
c         matrix ao or only the pattern., i.e.,ia, and ja (job .ne.1)
c
c ipos  = starting position in ao, jao of the transposed matrix.
c         the iao array takes this into account (thus iao(1) is set to ipos.)
c         Note: this may be useful if one needs to append the data structure
c         of the transpose to that of A. In this case use for example
c                call csrcsc (n,1,n+2,a,ja,ia,a,ja,ia(n+2)) 
c	  for any other normal usage, enter ipos=1.
c a	= real array of length nnz (nnz=number of nonzero elements in input 
c         matrix) containing the nonzero elements.
c ja	= integer array of length nnz containing the column positions
c 	  of the corresponding elements in a.
c ia	= integer of size n+1. ia(k) contains the position in a, ja of
c	  the beginning of the k-th row.
c
c on return:
c ---------- 
c output arguments:
c ao	= real array of size nzz containing the "a" part of the transpose
c jao	= integer array of size nnz containing the column indices.
c iao	= integer array of size n+1 containing the "ia" index array of
c	  the transpose. 
c
c----------------------------------------------------------------------- 
c----------------- compute lengths of rows of transp(A) ----------------
      do 1 i=1,n+1
         iao(i) = 0
 1    continue
      do 3 i=1, n
         do 2 k=ia(i), ia(i+1)-1 
            j = ja(k)+1
            iao(j) = iao(j)+1
 2       continue 
 3    continue
c---------- compute pointers from lengths ------------------------------
      iao(1) = ipos 
      do 4 i=1,n
         iao(i+1) = iao(i) + iao(i+1)
 4    continue
c--------------- now do the actual copying ----------------------------- 
      do 6 i=1,n
         do 62 k=ia(i),ia(i+1)-1 
            j = ja(k) 
            next = iao(j)
            if (job .eq. 1)  ao(next) = a(k)
            jao(next) = i
            iao(j) = next+1
 62      continue
 6    continue
c-------------------------- reshift iao and leave ---------------------- 
      do 7 i=n,1,-1
         iao(i+1) = iao(i)
 7    continue
      iao(1) = ipos
c--------------- end of csrcsc ----------------------------------------- 
c-----------------------------------------------------------------------
      end
c-----------------------------------------------------------------------
