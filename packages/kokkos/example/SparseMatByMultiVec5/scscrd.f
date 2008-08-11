      subroutine scscrd ( nmax, nzmax, job, iounit, a, ja, ia, rhs,
     &                    rhsptr, rhsind, nrhs, xexact, sguess, 
     &                    rhstyp, nrow, ncol,
     &                    nnz, title, key, type, ierr)
c-----------------------------------------------------------------------
c this subroutine reads  a boeing/harwell matrix. handles right hand 
c sides in full format only (no sparse right hand sides).
c-----------------------------------------------------------------------
c on entry:
c---------
c nmax 	 =  max dimension  allowed for matrix. The array ia should 
c	    be of length at least ncol+1 (see below) if job.gt.0
c nzmax	 = max number of nonzeros elements allowed. the arrays a, 
c          and ja should be of length equal to nnz (see below) if these
c          arrays are to be read (see job).
c          
c job	 = integer to indicate what is to be read. (note: job is an
c          input and output parameter, it can be modified on return)
c          job = 0    read the values of ncol, nrow, nnz title, key,
c                     type and return. matrix is not read and arrays
c                     a, ja, ia, rhs are not touched.
c          job = 1    read srtucture only, i.e., the arrays ja and ia.
c          job = 2    read matrix including values, i.e., a, ja, ia
c          job = 3    read matrix and right hand sides: a,ja,ia,rhs.
c		      rhs may contain initial guesses and exact 
c                     solutions appended to the actual right hand sides.
c		      this will be indicated by the output parameter
c                     guesol [see below]. 
c                     
c nrhs   = integer. nrhs is an input as well as ouput parameter.
c          at input nrhs contains the total length of the array rhs.
c          See also ierr and nrhs in output parameters.
c
c iounit = logical unit number where to read the matrix from.
c
c on return:
c---------- 
c job    = on return job may be modified to the highest job it could
c          do: if job=2 on entry but no matrix values are available it
c          is reset to job=1 on return. Similarly of job=3 but no rhs 
c          is provided then it is rest to job=2 or job=1 depending on 
c          whether or not matrix values are provided.
c          Note that no error message is triggered (i.e. ierr = 0 
c          on return in these cases. It is therefore important to
c          compare the values of job on entry and return ).
c
c a	 = the a matrix in the a, ia, ja (column) storage format
c ja 	 = row number of element a(i,j) in array a.
c ia     = pointer  array. ia(i) points to the beginning of column i.
c
c rhs    = real array of size nrow + 1 if available (see job)
c
c nrhs   = integer containing the number of right-hand sides found
c          each right hand side may be accompanied with an intial guess
c          and also the exact solution.
c
c rhstyp = a 3-character string indicating the form of the rhs
c          (1st-character), whether an initial guess 
c          (2-nd character) and / or the exact solution (3-rd
c          character) is provided with the right hand side.
c       if the first character of rhstyp is 'F' then the rhs is in
c          full storage.  if it is 'M' then the rhs is stored in the
c          same format as the matrix.
c	   if the second character of rhstyp is 'G' it means that an
c          an intial guess is provided for each right-hand side.
c          These are put sguess.
c	   if the third character of rhstyp is 'X' it means that an
c          exact solution is provided for each right-hand side.
c          These are in xexact.
c
c nrow   = number of rows in matrix
c ncol	 = number of columns in matrix 
c nnz	 = number of nonzero elements in A. This info is returned
c          even if there is not enough space in a, ja, ia, in order
c          to determine the minimum storage needed. 
c
c title  = character*72 = title of matrix test ( character a*72). 
c key    = character*8  = key of matrix 
c type   = charatcer*3  = type of matrix.
c          for meaning of title, key and type refer to documentation 
c          Harwell/Boeing matrices.
c
c ierr   = integer used for error messages 
c         * ierr  =  0 means that  the matrix has been read normally. 
c         * ierr  =  1 means that  the array matrix could no be read 
c         because ncol+1 .gt. nmax
c         * ierr  =  2 means that  the array matrix could no be read 
c         because nnz .gt. nzmax 
c         * ierr  =  3 means that  the array matrix could no be read 
c         because both (ncol+1 .gt. nmax) and  (nnz .gt. nzmax )
c         * ierr  =  4 means that  the right hand side (s) initial 
c         guesse (s) and exact solution (s)   could  not be
c         read because they are stored in sparse format (not handled
c         by this routine ...) 
c         * ierr  =  5 means that the right-hand-sides, initial guesses
c         and exact solutions could not be read because the length of 
c         rhs as specified by the input value of nrhs is not 
c         insufficient to store them. The rest of the matrix may have
c         been read normally.
c 
c Notes:
c-------
c 1) The file inout must be open (and possibly rewound if necessary)
c    prior to calling scscrd.
c 2) Refer to the documentation on the Harwell-Boeing formats
c    for details on the format assumed by scscrd.
c    We summarize the format here for convenience.
c  
c    a) all lines in inout are assumed to be 80 character long.
c    b) the file consists of a header followed by the block of the 
c       column start pointers followed by the block of the
c       row indices, followed by the block of the real values and
c       finally the numerical values of the right-hand-side if a 
c       right hand side is supplied. 
c    c) the file starts by a header which contains four lines if no
c       right hand side is supplied and five lines otherwise.
c       * first line contains the title (72 characters long) followed by
c         the 8-character identifier (name of the matrix, called key)
c        [ A72,A8 ]
c       * second line contains the number of lines for each
c         of the following data blocks (4 of them) and the total number 
c         of lines excluding the header.
c        [5i4]
c       * the third line contains a three character string identifying
c         the type of matrices as they are referenced in the Harwell
c         Boeing documentation [e.g., rua, rsa,..] and the number of
c         rows, columns, nonzero entries.
c         [A3,11X,4I14]
c       * The fourth line contains the variable fortran format
c         for the following data blocks.
c         [2A16,2A20] 
c       * The fifth line is only present if right-hand-sides are 
c         supplied. It consists of three one character-strings containing
c         the storage format for the right-hand-sides 
c         ('F'= full,'M'=sparse=same as matrix), an initial guess 
c         indicator ('G' for yes), an exact solution indicator 
c         ('X' for yes), followed by the number of right-hand-sides
c         and then the number of row indices. 
c         [A3,11X,2I14] 
c     d) The three following blocks follow the header as described 
c        above.
c     e) In case the right hand-side are in sparse formats then 
c        the fourth block uses the same storage format as for the matrix
c        to describe the NRHS right hand sides provided, with a column
c        being replaced by a right hand side.
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      character title*72, key*8, type*3, ptrfmt*16, indfmt*16,
     1       valfmt*20, rhsfmt*20, rhstyp*3
      integer totcrd, ptrcrd, indcrd, valcrd, rhscrd, nrow, ncol,
     1     nnz, neltvl, nrhs, nrhsix, nmax, nzmax
      integer ia (nmax+1), ja (nzmax), rhsind(*), rhsptr(*)
      real a(nzmax), rhs(*) , sguess(*), xexact(*)
c-----------------------------------------------------------------------
c
      read (iounit,10) title, key, totcrd, ptrcrd, indcrd, valcrd, 
     1     rhscrd, type, nrow, ncol, nnz, neltvl, ptrfmt, indfmt, 
     2     valfmt, rhsfmt
 10   format (a72, a8 / 5i14 / a3, 11x, 4i14 / 2a16, 2a20)
c
      if (rhscrd .gt. 0) read (iounit,11) rhstyp, nrhs ,nrhsix
 11   format (a3,11x,2i14)
c
c anything else to read ?
c
      if (job .le. 0) return
      ierr = 0
c     ---- check whether matrix is readable ------ 
      n = ncol
      if (ncol+1 .gt. nmax) ierr = 1
      if (nnz .gt. nzmax) ierr = ierr + 2
      if (ierr .ne. 0) return
c     ---- read pointer and row numbers ---------- 
      read (iounit,ptrfmt) (ia (i), i = 1, n+1)
      read (iounit,indfmt) (ja (i), i = 1, nnz)
c     --- reading values of matrix if required....
      if (job .le. 1)  return
c     --- and if available ----------------------- 
      if (valcrd .le. 0) then
	 job = 1
	 return
      endif
      read (iounit,valfmt) (a(i), i = 1, nnz)
c     --- reading rhs if required ---------------- 
      if (job .le. 2)  return
c     --- and if available ----------------------- 
      if ( rhscrd .le. 0) then
	 job = 2
	 return
      endif
c     
c     --- read right-hand-side.-------------------- 
c     
      if (nrhs.gt.0) then
        if (rhstyp(1:1) .eq. 'F') then 
c
c read right-hand-sides
c
          len = nrow*nrhs
          read(iounit,rhsfmt) (rhs(i), i = 1, len)
        else
          if (type(3:3).eq.'A') then
            read(iounit,ptrfmt)(rhsptr(i),i=1,nrhs+1)
            read(iounit,indfmt)(rhsind(i),i=1,nrhsix)
            read(iounit,rhsfmt)(rhs   (i),i=1,nrhsix)
          else
            write(*,*)'Cannot read unassembled matrices'
            STOP
          endif
        endif
      endif
c
c read initial guesses if available
c
      if (rhstyp(2:2) .eq. 'G') then
        len = nrow*nrhs
        read(iounit,rhsfmt) (sguess(i), i = 1, len)
      endif
c     
c read exact solutions if available
c
      if (rhstyp(3:3) .eq. 'X') then
        len = nrow*nrhs
        read(iounit,rhsfmt) (xexact(i), i = 1, len)
      endif
c     
      return
c--------- end of scscrd -----------------------------------------------
c----------------------------------------------------------------------- 
      end
