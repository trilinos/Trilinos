      subroutine   genbtf   ( nrows , ncols , colstr, rowidx, rowstr,
     $                        colidx, w     , rnto  , cnto  , nhrows,
     $                        nhcols, hrzcmp, nsrows, sqcmpn, nvrows,
     $                        nvcols, vrtcmp, rcmstr, ccmstr, msglvl,
     $                        output )

c     ==================================================================
c     ==================================================================
c     ====  genbtf -- find the block triangular form (dulmadge-     ====
c     ====            mendelson decomposition) of a general         ====
c     ====            rectangular sparse matrix                     ====
c     ==================================================================
c     ==================================================================
c  
c     created        sept. 14, 1990 (jgl)
c     last modified  oct. 4, 1990 (jgl)

c     algorithm by alex pothen and chin-ju fan
c     this code based on code from alex pothen, penn state university

c     ... input variables
c     -------------------
c
c        nrows  -- number of rows in matrix
c        ncols  -- number of columns in matrix
c        colstr, rowidx
c               -- adjacency structure of matrix, where each
c                  column of matrix is stored contiguously
c                  (column-wise representation)
c        rowstr, colidx
c               -- adjacency structure of matrix, where each
c                  row of matrix is stored contiguously
c                  (row-wise representation)
c                  (yes, two copies of the matrix)
c
c     ... temporary storage
c     ---------------------
c
c        w      -- integer array of length 5*nrows + 5*ncols
c
c     ... output variables:
c     ---------------------
c
c        rnto   -- the new to old permutation array for the rows
c        cotn   -- the old to new permutation array for the columns
c        nhrows, nhcols, hrzcmp 
c               -- number of rows, columns and connected components
c                  in the horizontal (underdetermined) block
c        nsrows, sqcmpn 
c               -- number of rows (and columns) and strong components
c                  in the square (exactly determined) block
c        nvrows, nvcols, vrtcmp 
c               -- number of rows, columns and connected components
c                  in the vertical (overdetermined) block
c        rcmstr -- index of first row in a diagonal block
c                  (component starting row)
c                  where
c                      (rcmstr(1), ..., rcmstr(hrzcmp)) give the 
c                           indices for the components in the
c                            horizontal block
c                      (rcmstr(hrzcmp+1), ..., rcmstr(hrzcmp+sqcmpn))
c                           give the indices for the components in the
c                           square block
c                      (rcmstr(hrzcmp+sqcmpn+1), ..., 
c                       rcmstr(hrzcmp+sqcmpn+vrtcmp)) give the indices
c                           for the components in the vertical block
c                       rcmstr(hrzcmp+sqcmpn+vrtcmp+1) is equal to
c                           nrows+1 for convenience
c        ccmstr -- index of first column in a diagonal block
c                  (component starting column)
c                  where
c                      (ccmstr(1), ..., ccmstr(hrzcmp)) give the 
c                           indices for the components in the
c                            horizontal block
c                      (ccmstr(hrzcmp+1), ..., ccmstr(hrzcmp+sqcmpn))
c                           give the indices for the components in the
c                           square block, making this block itself
c                           in block lower triangular form
c                      (ccmstr(hrzcmp+sqcmpn+1), ..., 
c                       ccmstr(hrzcmp+sqcmpn+vrtcmp)) give the indices
c                           for the components in the vertical block
c                       ccmstr(hrzcmp+sqcmpn+vrtcmp+1) is equal to
c                           ncols+1 for convenience
c
c                  note -- if the matrix has entirely empty rows,
c                          these rows will be placed in the vertical
c                          block, each as a component with one row
c                          and zero columns.  similarly, entirely
c                          empty columns will appear in the horizontal
c                          block, each as a component with no rows and
c                          one column.
c
c        msglvl -- message level
c                  = 0 -- no output
c                  = 1 -- timing and summary output
c                  = 2 -- adds final permutation vectors
c                  = 3 -- adds intermediate copies of vectros as
c                         debugging output
c     
c        output -- fortran unit number for printable output
c
c     efficiency note:
c     ----------------

c        although it is not required by this code that the number
c        of rows be larger than the number of columns, the first
c        phase (the matching) will be faster in this case.  thus,
c        in cases where the number of columns is substantially
c        larger than the number of rows, it will probably be more
c        efficient to apply this algorithm to the transpose of the
c        matrix.  since the matrix is required with both row and
c        column representations, applying the algorithm to the
c        transposed matrix can be achieved simply by interchanging
c        appropriate parameters in the call to  genbtf.
c
c     ==================================================================

c     --------------
c     ... parameters
c     --------------

      integer         nrows , ncols , nhrows, nhcols, hrzcmp, nsrows,
     $                sqcmpn, nvrows, nvcols, vrtcmp, msglvl, output

      integer         colstr (ncols+1), rowidx (*),
     $                rowstr (nrows+1), colidx (*),
     $                w      (*)      ,
     $                cnto   (ncols)  , rnto   (nrows),
     $                rcmstr (nrows+1), ccmstr (ncols+1)

c     -------------------
c     ... local variables
c     -------------------

      integer         cmk   , cmbase, cnbase, cst   , cw1   , cw2   ,
     $                cw3   , i     , hindex, ncompn, nscols, rmk   , 
     $                rnbase, rst   , rw1   , rw2   , rw3   , sqindx,
     $                vindex

      real            timeh , timem , times , timev , tmstrt

      real            tarray (2)

      real            etime

c     ==================================================================

c     --------------
c     ... initialize
c     --------------

      vindex = -1
      sqindx = -2
      hindex = -3

      cmk = 1
      cst = cmk + ncols
      rmk = cst + ncols
      rst = rmk + nrows
      rw1 = rst + nrows
      rw2 = rw1 + nrows
      rw3 = rw2 + nrows
      cw1 = rw3 + nrows
      cw2 = cw1 + ncols
      cw3 = cw2 + ncols
      call izero ( cw3 + ncols - 1, w, 1 )

c     ---------------------------------------
c     ... algorithm consists of three stages:
c         1.  find a maximum matching
c         2.  find a coarse decomposition
c         3.  find a fine decomposition
c     ---------------------------------------

c     -----------------------------
c     ... find the maximum matching
c     -----------------------------
      
      if  ( msglvl .ge. 1 )  then
         tmstrt = etime ( tarray )
      endif

      call maxmatch ( nrows , ncols , colstr, rowidx, w(cw1), w(cmk), 
     $                w(rw2), w(cw2), w(cw3), w(rst), w(cst) )

      do 100 i = 1, nrows
         w (rmk + i - 1) = sqindx
 100  continue

      do 200 i = 1, ncols
         w (cmk + i - 1) = sqindx
 200  continue

      if  ( msglvl .ge. 1 )  then
         timem = etime ( tarray ) - tmstrt
         if  ( msglvl .ge. 3 )  then
            call prtivs ( 'rowset', nrows, w(rst), output )
            call prtivs ( 'colset', ncols, w(cst), output )
         endif
      endif

c     ------------------------------------------------------------
c     ... coarse partitioning -- divide the graph into three parts
c     ------------------------------------------------------------

c     --------------------------------------------------------
c     ... depth-first search from unmatched columns to get the
c         horizontal submatrix
c     --------------------------------------------------------

      if  ( msglvl .ge. 1 )  then
         tmstrt = etime ( tarray )
      endif

      call rectblk ( nrows , ncols , hindex, sqindx, colstr, rowidx, 
     $               w(cst), w(rst), w(cw1), w(cw2), w(cmk), w(rmk), 
     $               nhrows, nhcols ) 

      if  ( msglvl .ge. 1 )  then
         timeh = etime ( tarray ) - tmstrt
         if  ( msglvl .ge. 3 )  then
            write ( output, * ) '0nhrows, nhcols', nhrows, nhcols
         endif
      endif

c     -----------------------------------------------------
c     ... depth-first search from unmatched rows to get the
c         vertical submatrix
c     -----------------------------------------------------

      if  ( msglvl .ge. 1 )  then
         tmstrt = etime ( tarray )
      endif

      tmstrt = etime ( tarray )

      call rectblk ( ncols , nrows , vindex, sqindx, rowstr, colidx, 
     $               w(rst), w(cst), w(rw1), w(rw2), w(rmk), w(cmk), 
     $               nvcols, nvrows )  

      if  ( msglvl .ge. 1 )  then
         timev = etime ( tarray ) - tmstrt
         if  ( msglvl .ge. 3 )  then
            write ( output, * ) '0nvrows, nvcols', nvrows, nvcols
         endif
      endif

c     ----------------------------------------
c     ... the square submatrix is what is left
c     ----------------------------------------

      nscols = ncols - nhcols - nvcols
      nsrows = nrows - nhrows - nvrows

      if  ( msglvl .ge. 1 )  then
         call corsum ( timem , timeh , timev , nhrows, nhcols, nsrows, 
     $                 nscols, nvrows, nvcols, output )
      endif

c     ----------------------------------------------
c     ... begin the fine partitioning and create the
c         new to old permutation vectors
c     ----------------------------------------------

c     ---------------------------------------------------------
c     ... find connected components in the horizontal submatrix 
c     ---------------------------------------------------------

      if  ( nhcols .gt. 0 )  then

         if  ( msglvl .ge. 1 )  then
            tmstrt = etime ( tarray )
         endif

         cmbase  = 0
         rnbase  = 0
         cnbase  = 0

         call concmp ( cmbase, cnbase, rnbase, hindex, ncols , nrows ,
     $                 nhcols, nhrows, colstr, rowidx, rowstr, colidx,
     $                 w(rw1), w(cw1), w(cw2), w(rw2), w(rw3), w(cw3), 
     $                 w(rmk), w(cmk), rcmstr, ccmstr, rnto  , cnto  , 
     $                 hrzcmp )
  
         if  ( msglvl .ge. 1 )  then
            timeh = etime ( tarray ) - tmstrt
            if  ( msglvl .ge. 3 )  then
               write ( output, * ) '0hrzcmp', hrzcmp
               call prtivs ( 'rcmstr', hrzcmp + 1, rcmstr, output )
               call prtivs ( 'ccmstr', hrzcmp + 1, ccmstr, output )
               call prtivs ( 'rnto', nrows, rnto, output )
               call prtivs ( 'cnto', ncols, cnto, output )
            endif
         endif

      else

         hrzcmp = 0
         timeh  = 0.0

      endif

      if  ( nsrows .gt. 0 ) then

         if  ( msglvl .ge. 1 )  then
            tmstrt = etime ( tarray )
         endif

c        -----------------------------------------------------------
c        ... find strongly connected components in square submatrix,
c            putting this block into block lower triangular form.
c        -----------------------------------------------------------

         call mmc13e ( nrows , ncols , nhcols, nhrows, nsrows, sqindx,
     $                 hrzcmp, rowstr, colidx, w(cst), w(rw1), w(rw2),
     $                 w(cw1), w(cw2), w(cmk), ccmstr, rcmstr, cnto  ,  
     $                 rnto  , sqcmpn )

         if  ( msglvl .ge. 1 )  then
            call strchk ( nrows , ncols , colstr, rowidx, nhrows,
     $                    nhcols, nsrows, rnto  , cnto  , w(cst),
     $                    w(rst), output )

         endif

         if  ( msglvl .ge. 1 )  then
            times = etime ( tarray ) - tmstrt
            if  ( msglvl .ge. 3 )  then
               ncompn = hrzcmp + sqcmpn + 1
               write ( output, * ) '0sqcmpn', sqcmpn
               call prtivs ( 'rcmstr', ncompn, rcmstr, output )
               call prtivs ( 'ccmstr', ncompn, ccmstr, output )
               call prtivs ( 'rnto', nrows, rnto, output )
               call prtivs ( 'cnto', ncols, cnto, output )
            endif
         endif

      else
         
         sqcmpn = 0
         times  = 0.0

      endif

      if  ( nvrows .gt. 0 )  then

         cmbase = hrzcmp + sqcmpn
         rnbase = nhrows + nscols
         cnbase = nhcols + nscols

c        ---------------------------------------------------
c        ... find connected components in vertical submatrix
c        ---------------------------------------------------

         if  ( msglvl .ge. 1 )  then
            tmstrt = etime ( tarray )
         endif

         call concmp ( cmbase, rnbase, cnbase, vindex, nrows , ncols ,
     $                 nvrows, nvcols, rowstr, colidx, colstr, rowidx,
     $                 w(cw1), w(rw1), w(rw2), w(cw2), w(cw3), w(rw3), 
     $                 w(cmk), w(rmk), ccmstr, rcmstr, cnto  , rnto  , 
     $                 vrtcmp )

         if  ( msglvl .ge. 1 )  then

            timev = etime ( tarray ) - tmstrt

            if  ( msglvl .ge. 2 )  then
               call prtivs ( 'rnto', nrows, rnto, output )
               call prtivs ( 'cnto', ncols, cnto, output )

               if  ( msglvl .ge. 3 )  then
                  ncompn = hrzcmp + sqcmpn + vrtcmp + 1
                  write ( output, * ) '0vrtcmp', vrtcmp
                  call prtivs ( 'rcmstr', ncompn, rcmstr, output )
                  call prtivs ( 'ccmstr', ncompn, ccmstr, output )
               endif

            endif
         endif

      else
         
         vrtcmp = 0
         timev  = 0.0

      endif

      if  ( msglvl .ge. 1 )  then
         call finsum ( timeh , times , timev , hrzcmp, sqcmpn, 
     $                 vrtcmp, ccmstr, rcmstr, output )
      endif

      return

      end
