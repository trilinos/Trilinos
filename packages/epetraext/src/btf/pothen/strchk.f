      subroutine   strchk   ( nrows , ncols , colstr, rowidx, nhrows, 
     $                        nhcols, nsrows, rnto  , cnto  , colset,
     $                        rowset, output )

c     ==================================================================
c     ==================================================================
c     ====  strchk -- check that square block has nonzero diagonal  ====
c     ==================================================================
c     ==================================================================

c     ... for debugging purposes only 
c         created by john lewis, bcs, sept. 18, 1990

c     --------------
c     ... parameters
c     --------------

      integer       nrows , ncols , nhrows, nhcols, nsrows, output

      integer       colstr (ncols+1), rowidx (*),
     $              rnto (nrows), cnto (ncols), colset (ncols),
     $              rowset (nrows)

c     -------------------
c     ... local variables
c     -------------------

      integer       i, row, col, xi

      logical       match

c     ==================================================================

      do 300 i = 1, nsrows

         row = rnto (nhrows + i)
         col = cnto (nhcols + i)

         match = .false.

         do 200 xi = colstr (col), colstr (col+1) - 1
            match = match .or. (rowidx(xi) .eq. row)
 200     continue

         if  (.not. match  .or. (rowset(row) .ne. col)  .or.
     $                          (colset(col) .ne. row) )  then
            write (output, *) ' failure in matching, row, col, ',
     $                         'rowset(row), colset (col)',
     $                        row, col, rowset(row), colset (col)
         endif

 300  continue

      return

      end

