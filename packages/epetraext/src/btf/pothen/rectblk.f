      subroutine   rectblk   ( nrows , ncols , marked, unmrkd, colstr,
     $                         rowidx, colset, rowset, prevcl, tryrow, 
     $                         colmrk, rowmrk, nhrows, nhcols )

c     ==================================================================
c     ==================================================================
c     ====  rectblk -- find rectangular portion of matrix by        ====
c     ====             depth-first search                           ====
c     ==================================================================
c     ==================================================================

c     original -- alex pothen and chin-ju fan, penn state, 1988
c     bcs modifications, john lewis, sept. 1990

c     use a depth-first serch to find all the rows and columns, which
c     can be reached via alternating paths beginning from all the
c     unmatched columns.  comments and names describe use of code
c     for finding the 'horizontal' block.  the same code is used
c     to find the vertical block by performing exactly the same
c     operations on the transpose of the matrix.
c
c     input variables:
c
c         nrows    -- number of rows
c         ncols    -- number of columns
c         marked   -- value to store in marker vectors to indicate
c                     that row/column has been reached and is
c                     therefore in the horizontal block
c         unmrkd   -- initial value of marker vectors, indicating
c                     that row or column is free to be chosen
c         colstr, 
c         rowidx   -- adjacency structure of graph
c         colset   -- maximum matching for columns
c         rowset   -- maximum matching for rows
c
c    output variables:
c
c         nhrows  -- number of rows in horizontal block
c         nhcols  -- number of columns in horizontal block 
c         rowmrk, 
c         colmrk  -- row and column marker vectors.  
c                    = unmrkd --> row/column is in neither the
c                                  horizontal or vertical block yet
c                    = marked --> row/column has been reached via
c                                 search in this routine and lies
c                                 in the horizontal block
c                    = neither --> row/column is not free for use.
c                                  it was found to lie in another
c                                  block.
c                                  
c    working variables:
c
c         tryrow -- tryrow (col) is a pointer into rowidx to the
c                   next row to be explored from col 'col' in
c                   the search.
c         prevcl -- pointer toward the root of the search from 
c                   column to column.
c
c     ==================================================================

c     --------------
c     ... parameters
c     --------------

      integer        nrows, ncols, marked, unmrkd, nhcols, nhrows

      integer        colstr (nrows+1), rowidx (*), rowset (nrows),
     $               colset (ncols)

      integer        prevcl (ncols), tryrow (ncols), colmrk (ncols),
     $               rowmrk (nrows)

c     -------------------
c     ... local variables
c     -------------------

      integer        col, fromc, nextcl, nextrw, p, row, xrow

c     ==================================================================

      nhcols = 0
      nhrows = 0

      do 300 p = 1, ncols

c        -----------------------------------------------------------
c        ... find an unmatched column to start the alternating path.
c        -----------------------------------------------------------

         if  ( colset (p) .eq. 0 )  then

            fromc = p

c           ---------------------------------------------
c           ... path starts from unmatched column "fromc"
c               put fromc into horizontal set "hc"
c               indicate fromc is the root of the path.
c           ---------------------------------------------

            nhcols         = nhcols + 1
            colmrk (fromc) = marked
            tryrow (fromc) = colstr (fromc)
            prevcl (fromc) = 0
            col            =  fromc

c           ------------------------------------------------------
c           ... main depth-first search loop begins here.
c               Each time through take a step forward if possible
c               or backtrack if not. quit when we backtrack to the
c               beginning of the search.
c           ------------------------------------------------------
c     
c           ------------------------------------------------
c           ... look for a forward step from column 'col' to
c               an unmarked row.
c           ------------------------------------------------

 100        nextrw = tryrow (col)
            do 200 xrow = nextrw, colstr (col + 1) - 1

               if  ( rowmrk (rowidx (xrow)) .eq. unmrkd )  then

c                 ---------------------------------------------------
c                 ... take a double forward step from 'col' to 'row'
c                     and then via matching edge from 'row' to column
c                     'nextcl'.  ('row' must be matched since 
c                     otherwise we have found an augmenting path
c                     and the maximum matching wasn't matching.)
c                 ---------------------------------------------------

                  tryrow (col) = xrow + 1
                  row          = rowidx (xrow)
                  rowmrk (row) = marked
                  nhrows       = nhrows + 1

                  nextcl       = rowset (row)
                  if  ( nextcl .eq. 0 )  then
                     write (6, 60000)
60000                format (' max matching is wrong -- augmenting ',
     $                         'path found')
                     stop
                  endif
      
                  nhcols          = nhcols + 1
                  colmrk (nextcl) = marked
                  prevcl (nextcl) = col
                  tryrow (nextcl) = colstr (nextcl)
                  col             = nextcl
                  go to 100
               endif
      
 200        continue

c           ------------------------------------------------
c           ... no forward step: backtrack.  if we backtrack
c               all the way, we have completed all searchs
c               beginning at column 'p'.
c           ------------------------------------------------

            col = prevcl (col)
            if  ( col .ne. 0 )  then
               go to 100
            endif

         endif

  300 continue

      return

      end

