      subroutine   mmc13e   ( nrows , ncols , nhcols, nhrows, nscols,
     $                        sqindx, hrzcmp, rowstr, colind, colset, 
     $                        trycol, cbegin, lowlnk, prev  , colmrk, 
     $                        ccmstr, rcmstr, cnto  , rnto  , sqcmpn )

c     ==================================================================
c     ==================================================================
c     ====  mmc13e -- lower block triangular form of square matrix  ====
c     ==================================================================
c     ==================================================================

c     mmc13e :   modified from harwell mc13e by alex pothen and
c                chin-ju fan
c     bcs modifications, john lewis, sept. 1990

c     finds the lower block triangular form of the square submatrix
c     in the general block triangular form.  the square submatrix
c     consists entirely of matched rows and columns.  therefore,
c     with each row matched to its matching column, the submatrix
c     has a nonzero diagonal, as required by duff's algorithm.
c
c     from a graph-theoretic standard, this is the same as considering
c     the subgraph induced by sr and sc, if non-matching edges
c     are directed from rows to columns, and matching edges are shrunk 
c     into single vertices, the resulting directed graph has strongly 
c     connected components.
c
c     mmc13e uses Tarjan's algorithm to find the strongly connected
c     components by depth-first search. All the pairs have been visited
c     will be labeled in the order they are visited, and associated a
c     lowlink for each vertex, stored in the stack - lowlnk.
c
c     input variables :
c
c        nrows  -- number of rows in matrix
c        ncols  -- number of columns in matrix
c        nhcols -- number of columns in horizontal (underdetermined)
c                  partition
c        nhrows -- number of rows in horizontal (underdetermined)
c                  partition
c        nscols -- number of rows and columns in square partition
c        sqindx -- index for SR and SC, for rows and columns
c                  in the square partition
c        hrzcmp -- number of components in vertical partition
c        rowstr, colind
c               -- the adjacency structure, stored by rows
c        colset -- the row matched to a column (if any)
c
c     output variables :
c
c        sqcmpn -- number of components in the square partition
c        ccmstr -- global component start vector
c        rcmstr -- global component start vector
c        cnto   -- new to old mapping for columns
c        rnto   -- new to old mapping for rows
c
c     working variables  :
c
c        trycol -- pointer to next unsearched column for this row 
c        cbegin -- is the beginning of the component.
c        colmrk -- column mark vector.
c                  on input, is negative for all columns
c                            = sqindx for columns in sc
c                  used temporarily as a stack to
c                  store the depth-first numbering for each pair.
c                  that is, is the position of pair i in the stack
c                  if it is on it, is the permuted order of pair i for
c                  those pairs whose final position has been found and
c                  is otherwise zero for columns in sc and negative
c                  for all other columns.
c                  on output, is restored to original values
c        lowlnk -- stores the lowlink for each pair.
c                  is the smallest stack position of any pair to which
c                  a path from pair i has been found. it is set to n+1
c                  when pair i is removed from the stack.
c   	 prev   -- is the pair at the end of the path when pair i was 
c                  placed on the stack
c        
c     ==================================================================

c     --------------
c     ... parameters
c     --------------

      integer          nrows , ncols , nhcols, nhrows, nscols, 
     $                 sqindx, hrzcmp, sqcmpn

      integer          rowstr (nrows+1), colind (*), colset (ncols)

      integer          trycol (nrows), cbegin (nscols), lowlnk (ncols),
     $                 prev   (ncols) 

      integer          colmrk (ncols), cnto (ncols), rnto (nrows),
     $                 ccmstr (*)    , rcmstr (*)

c     -------------------
c     ... local variables
c     -------------------

      integer          cmpbeg, col   , compnt, count , passes, fcol  , 
     $                 fnlpos, frow  , rootcl, j     , pair  , scol  ,
     $                 stackp, xcol

c     ==================================================================

c     fnlpos  is the number of pairs whose positions in final ordering
c             have been found.
c     sqcmpn  is the number of components that have been found.
c     count   is the number of pairs on the stack  (stack pointer)

c     ------------------------------------------------------
c     ... initialization for columns in the square partition
c     ------------------------------------------------------

      fnlpos = 0
      sqcmpn = 0

      do 100 col = 1, ncols
         if  ( colmrk (col) .eq. sqindx )  then
            colmrk (col) = 0
         endif
 100  continue

      do 200 j = 1, nrows
         trycol (j) = rowstr (j)
 200  continue

c     ----------------------------
c     ... look for a starting pair
c     ----------------------------

      do 700 rootcl = 1, ncols

         if  ( colmrk (rootcl) .eq. 0 )  then

c           -------------------------------------
c           ... put pair (rootcl, colset(rootcl))
c               at beginning of stack
c           -------------------------------------

            fcol            = rootcl
            count           = 1
            lowlnk (fcol)   = count
            colmrk (fcol)   = count
            cbegin (nscols) = fcol

c           --------------------------------------------
c           ... the body of this loop puts a new pair on
c               the stack or backtracks
c           --------------------------------------------

            do 600 passes = 1, 2*nscols - 1

               frow = colset (fcol)
               
c              -------------------------------
c              ... have all edges leaving pair  
c                  (frow,fcol)  been searched?
c              -------------------------------

               if  ( trycol (frow) .gt. 0 )  then

c                 -----------------------------------------------
c                 ... look at edges leaving from row "frow" until
c                     we find a new column "scol" that has not
c                     yet been encountered or until all edges are
c                     exhausted.
c                 -----------------------------------------------

                  do 300 xcol = trycol (frow), rowstr (frow+1)-1
                     
                     scol = colind (xcol)
                     if  ( colmrk (scol) .eq. 0 )  then

c                       --------------------------------------
c                       ... put new pair  (scol, colset(scol))
c                           on the stack
c                       --------------------------------------

                        trycol (frow)           = xcol + 1
                        prev (scol)             = fcol
                        fcol                    = scol
                        count                   = count + 1 
                        lowlnk (fcol)           = count 
                        colmrk (fcol)           = count 
                        cbegin (nscols+1-count) = fcol 
                        go to 600

                     else
     $               if  ( colmrk (scol) .gt. 0 )  then

c                       -------------------------------------------
c                       ... has scol been on stack already?  then
c                           update value of low (fcol) if necessary
c                       -------------------------------------------

                        if  (lowlnk (scol) .lt. lowlnk (fcol))  then
                           lowlnk (fcol) = lowlnk (scol)
                        endif
                     endif
 300              continue

c                 ----------------------------------------
c                 ... there are no more edges leaving frow
c                 ----------------------------------------

                  trycol (frow) = -1

               endif

c              ------------------------------
c              
c              ------------------------------

               if  ( lowlnk (fcol) .ge. colmrk (fcol) )  then
c
c                 -----------------------------------------------------
c                 ... is  frow  the root of a block?  if so, we have 
c                     found a component.  order the nodes in this
c                     block by starting at the top of the stack and
c                     working down to the root of the block
c                 -----------------------------------------------------

                  sqcmpn = sqcmpn + 1
                  cmpbeg = fnlpos + 1
                  
                  do 400 stackp = nscols + 1 - count, nscols
                     pair          = cbegin (stackp)
                     fnlpos        = fnlpos + 1
                     colmrk (pair) = fnlpos
                     count         = count-1
                     lowlnk (pair) = nscols + 1
                     if  ( pair .eq. fcol )  go to 500
 400              continue

c                 -------------------------------------------------------
c                 ... record the starting position for the new component
c                 -------------------------------------------------------

 500              cbegin (sqcmpn) = cmpbeg 

c                 --------------------------------------------
c                 ... are there any pairs left on the stack.
c                     if so, backtrack.
c                     if not, have all the pairs been ordered?
c                 --------------------------------------------

                  if  ( count .eq. 0 )  then
                     if  ( fnlpos .lt. nscols )  then
                        go to 700
                     else
                        go to 800
                     endif
                  endif
                  
               endif
c
c              --------------------------------------
c              ... backtrack to previous pair on path
c              --------------------------------------

               scol = fcol
               fcol = prev (fcol)
               if  ( lowlnk (scol) .lt. lowlnk (fcol) )  then
                  lowlnk (fcol) = lowlnk (scol)
               endif

 600         continue
            
         endif

 700  continue

c     ----------------------------------------
c     ... put permutation in the required form
c     ----------------------------------------

 800  do 900 compnt = 1, sqcmpn
         ccmstr (compnt + hrzcmp) = (cbegin (compnt) + nhcols)
         rcmstr (compnt + hrzcmp) = (cbegin (compnt) + nhcols) -
     $                              (nhcols - nhrows)
 900  continue

      ccmstr (hrzcmp + sqcmpn + 1)  = nhcols + nscols + 1
      rcmstr (hrzcmp + sqcmpn + 1)  = nhrows + nscols + 1

c     ------------------------------------------------------
c     ... note that columns not in the square partition have
c         colmrk set negative.  diagonal entries in the
c         square block all correspond to matching pairs.
c     ------------------------------------------------------

      do 1000 col = 1, ncols
         j = colmrk (col)
         if  ( j .gt. 0 )  then
            cnto (nhcols + j) = col
            rnto (nhrows + j) = colset (col)
            colmrk (col)      = sqindx
         endif
 1000 continue

      return

      end

