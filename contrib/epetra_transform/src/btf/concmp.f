      subroutine   concmp   ( cmbase, rnbase, cnbase, vindex, nrows ,
     $                        ncols , nvrows, nvcols, rowstr, colidx,
     $                        colstr, rowidx, predrw, nextrw, predcl,    ,
     $                        nextcl, ctab  , rtab  , colmrk, rowmrk,
     $                        cmclad, cmrwad, cnto  , rnto  , numcmp ) 

c     ==================================================================
c     ==================================================================
c     ====  concmp -- find the connected components in the          ====
c     ====            vertical (horizontal) block                   ====
c     ==================================================================
c     ==================================================================

c     original -- alex pothen and chin-ju fan, penn state, 1988
c     bcs modifications, john lewis, sept. 19, 1990

c     concmp:  find the connected components in the subgraph spanned
c              by the rows and columns in the vertical block.  the
c              same subroutine is used to find the connected
c              components in the horizontal block -- the transpose
c              of the matrix is used for that case.
c
c     input variables:
c
c         cmbase -- the number of components found in previous fine
c                   analysis of the coarse partition
c         rnbase -- the number of rows in earlier numbered partitions
c                   (0 for the horizontal block, nhrows+nsrows for
c                    the vertical partition)
c         cnbase -- the number of columns in earlier numbered partitions
c         vindex -- used to check whether the nodes belong in the
c                   vertical block
c         nrows  -- number of rows in the matrix 
c         ncols  -- number of columns in the matrix 
c         nvrows -- number of rows in the vertical block
c         nvcols -- number of columns in the vertical block
c         rowstr, colidx
c               -- the adjacency structure of the matrix using
c                  row-wise storage
c         colstr, rowidx
c               -- the adjacency structure of the matrix using
c                  column-wise storage
c
c     output variables:
c
c        numcmp  -- number of connected components
c        colmrk  -- initially,                        
c                    colmrk(i) = vindex if i belongs to vc.
c                              < 0 otherwise.
c                    during execution, 
c                    colmrk(i) = j, if i belongs to the jth component.
c                    after execution, original values restored
c        rowmrk -- initially,                        
c                    rowmrk(i) = vindex if i belongs to vr.
c                              < 0  otherwise.
c                    during execution, 
c                    rowmrk(i) = j, if i belongs to the jth component.
c                              < 0 otherwise.
c                    after execution, original values restored
c        cmclad, cmrwad 
c               -- the address (in the new ordering) of the 
c                  first column/row in each component,
c        cnto   -- the new to old mapping for the columns
c        rnto   -- the new to old mapping for the rows
c
c     working variables:
c
c        predrw, predcl
c               -- the path stack --
c                     predrw(i) = j means that we have in the path an
c                                   edge leaving from row node j to
c                                   column node i.
c                     predcl(i) = j means that we have in the path an 
c                                   edge leaving from column node j to
c                                   row node i.
c        nextcl -- nextcl(i) is index of first unsearched edge leaving
c                      from column node i.
c        nextrw -- nextrw(i) is index of first unsearched edge leaving
c                      from row node i.
c
c        ctab, rtab
c               -- temporary copy of the address (in the new ordering)
c                  of the first column/row in each component
c
c     ==================================================================

c     --------------
c     ... parameters
c     --------------

      integer         cmbase, rnbase, cnbase, vindex, nrows , ncols ,
     $                nvrows, nvcols, numcmp

      integer         colstr (nrows+1), rowstr (ncols+1), rowidx (*),
     $                colidx (*)

      integer         predrw (ncols), nextrw (nrows),
     $                predcl (nrows), nextcl (ncols),
     $                cmclad (ncols), cmrwad (nrows),
     $                colmrk (ncols), rowmrk (nrows),
     $                ctab   (*)    , rtab (*),
     $                cnto  (ncols) , rnto (nrows)

c     -------------------
c     ... local variables
c     -------------------

      integer         col, compn, p, cn, rn, row, xcol, xrow

c     ==================================================================

c     initialization
c     cn -- the number of the scanned column node
c     rn -- the number of the scanned row node

      cn     = 0
      rn     = 0
      numcmp = 0

c     ----------------------------------------------------------------
c     ... number of vertical rows > number of vertical columns.
c         start each search for a connected component with an unmarked
c         row in the vertical block.
c     ----------------------------------------------------------------


      do 500 p = 1, nrows

         if  ( rowmrk (p) .eq. vindex )  then

            row = p

c           --------------------------------------------------------
c           ... update the value of the current working component
c               put 'row' into the new component as the root of path
c           --------------------------------------------------------

            numcmp                   = numcmp + 1
            ctab (numcmp)            = cnbase + 1 + cn
            rtab (numcmp)            = rnbase + 1 + rn
            cmclad (cmbase + numcmp) = ctab (numcmp)
            cmrwad (cmbase + numcmp) = rtab (numcmp)
            rowmrk (row)             = numcmp
            rn                       = rn + 1
            nextrw (row)             = rowstr (row)
            predcl (row)             = 0

c           ------------------------------------------
c           ... from row node to col node --
c               try to find a forward step if possible
c               else backtrack
c           ------------------------------------------

 100        do 200 xcol = nextrw (row), rowstr (row + 1) -1
               col  = colidx (xcol)

               if  ( colmrk (col) .eq. vindex )  then

c                 ------------------------------------------------
c                 ... forward one step :
c                     find a forward step from row 'row' to column 
c                     'col'.  put 'col' into the current component
c                 ------------------------------------------------

                  nextrw (row) = xcol + 1
                  colmrk (col) = numcmp
                  cn           = cn + 1
                  nextcl (col) = colstr (col)
                  predrw (col) = row
                  go to 300

               endif
 200        continue
            
c           -----------------------------------------
c           ... backward one step  (back to col node)
c           -----------------------------------------

            nextrw (row) = rowstr (row + 1)
            col          = predcl (row)
            if  ( col .eq. 0 )  go to 500
                        
c           ------------------------------------------
c           ... from col node to row node
c               try to find a forward step if possible
c               else backtrack
c           ------------------------------------------

 300        do 400 xrow = nextcl (col), colstr (col + 1) - 1
               row = rowidx (xrow)
               if  ( rowmrk (row) .eq. vindex )  then

c                 --------------------------------------------------
c                 ... forward one step :
c                     find a forward step from column 'col' to
c                     row 'row'.  put row into the current component
c                 --------------------------------------------------

                  nextcl (col) = xrow + 1
                  rowmrk (row) = numcmp
                  rn           = rn + 1
                  nextrw (row) = rowstr (row)
                  predcl (row) = col
                  go to 100
               endif
 400        continue

c           -----------------------------------------
c           ... backward one step  (back to row node)
c           -----------------------------------------

            nextcl (col) = colstr (col + 1)
            row          = predrw (col)
            go to 100
            
         endif

  500 continue

c     ------------------------------------------------------------
c     ... generate the column and row permutations (cnto and rnto)
c         so that each component is numbered consecutively
c     ------------------------------------------------------------

      cmclad (cmbase + 1 + numcmp)  = cnbase + 1 + nvcols
      cmrwad (cmbase + 1 + numcmp)  = rnbase + 1 + nvrows

      do 600 col = 1, ncols
	 compn = colmrk (col)
	 if  ( compn .gt. 0 )  then
            cnto (ctab (compn)) = col
            ctab (compn)        = ctab (compn) + 1
            colmrk (col)        = vindex
         endif
  600 continue

      do 700 row = 1, nrows
	 compn = rowmrk (row)
	 if  ( compn .gt. 0 ) then
            rnto (rtab (compn)) = row
            rtab (compn)        = rtab (compn) + 1
            rowmrk (row)        = vindex
         endif
  700 continue

      return
      end

