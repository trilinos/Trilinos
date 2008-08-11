      subroutine   maxmatch   ( nrows , ncols , colstr, rowind, prevcl,
     $                          prevrw, marker, tryrow, nxtchp, rowset,
     $                          colset )
c
c     ==================================================================
c     ==================================================================
c     ====  maxmatch -- find maximum matching                       ====
c     ==================================================================
c     ==================================================================

c     maxmatch uses depth-first search to find an augmenting path from
c     each column node to get the maximum matching.
c
c     Alex Pothen and Chin-Ju Fan, Penn State University, 1988
c     last modifed: Alex Pothen July 1990 
c     last bcs modifications:  John Lewis, Sept. 1990
c 
c     input variables :
c
c        nrows -- number of row nodes in the graph.
c        ncols -- number of column nodes in the graph.
c        colstr, rowind -- adjacency structure of graph, stored by
c                          columns
c
c     output variables :      
c
c        rowset -- describe the matching.
c                  rowset (row) = col > 0 means column "col" is matched
c                                         to row "row"
c                               = 0       means "row" is an unmatched
c                                         node.
c        colset -- describe the matching.
c                  colset (col) = row > 0 means row "row" is matched to
c                                 column "col"
c                               = 0       means "col" is an unmatched
c                                         node.
c     Working variables :
c
c         prevrw (ncols) -- pointer toward the root of the depth-first 
c                           search from a column to a row.
c         prevcl (ncols) -- pointer toward the root of the depth-first 
c                           search from a column to a column.
c                           the pair (prevrw,prevcl) represent a
c                           matched pair.
c         marker (nrows) -- marker (row) <= the index of the root of the
c                           current depth-first search.  row has been
c                           visited in current pass when equality holds.
c         tryrow (ncols) -- tryrow (col) is a pointer into rowind to 
c                           the next row to be explored from column col
c                           in the depth-first search.
c         nxtchp (ncols) -- nxtchp (col) is a pointer into rowind to the
c                           next row to be explored from column col for
c                           the cheap assignment.  set to -1 when
c                           all rows have been considered for
c                           cheap assignment
c
c     ==================================================================

c     --------------
c     ... parameters
c     --------------

      integer        nrows, ncols

      integer        colstr (ncols+1), rowind (*), rowset (nrows),
     $               colset (ncols) 

      integer        prevrw (ncols), prevcl (ncols), tryrow (ncols),
     $               marker (nrows), nxtchp (ncols)

c     -------------------
c     ... local variables 
c     -------------------
c
      integer        nodec, col, nextrw, lastrw, xrow, row, nxtcol,
     $               prow, pcol
c
c     ==================================================================

      do 600 nodec = 1, ncols

c        --------------------------------------------------
c        ... initialize node 'col' as the root of the path.
c        --------------------------------------------------

         col          = nodec
         prevrw (col) = 0
         prevcl (col) = 0
         nxtchp (col) = colstr (col)

c        -----------------------------------------------------------
c        ... main loop begins here. Each time through, try to find a
c            cheap assignment from node col.
c        -----------------------------------------------------------

 100     nextrw = nxtchp (col)
         lastrw = colstr (col+1) - 1

         if  (nextrw .gt. 0 )  then

            do 200  xrow = nextrw, lastrw
               row = rowind (xrow)
               if  ( rowset (row) .eq. 0 )  go to 400
 200         continue

c           ------------------------------------------------
c           ... mark column when all adjacent rows have been
c               considered for cheap assignment.
c           ------------------------------------------------

            nxtchp (col)  = -1

         endif

c        ------------------------------------------------------------
c        ... Each time through, take a step forward if possible, or
c            backtrack if not .  Quit when backtracking takes us back 
c            to the beginning of the search.
c        ------------------------------------------------------------

         tryrow (col) = colstr (col)
         nextrw       = tryrow (col)
c$$$         lastrw = colstr (col+1) - 1

         if  ( lastrw .ge. nextrw )  then
            do 300 xrow = nextrw, lastrw
c              next line inserted by Alex Pothen, July 1990
c$$$               ii  = xrow
               row = rowind (xrow)
               if  ( marker (row) .lt. nodec )  then 
                  
c                 ---------------------------------------
c                 ... row is unvisited yet for this pass.
c                     take a forward step
c                 ---------------------------------------

                  tryrow (col) = xrow + 1
                  marker (row) = nodec  
                  nxtcol       = rowset (row)

                  if  ( nxtcol .lt. 0 )  then
                     go to 801
                  else
     $            if  ( nxtcol .eq. col )  then
                     go to 802
                  else
     $            if  ( nxtcol .gt. 0 )  then

c                    -----------------------------------------
c                    ... the forward step led to a matched row
c                        try to extend augmenting path from
c                        the column matched by this row.
c                    -----------------------------------------

                     prevcl (nxtcol) = col
                     prevrw (nxtcol) = row
                     tryrow (nxtcol) = colstr (nxtcol)
                     col             = nxtcol
                     go to 100

                  else

c                    -----------------
c                    ... unmatched row
c                    -----------------

                     go to 400

                  endif

               endif
 300        continue
         endif

c        ---------------------------------------------------
c        ... no forward step -- backtrack.
c            if we backtrack all the way, the search is done
c        ---------------------------------------------------
c
         col = prevcl (col)
         if  ( col .gt. 0 )  then
            go to 100
         else
            go to 600
         endif
 
c        ---------------------------------------------------
c        ... update the matching by alternating the matching
c            edge backward toward the root
c        ---------------------------------------------------

 400     rowset (row) = col
         prow         = prevrw (col)
         pcol         = prevcl (col)

 500         if  ( pcol .gt. 0 )  then
                if  ( rowset (prow) .ne. col ) go to 803
                rowset (prow) = pcol
                col           = pcol
                prow          = prevrw (pcol)
                pcol          = prevcl (pcol)
                go to 500
             endif

 600  continue

c     ------------------------------------------------------
c     ... compute the matching from the view of column nodes
c     ------------------------------------------------------

      do 700 row = 1, nrows
	col = rowset (row)
	if  ( col .gt. 0 )  then
           colset (col) = row
        endif
  700 continue

      return

c     -------------
c     ... bug traps
c     -------------

  801 write (6, 901)
  901 format (' bug in maxmatch : search reached a forbidden column')
      stop

  802 write (6, 902)
  902 format (' bug in maxmatch : search followed a matching edge')
      stop

  803 write (6, 903) col, row, row, rowset (row)
  903 format (' bug in maxmatch : pointer toward root disagrees with ',
     $        'matching.' /
     $        'prevcl (', i4, ')  = ', i4, ' but colset (', i4, ')  = ',
     $        i4)
      stop

      end

