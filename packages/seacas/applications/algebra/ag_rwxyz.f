C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details
C=======================================================================
      SUBROUTINE RWXYZ (NDBIN, NDBOUT, NDIM, NUMNP, NUMNPO,
     &                  IXNODE, CORD, CRDSCR)
C=======================================================================

C   --*** RWXYZ *** (ALGEBRA) Read and write database coordinates
C   --   Written by Amy Gilkey - revised 11/30/87
C   --   Modified for EXODUSIIV2 format 8/29/95
C   --
C   --RWXYZ reads and writes the coordinate array to the database.
C   --Deleted nodes are removed.
C   --
C   --Parameters:
C   --   NDBIN, NDBOUT - IN - the input and output database file
C   --   NDIM   - IN - the number of coordinates per node
C   --   NUMNP  - IN - the number of nodes
C   --   NUMNPO - IN - the number of nodes
C   --   IXNODE - IN - the indices of the output nodes (iff NUMNPO <> NUMNP)
C   --   CORD   - SCRATCH - coordinate I/O
C   --   CRDSCR - SCRATCH - coordinate I/O

      INTEGER NDBIN, NDBOUT
      INTEGER NDIM
      INTEGER NUMNP
      INTEGER IXNODE(*)
      REAL    CORD(NUMNP,NDIM)
      REAL    CRDSCR(NUMNPO,NDIM)

      if (ndim .eq. 2) then
         CALL EXGCOR(ndbin, cord(1,1), cord(1,2), rdum, ierr)
      else if (ndim .eq. 3) then
         CALL EXGCOR(ndbin, cord(1,1), cord(1,2),
     &               cord(1,3), ierr)
      else
         call prterr('FATAL', 'Illegal model dimension')
         RETURN
      end if

      IF ((NUMNPO .GT. 0) .AND. (NDIM .GT. 0)) THEN
         IF (NUMNP .EQ. NUMNPO) THEN
           if (ndim .eq. 2) then
              CALL EXPCOR(ndbout, cord(1,1), cord(1,2), rdum, ierr)
           else if (ndim .eq. 3) then
              CALL EXPCOR(ndbout, cord(1,1), cord(1,2),
     &                    cord(1,3), ioerr)
           else
              call prterr('FATAL', 'Illegal model dimension')
              RETURN
           end if
         ELSE
           do 20 idim=1, ndim
              do 10 ix=1, numnpo
                 crdscr(ix,idim) = cord(ixnode(ix),idim)
 10           continue
 20        continue
           if (ndim .eq. 2) then
              CALL EXPCOR(ndbout, crdscr(1,1), crdscr(1,2),
     &                    rdum, ioerr)
           else if (ndim .eq. 3) then
              CALL EXPCOR(ndbout, crdscr(1,1), crdscr(1,2),
     &                    crdscr(1,3), ioerr)
           else
              call prterr('FATAL', 'Illegal model dimension')
           end if
         END IF
      ELSE
         continue
      END IF

      RETURN
      END
