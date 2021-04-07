C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE MATEDG (LENF, IELBST, IEDSET, NEDGES)
C=======================================================================

C   --*** MATEDG *** (MESH) Delete edges of selected element block
C   --   Written by Amy Gilkey - revised 07/06/87
C   --
C   --MATEDG deletes all edges that are of a selected element block.
C   --
C   --Parameters:
C   --   LENF - IN - the cumulative face counts by element block
C   --   IELBST - IN - the element block status (>0 if selected)
C   --   IEDSET - IN/OUT - the edge line set;
C   --      (0) = face defining edge; 0 to delete edge
C   --   NEDGES - IN/OUT - the number of lines in the edge set
C   --
C   --Common Variables:
C   --   Uses NELBLK of /DBNUMS/

      include 'debug.blk'
      include 'dbnums.blk'

      INTEGER LENF(0:NELBLK)
      INTEGER IELBST(NELBLK)
      INTEGER IEDSET(0:2,*)

      nhid = 0
      DO 120 IEDG = 1, NEDGES

         IFAC = IEDSET(0,IEDG)
         IF (IFAC .EQ. 0) GOTO 120

C      --Find the face element block
         DO 100 IELB = 1, NELBLK
            IF (IFAC .LE. LENF(IELB)) GOTO 110
  100    CONTINUE
  110    CONTINUE

C      --Delete edge if face is of a selected element block

         IF (IELBST(IELB) .GT. 0) THEN
            IEDSET(0,IEDG) = 0
            nhid = nhid + 1
         END IF
  120 CONTINUE
      if ((cdebug .eq. 'HIDDEN') .and. (idebug .ge. 1))
     &   write (*, '(1x,a,i5)') 'edges in selected block =', nhid

      RETURN
      END
