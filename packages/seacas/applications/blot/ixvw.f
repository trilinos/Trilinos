C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      INTEGER FUNCTION IXVW (INCEMP, IVWNUM)
C=======================================================================

C   --*** IXVW *** (MESH) Return the index of the Nth defined view
C   --   Written by Amy Gilkey - revised 10/08/87
C   --
C   --IXVW returns the index of IVWNUMth defined view, which may include
C   --empty views.  If there are not IVWNUM defined views, zero is returned.
C   --
C   --Parameters:
C   --   INCEMP - IN - a view with an empty mode are returned iff true
C   --   IVWNUM - IN - the view number
C   --
C   --Common Variables:
C   --   Uses MSHDEF of /MSHOPT/

      include 'mshopt.blk'

      LOGICAL INCEMP

      IVW = 0
      DO 100 IXVW = 1, 4
         IF (MSHDEF(IXVW) .NE. 'NONE') THEN
            IF (INCEMP .OR. (MSHDEF(IXVW) .NE. 'EMPTY')) THEN
               IVW = IVW + 1
               IF (IVW .GE. IVWNUM) GOTO 110
            END IF
         END IF
  100 CONTINUE
      IXVW = 0
  110 CONTINUE

      RETURN
      END
