C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      INTEGER FUNCTION NDEFVW (INCEMP)
C=======================================================================

C   --*** NDEFVW *** (MESH) Return the number of defined views
C   --   Written by Amy Gilkey - revised 10/08/87
C   --
C   --NDEFVW returns the number of defined views, which may include
C   --empty views.
C   --
C   --Parameters:
C   --   INCEMP - IN - empty modes are included in the number of modes
C   --      iff true
C   --
C   --Common Variables:
C   --   Uses MSHDEF of /MSHOPT/

      include 'mshopt.blk'

      LOGICAL INCEMP

      INTEGER NUMMOD

      NDEFVW = 4 - NUMMOD (MSHDEF, ' ', 'NONE', ' ')
      IF (.NOT. INCEMP)
     &   NDEFVW = NDEFVW - NUMMOD (MSHDEF, ' ', 'EMPTY', ' ')

      RETURN
      END
