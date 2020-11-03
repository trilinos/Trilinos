C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE ELESTF (NLNKF, LINKF1, XN, YN, ZN)
C=======================================================================

C   --*** ELESTF *** (DETOUR) Plot element symbol for state for a face
C   --   Written by Amy Gilkey, revised 10/21/87
C   --   D. P. Flanagan, 12/08/83
C   --
C   --ELESTF plots the element symbol for the state for a face.
C   --
C   --Parameters:
C   --   NLNKF - IN - the number of nodes per face
C   --   LINKF1 - IN - the connectivity for the face
C   --   XN, YN, ZN - IN - the nodal coordinates

      INTEGER LINKF1(*)
      INTEGER NLNKF
      REAL XN(*), YN(*), ZN(*)

      N2 = LINKF1(NLNKF)
      DO 100 ILINK = 1, NLNKF
         N1 = LINKF1(ILINK)
         CALL MPD2VC (1, XN(N1), YN(N1), XN(N2), YN(N2))
         N2 = N1
  100 CONTINUE
      IF (NLNKF .EQ. 4) THEN
         N1 = LINKF1(1)
         N2 = LINKF1(3)
         CALL MPD2VC (1, XN(N1), YN(N1), XN(N2), YN(N2))
         N1 = LINKF1(2)
         N2 = LINKF1(4)
         CALL MPD2VC (1, XN(N1), YN(N1), XN(N2), YN(N2))
      END IF

      RETURN
      END
