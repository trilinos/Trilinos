C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      INTEGER FUNCTION PLTITL(REALN)

      IF (PLTFRC(REALN).EQ.0.) THEN
         PLTITL = INT(REALN)
         RETURN

      END IF

      AREAL = ABS(REALN)
      PLTITL = INT(AREAL)
      IF (REALN.LT.0.) THEN
         PLTITL = -PLTITL - 1
      END IF

      RETURN

      END
