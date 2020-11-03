C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE NORMND (NPNODE, BMESUR, RMAX)
C***********************************************************************

C  SUBROUTINE NORMND = NORMALIZES A NODE VARIABLE

C***********************************************************************

      DIMENSION BMESUR(NPNODE)

      BMIN = BMESUR(1)
      BMAX = BMESUR(1)
      DO 100 NODE = 2, NPNODE
         BMAX = AMAX1 (BMESUR(NODE), BMAX)
         BMIN = AMIN1 (BMESUR(NODE), BMIN)
  100 CONTINUE

      BMAX = BMAX - BMIN
      DO 110 NODE = 1, NPNODE
         BMESUR(NODE) = BMESUR(NODE) - BMIN
  110 CONTINUE

C  RMAX = MAXIMUM RATIO FOR PLATEAU VALUES

      DO 120 NODE = 1, NPNODE
         IF (BMESUR (NODE) .GE. (BMAX * RMAX)) THEN
            BMESUR(NODE) = 1.0
         ELSE
            BMESUR (NODE) = BMESUR(NODE) / (BMAX * RMAX)
         ENDIF
  120 CONTINUE

      RETURN

      END
