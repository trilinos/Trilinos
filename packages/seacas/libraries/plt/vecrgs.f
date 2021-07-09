C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE VECRGS(ND,V,VMAX,VMIN)
      DIMENSION V(*)

      IF (ND.GE.1) THEN

         VMAX = V(1)
         VMIN = V(1)
         I = 1
   10    CONTINUE
         IF ((I.LE.ND)) THEN
            T = V(I)
            VMAX = MAX(VMAX,T)
            VMIN = MIN(VMIN,T)
            I = I + 1
            GO TO 10

         END IF

      END IF

      END
