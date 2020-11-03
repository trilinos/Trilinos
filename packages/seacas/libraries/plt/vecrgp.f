C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE VECRGP(ND,V,VMAX,VMIN)
      DIMENSION V(*)

      IF (ND.LT.1) THEN
         RETURN

      END IF

      VMAX = V(1)
      DO 3070 I = 1,ND
         VMAX = MAX(VMAX,V(I))
 3070 CONTINUE
      IF (VMAX.LT.0) THEN
         VMAX = 1.
         VMIN = .1
         RETURN

      END IF

      VMIN = VMAX
      DO 3090 I = 1,ND
         IF (V(I).GT.0.) THEN
            VMIN = MIN(VMIN,V(I))
         END IF

 3090 CONTINUE
      IF (VMIN.EQ.VMAX) THEN
         VMIN = .1*VMAX
      END IF

      RETURN

      END
