C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE COMSRT (MXND, MXCORN, MXPICK, MLN, LNODES, LCORN,
     &   NCORN, ICOMB, ITYPE, NPICK)
C***********************************************************************

C  SUBROUTINE COMSRT = THIS SUBROUTINE GETS ALL THE COMBINATIONS
C                      POSSIBLE OF CORNERS, SIDES, AND DISSECION
C                      NODES

C***********************************************************************

C  VARIABLES USED:
C     ICOMB  = THE DIFFERENT VALID PRIMITIVE COMBINATIONS
C                 ICOMB (I, J) WHERE I IS THE COMBINATION NUMBER
C                              AND J IS THE CORNER NUMBER.
C                 THE VALUE OF ICOMB (I,J) IS 1 FOR CORNER AND 0 FOR
C                 FOR A SIDE INTERPRETATION
C     ITYPE  = THE TYPE OF PRIMITIVE OR NUMBER OF CORNERS IN THIS
C              COMBINATION.
C                  = 0 FOR A CIRCLE
C                  = 1 FOR A TEARDROP
C                  = 2 FOR A SEMICIRCLE
C                  = 3 FOR A TRIANGLE
C                  = 4 FOR A RECTANGLE
C                  = >4 OTHERWISE

C***********************************************************************

      DIMENSION LNODES (MLN, MXND), ICOMB (MXCORN, MXPICK)
      DIMENSION ITYPE (MXPICK), LCORN (MXCORN)

      NPICK = 1

      DO 100 I = 1, MXPICK
         ITYPE (I) = 0
  100 CONTINUE

      DO 150 I = 1, NCORN

C  PUT PURE CORNER AND CORNER/SIDE DESIGNATIONS IN ORDER

         ITEST = LNODES (6, LCORN (I))
         IF (ITEST .LE. 2) THEN
            DO 110 J = 1, NPICK
               ICOMB (I, J) = 1
               ITYPE (J) = ITYPE (J) + 1
  110       CONTINUE
            IF (ITEST .EQ. 2) THEN
               DO 130 J = NPICK + 1, NPICK * 2
                  DO 120 K = 1, I - 1
                     ICOMB (K, J) = ICOMB (K, J - NPICK)
  120             CONTINUE
                  ITYPE (J) = ITYPE (J - NPICK) - 1
                  ICOMB (I, J) = 0
  130          CONTINUE
               NPICK = NPICK * 2
            ENDIF

C  PUT PURE SIDE AND SIDE/DISSECTIONS DESIGNATIONS IN ORDER

         ELSEIF (ITEST .LE. 4) THEN
            DO 140 J = 1, NPICK
               ICOMB (I, J) = 0
  140       CONTINUE
         ENDIF

  150 CONTINUE

      RETURN

      END
