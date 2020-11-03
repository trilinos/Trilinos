C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE MNMXK (NPELEM, NPNODE, NNXK, NXK, XN, YN, CENTK, KKK,
     &   XMIN, XMAX, YMIN, YMAX)
C***********************************************************************

C  SUBROUTINE MNMXK = FINDS MIN AND MAX DIMENSIONS FOR FLAGGED ELEMENTS

C**********************************************************************

      DIMENSION NXK (NNXK, NPELEM), CENTK (2, NPELEM)
      DIMENSION XN (NPNODE), YN (NPNODE)

C  FIND THE FIRST ELEMENT TO BE PLOTTED

      DO 150 I = 1, KKK
         IF (NXK (1, I) .GT. 0) THEN
            JX1 = I
            JX2 = I
            JY1 = I
            JY2 = I

C  COMPARE CENTERS TO GET MIN AND MAX ELEMENTS

            DO 100 J = I + 1, KKK
               IF (NXK (1, J) .GT. 0) THEN
                  IF (CENTK (1, J) .LT. CENTK (1, JX1))JX1 = J
                  IF (CENTK (1, J) .GT. CENTK (1, JX2))JX2 = J
                  IF (CENTK (2, J) .LT. CENTK (2, JY1))JY1 = J
                  IF (CENTK (2, J) .GT. CENTK (2, JY2))JY2 = J
               ENDIF
  100       CONTINUE

C  FIND CORRECT MIN AND MAX FROM NODES OF MIN AND MAX ELEMENTS

            XMIN = XN (NXK (1, JX1))
            DO 110 K = 2, NNXK
               IF (NXK (K, JX1) .GT. 0)
     &            XMIN = AMIN1 (XMIN, XN (NXK (K, JX1)))
  110       CONTINUE
            XMAX = XN (NXK (1, JX2))
            DO 120 K = 2, NNXK
               IF (NXK (K, JX2) .GT. 0)
     &            XMAX = AMAX1 (XMAX, XN (NXK (K, JX2)))
  120       CONTINUE
            YMIN = YN (NXK (1, JY1))
            DO 130 K = 2, NNXK
               IF (NXK (K, JY1) .GT. 0)
     &            YMIN = AMIN1 (YMIN, YN (NXK (K, JY1)))
  130       CONTINUE
            YMAX = YN (NXK (1, JY2))
            DO 140 K = 2, NNXK
               IF (NXK (K, JY2) .GT. 0)
     &            YMAX = AMAX1 (YMAX, YN (NXK (K, JY2)))
  140       CONTINUE
            RETURN
         ENDIF
  150 CONTINUE

C  RETURN WITH DEFAULT MINS AND MAXS

      XMIN = 0.
      XMAX = 1.
      YMIN = 0.
      YMAX = 1.
      RETURN

      END
