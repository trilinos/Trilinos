C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      LOGICAL FUNCTION PLTG2D(XG,YG,XD,YD)
      DIMENSION UMAP(14)
      REAL TYPE(1)

      PLTG2D = .FALSE.
      CALL PLTGTG(27,UMAP)
      CALL PLTGTG(9,TYPE)
      IF (TYPE(1).EQ.1.) THEN
         XD = XG*UMAP(1) + YG*UMAP(3) + UMAP(5)
         YD = XG*UMAP(2) + YG*UMAP(4) + UMAP(6)

      ELSE IF (TYPE(1).EQ.2.) THEN
         IF (XG.LE.0.) THEN
            CALL PLTFLU
            CALL SIORPT('PLTG2D',
     *                 'Cannot convert a negative number on log x axis.'
     *                  ,2)
            RETURN

         END IF

         XD = LOG10(XG)*UMAP(1) + YG*UMAP(3) + UMAP(5)
         YD = LOG10(XG)*UMAP(2) + YG*UMAP(4) + UMAP(6)

      ELSE IF (TYPE(1).EQ.3.) THEN
         IF (YG.LE.0.) THEN
            CALL PLTFLU
            CALL SIORPT('PLTG2D',
     *                 'Cannot convert a negative number on log y axis.'
     *                  ,2)
            RETURN

         END IF

         XD = XG*UMAP(1) + LOG10(YG)*UMAP(3) + UMAP(5)
         YD = XG*UMAP(2) + LOG10(YG)*UMAP(4) + UMAP(6)

      ELSE IF (TYPE(1).EQ.4.) THEN
         IF (XG.LE.0.) THEN
            CALL PLTFLU
            CALL SIORPT('PLTG2D',
     *                 'Cannot convert a negative number on log x axis.'
     *                  ,2)
            RETURN

         END IF

         IF (YG.LE.0.) THEN
            CALL PLTFLU
            CALL SIORPT('PLTG2D',
     *                 'Cannot convert a negative number on log y axis.'
     *                  ,2)
            RETURN

         END IF

         XD = LOG10(XG)*UMAP(1) + LOG10(YG)*UMAP(3) + UMAP(5)
         YD = LOG10(XG)*UMAP(2) + LOG10(YG)*UMAP(4) + UMAP(6)
      END IF

      PLTG2D = .TRUE.
      RETURN

      END
