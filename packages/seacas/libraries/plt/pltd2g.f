C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      LOGICAL FUNCTION PLTD2G(XD,YD,XG,YG)
      DIMENSION UMAP(14), TYPE(1)

      CALL PLTGTG(27,UMAP)
      CALL PLTGTG(9,TYPE)
      PLTD2G = .FALSE.
      TOP = XD*UMAP(4) - YD*UMAP(3) - UMAP(5)*UMAP(4) + UMAP(6)*UMAP(3)
      BOTTOM = UMAP(1)*UMAP(4) - UMAP(2)*UMAP(3)
      XG = TOP/BOTTOM
      TOP = XD*UMAP(2) - YD*UMAP(1) - UMAP(5)*UMAP(2) + UMAP(6)*UMAP(1)
      BOTTOM = UMAP(3)*UMAP(2) - UMAP(4)*UMAP(1)
      YG = TOP/BOTTOM
      IF (TYPE(1).EQ.2.) THEN
         IF (XD.LE.0.) THEN
            CALL PLTFLU
            CALL SIORPT('PLTD2G',
     *                 'Cannot convert a negative number on log x axis.'
     *                  ,2)
            RETURN

         END IF

         XG = 10.**XG

      ELSE IF (TYPE(1).EQ.3.) THEN
         IF (YD.LE.0.) THEN
            CALL PLTFLU
            CALL SIORPT('PLTD2G',
     *                 'Cannot convert a negative number on log y axis.'
     *                  ,2)
            RETURN

         END IF

         YG = 10.**YG

      ELSE IF (TYPE(1).EQ.4.) THEN
         IF (XD.LE.0.) THEN
            CALL PLTFLU
            CALL SIORPT('PLTD2G',
     *                 'Cannot convert a negative number on log x axis.'
     *                  ,2)
            RETURN

         END IF

         IF (YD.LE.0.) THEN
            CALL PLTFLU
            CALL SIORPT('PLTD2G',
     *                 'Cannot convert a negative number on log y axis.'
     *                  ,2)
            RETURN

         END IF

         XG = 10.**XG
         YG = 10.**YG
      END IF

      PLTD2G = .TRUE.
      RETURN

      END
