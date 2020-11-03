C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE PLTCG2(N,XV,YV,NO,XVO,YVO,C1,C2)
      INTEGER N
      REAL XV(*),YV(*)
      INTEGER NO
      REAL XVO(*),YVO(*)
      REAL C1(2),C2(2)
      REAL S(2),P(2)
      CHARACTER*6 SUBNAM
      PARAMETER (SUBNAM='PLTCG2')

      NMAX = NO
      NO = 0
      S(1) = XV(N)
      S(2) = YV(N)
      DO 2000 I = 1,N
         P(1) = XV(I)
         P(2) = YV(I)
         CP = (C2(1)-C1(1))* (P(2)-C1(2)) - (C2(2)-C1(2))* (P(1)-C1(1))
         IF (CP.GE.0.) THEN
            CP = (C2(1)-C1(1))* (S(2)-C1(2)) -
     *           (C2(2)-C1(2))* (S(1)-C1(1))
            IF (CP.GE.0.) THEN
               NO = NO + 1
               IF (NO.GT.NMAX) THEN
                  CALL PLTFLU
                  CALL SIORPT(SUBNAM,
     *   'Not enough space for output vertices; polygon clip unfinished'
     *                        ,2)
                  RETURN

               END IF

               XVO(NO) = P(1)
               YVO(NO) = P(2)

            ELSE
               FP = (S(2)-C1(2))* (C2(1)-C1(1)) -
     *              (S(1)-C1(1))* (C2(2)-C1(2))
               FQ = (P(2)-C1(2))* (C2(1)-C1(1)) -
     *              (P(1)-C1(1))* (C2(2)-C1(2))
               XL = FQ/ (FQ-FP)
               NO = NO + 1
               IF (NO.GT.NMAX) THEN
                  CALL PLTFLU
                  CALL SIORPT(SUBNAM,
     *   'Not enough space for output vertices; polygon clip unfinished'
     *                        ,2)
                  RETURN

               END IF

               XVO(NO) = P(1) + XL* (S(1)-P(1))
               YVO(NO) = P(2) + XL* (S(2)-P(2))
               NO = NO + 1
               IF (NO.GT.NMAX) THEN
                  CALL PLTFLU
                  CALL SIORPT(SUBNAM,
     *   'Not enough space for output vertices; polygon clip unfinished'
     *                        ,2)
                  RETURN

               END IF

               XVO(NO) = P(1)
               YVO(NO) = P(2)
            END IF

         ELSE
            CP = (C2(1)-C1(1))* (S(2)-C1(2)) -
     *           (C2(2)-C1(2))* (S(1)-C1(1))
            IF (CP.GE.0.) THEN
               FP = (S(2)-C1(2))* (C2(1)-C1(1)) -
     *              (S(1)-C1(1))* (C2(2)-C1(2))
               FQ = (P(2)-C1(2))* (C2(1)-C1(1)) -
     *              (P(1)-C1(1))* (C2(2)-C1(2))
               XL = FQ/ (FQ-FP)
               NO = NO + 1
               IF (NO.GT.NMAX) THEN
                  CALL PLTFLU
                  CALL SIORPT(SUBNAM,
     *   'Not enough space for output vertices; polygon clip unfinished'
     *                        ,2)
                  RETURN

               END IF

               XVO(NO) = P(1) + XL* (S(1)-P(1))
               YVO(NO) = P(2) + XL* (S(2)-P(2))
            END IF

         END IF

         S(1) = P(1)
         S(2) = P(2)
 2000 CONTINUE
      RETURN

      END
