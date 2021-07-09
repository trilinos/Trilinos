C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE PLTLI1(PLL,PUR,N,XV,YV,NO,XVO,YVO)
      CHARACTER*6 SUBNAM
      PARAMETER (SUBNAM='PLTLI1')
      DIMENSION PLL(2),PUR(2),XV(*),YV(*),XVO(*),YVO(*)
      DIMENSION P(2),S(2)
      LOGICAL INSIDE

      NMAX = NO
      NO = 0
      S(1) = XV(N)
      S(2) = YV(N)
      DO 2220 I = 1,N
         P(1) = XV(I)
         P(2) = YV(I)
         INSIDE = P(2) .GE. PLL(2)
         IF (INSIDE) THEN
            INSIDE = S(2) .GE. PLL(2)
            IF (INSIDE) THEN
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
               TEMP = PUR(1) - PLL(1)
               FP = (S(2)-PLL(2))*TEMP
               FQ = (P(2)-PLL(2))*TEMP
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
            INSIDE = S(2) .GE. PLL(2)
            IF (INSIDE) THEN
               TEMP = PUR(1) - PLL(1)
               FP = (S(2)-PLL(2))*TEMP
               FQ = (P(2)-PLL(2))*TEMP
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
 2220 CONTINUE
      RETURN

      END
