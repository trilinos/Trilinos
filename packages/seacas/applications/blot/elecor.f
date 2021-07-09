C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE ELECOR (NDIM, NELBLK, LEN, NLNK, LINK,
     &   XN, YN, ZN, XE, YE, ZE)
C=======================================================================

C   --*** ELECOR *** (BLOT) Calculate element centers
C   --   Written by Amy Gilkey, revised 10/23/87
C   --
C   --ELECOR calculates the element centers from the nodal coordinates.
C   --
C   --Parameters:
C   --   NDIM - IN - the number of dimensions
C   --   NELBLK - IN - the number of element blocks
C   --   LEN - IN - the cumulative element counts by element block
C   --   NLNK - IN - the number of nodes per element
C   --   LINK - IN - the connectivity array; connectivity all zero if
C   --      element is undefined
C   --   XN, YN, ZN - IN - the nodal coordinates
C   --   XE, YE, ZE - OUT - the element center coordinates

      INTEGER LEN(0:*), LINK(*)
      INTEGER NLNK(*)
      REAL XN(*), YN(*), ZN(*)
      REAL XE(*), YE(*), ZE(*)

      DO 170 IELB = 1, NELBLK
         IF (NLNK(IELB) .LE. 0) GOTO 170

         DIVLNK = 1.0 / NLNK(IELB)

         IXL0 = IDBLNK (IELB, 0, LEN, NLNK) - 1
         IF (NDIM .EQ. 2) THEN
            IF (NLNK(IELB) .NE. 8) THEN
               DO 110 IEL = LEN(IELB-1)+1, LEN(IELB)
                  X0 = 0.0
                  Y0 = 0.0
                  IF (LINK(IXL0+1) .NE. 0) THEN
                     DO 100 I = 1, NLNK(IELB)
                        X0 = X0 + XN(LINK(IXL0+I))
                        Y0 = Y0 + YN(LINK(IXL0+I))
  100                CONTINUE
                  END IF
                  XE(IEL) = X0 * DIVLNK
                  YE(IEL) = Y0 * DIVLNK
                  IXL0 = IXL0 + NLNK(IELB)
  110          CONTINUE

            ELSE IF (NLNK(IELB) .EQ. 8) THEN
C            --Note that 8-node elements are numbered with corners at 1-3-5-7
               DO 140 IEL = LEN(IELB-1)+1, LEN(IELB)
                  X0 = 0.0
                  Y0 = 0.0
                  IF (LINK(IXL0+1) .NE. 0) THEN
                     DO 120 I = 1, NLNK(IELB), 2
                        X0 = X0 - 0.25 * XN(LINK(IXL0+I))
                        Y0 = Y0 - 0.25 * YN(LINK(IXL0+I))
  120                CONTINUE
                     DO 130 I = 2, NLNK(IELB), 2
                        X0 = X0 + 0.5 * XN(LINK(IXL0+I))
                        Y0 = Y0 + 0.5 * YN(LINK(IXL0+I))
  130                CONTINUE
                  END IF
                  XE(IEL) = X0
                  YE(IEL) = Y0
                  IXL0 = IXL0 + NLNK(IELB)
  140          CONTINUE
            END IF

         ELSE IF (NDIM .EQ. 3) THEN
            DO 160 IEL = LEN(IELB-1)+1, LEN(IELB)
               X0 = 0.0
               Y0 = 0.0
               Z0 = 0.0
               IF (LINK(IXL0+1) .NE. 0) THEN
                  DO 150 I = 1, NLNK(IELB)
                     X0 = X0 + XN(LINK(IXL0+I))
                     Y0 = Y0 + YN(LINK(IXL0+I))
                     Z0 = Z0 + ZN(LINK(IXL0+I))
  150             CONTINUE
               END IF
               XE(IEL) = X0 * DIVLNK
               YE(IEL) = Y0 * DIVLNK
               ZE(IEL) = Z0 * DIVLNK
               IXL0 = IXL0 + NLNK(IELB)
  160       CONTINUE
         END IF
  170 CONTINUE

      RETURN
      END
