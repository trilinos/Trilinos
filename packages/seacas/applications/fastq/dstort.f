C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE DSTORT (X1, X2, X3, X4, Y1, Y2, Y3, Y4, VALUE)
C***********************************************************************

C  SUBROUTINE DSTORT = CALCULATES A DISTORTION METRIC FOR AN ELEMENT
C                    USING THE IDEAS IN THE PAPER BY ODDY, 1988.

C***********************************************************************

C  SETUP THE JACOBIAN MATRIX

      XJ11 = (X1 * .125) + (X2 * .375) - (X3 * .375) - (X4 * .125)
      XJ12 = (Y1 * .125) + (Y2 * .375) - (Y3 * .375) - (Y4 * .125)
      XJ21 = - (X1 * .375) + (X2 * .375) + (X3 * .125) - (X4 * .125)
      XJ22 = - (Y1 * .375) + (Y2 * .375) + (Y3 * .125) - (Y4 * .125)

C  NORMALIZE THE JACOBIAN WITH RESPECT TO THE ELEMENT SIZE

      DETERM = (XJ11 * XJ22) - (XJ12 * XJ21)
      IF (DETERM .LE. 0.) THEN
         VALUE = 1.0E10
         RETURN
      ENDIF
      FACTOR = 1. / SQRT (DETERM)
      XJ11 = XJ11 * FACTOR
      XJ12 = XJ12 * FACTOR
      XJ21 = XJ21 * FACTOR
      XJ22 = XJ22 * FACTOR

C  NOW USE THE SECOND INVARIANT OF GREEN'S STRAIN

      C11 = XJ11*XJ11 + XJ21*XJ21
      C12 = XJ11*XJ12 + XJ21*XJ22
      C22 = XJ12*XJ12 + XJ22*XJ22

      VALUE = C11**2 + 2.*(C12**2) + C22**2 -
     &   (.5 * (C11+C22)**2 )
      VALUE = AMAX1 (VALUE, 0.)

      RETURN

      END
