C $Id: dstort.f,v 1.2 2000/11/13 15:39:04 gdsjaar Exp $
C $Log: dstort.f,v $
C Revision 1.2  2000/11/13 15:39:04  gdsjaar
C Cleaned up unused variables and labels.
C
C Removed some real to int conversion warnings.
C
C Revision 1.1.1.1  1990/11/30 11:06:27  gdsjaar
C FASTQ Version 2.0X
C
c Revision 1.1  90/11/30  11:06:25  gdsjaar
c Initial revision
c 
C
CC* FILE: [.PAVING]DSTORT.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE DSTORT (X1, X2, X3, X4, Y1, Y2, Y3, Y4, VALUE)
C***********************************************************************
C
C  SUBROUTINE DSTORT = CALCULATES A DISTORTION METRIC FOR AN ELEMENT
C                    USING THE IDEAS IN THE PAPER BY ODDY, 1988.
C
C***********************************************************************
C
C  SETUP THE JACOBIAN MATRIX
C
      XJ11 = (X1 * .125) + (X2 * .375) - (X3 * .375) - (X4 * .125)
      XJ12 = (Y1 * .125) + (Y2 * .375) - (Y3 * .375) - (Y4 * .125)
      XJ21 = - (X1 * .375) + (X2 * .375) + (X3 * .125) - (X4 * .125)
      XJ22 = - (Y1 * .375) + (Y2 * .375) + (Y3 * .125) - (Y4 * .125)
C
C  NORMALIZE THE JACOBIAN WITH RESPECT TO THE ELEMENT SIZE
C
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
C
C  NOW USE THE SECOND INVARIANT OF GREEN'S STRAIN
C
      C11 = XJ11*XJ11 + XJ21*XJ21
      C12 = XJ11*XJ12 + XJ21*XJ22
      C22 = XJ12*XJ12 + XJ22*XJ22
C
      VALUE = C11**2 + 2.*(C12**2) + C22**2 -
     &   (.5 * (C11+C22)**2 )
      VALUE = AMAX1 (VALUE, 0.)
C
      RETURN
C
      END
