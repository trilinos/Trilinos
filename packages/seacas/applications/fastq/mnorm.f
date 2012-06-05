C $Id: mnorm.f,v 1.1 1990/11/30 11:12:25 gdsjaar Exp $
C $Log: mnorm.f,v $
C Revision 1.1  1990/11/30 11:12:25  gdsjaar
C Initial revision
C
C
CC* FILE: [.QMESH]MNORM.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE MNORM (MXND, XN, YN, NXL, LLL, STDLEN)
C***********************************************************************
C
C  SUBROUTINE MNORM = FINDS THE AVERAGE LENGTH OF THOSE LINES NOT MUCH
C                     LONGER THAN THE AVERAGE
C
C***********************************************************************
C
      DIMENSION NXL (2, 3 * MXND), XN (MXND), YN (MXND)
C
      STDLEN = 0.
      NUML = 0
      S = 0.0
      DO 100 L = 1, LLL
         N1 = NXL (1, L)
         N2 = NXL (2, L)
         IF (N1 .GT. 0) THEN
            D = SQRT ((XN (N1) - XN (N2)) **2 + (YN (N1) - YN (N2)) **2)
            S = S + D
            NUML = NUML + 1
         ENDIF
  100 CONTINUE
C
      IF (NUML .LE. 0) RETURN
      TOL = 1.25 * S / FLOAT (NUML)
      NUML = 0
      S = 0.0
      DO 110 L = 1, LLL
         N1 = NXL (1, L)
         N2 = NXL (2, L)
         IF (N1 .GT. 0) THEN
            D = SQRT ((XN (N1) - XN (N2)) **2 + (YN (N1) - YN (N2)) **2)
            IF (D .LT. TOL) THEN
               S = S + D
               NUML = NUML + 1
            ENDIF
         ENDIF
  110 CONTINUE
      STDLEN = S / FLOAT (NUML)
C
      RETURN
C
      END
