C $Id: elipse.f,v 1.2 1991/03/21 15:44:37 gdsjaar Exp $
C $Log: elipse.f,v $
C Revision 1.2  1991/03/21 15:44:37  gdsjaar
C Changed all 3.14159... to atan2(0.0, -1.0)
C
c Revision 1.1.1.1  1990/11/30  11:06:32  gdsjaar
c FASTQ Version 2.0X
c
c Revision 1.1  90/11/30  11:06:31  gdsjaar
c Initial revision
c 
C
CC* FILE: [.MAIN]ELIPSE.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      REAL FUNCTION ELIPSE (A7, A8, A2, ANG)
C***********************************************************************
C
C  FUNCTION ELIPSE = CALCULATES THE ANGULAR EQUATION ERROR WHEN FINDING
C                    AN ELIPSE PATTERN
C
C***********************************************************************
C
      PI = ATAN2(0.0, -1.0)
      A4 = A8 - ANG
      A5 = A7 - ANG
      A3 = A2 - A4
      A6 = PI - A5 - A2
      ELIPSE = SIN(A4) * SIN(A6) - SIN(A5) * SIN (A3)
      RETURN
C
      END
