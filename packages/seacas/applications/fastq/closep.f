C $Id: closep.f,v 1.1 1990/11/30 11:05:03 gdsjaar Exp $
C $Log: closep.f,v $
C Revision 1.1  1990/11/30 11:05:03  gdsjaar
C Initial revision
C
C
CC* FILE: [.MAIN]CLOSEP.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE CLOSEP (MP, N15, X, Y, IPOINT, COOR, LINKP, JJ)
C***********************************************************************
C
C  SUBROUTINE CLOSE = FINDS THE CLOSEST EXISTING POINT TO THE MOUSE
C
C***********************************************************************
C
C  SUBROUTINE CALLED BY:
C     INPUT  = INPUTS MESH DEFINITIONS FROM THE LIGHT TABLE
C
C***********************************************************************
C
C  VARIABLES USED:
C     X      = THE X LOCATION IN USER COORDINATES
C     Y      = THE Y LOCATION IN USER COORDINATES
C     POINT  = ARRAY OF VALUES DEFINING A POINT
C               (I, 1) = THE NUMBER OF THE POINT
C               (I, 2) = THE X COORDINATE OF THE POINT
C               (I, 3) = THE Y COORDINATE OF THE POINT
C               (I, 4) = THE BOUNDARY FLAG OF THE POINT
C     I      = THE NUMBER OF THE CLOSEST POINT FOUND
C     K      = THE NUMBER OF POINTS IN THE DATABASE
C
C***********************************************************************
C
      DIMENSION IPOINT (MP), COOR (2, MP), LINKP (2, MP)
C
      LOGICAL ADDLNK
C
      ADDLNK = .FALSE.
      DMIN = 100000.
C
      DO 100 I = 1, N15
         CALL LTSORT (MP, LINKP, I, IPNTR, ADDLNK)
         IF (IPNTR .GT. 0) THEN
            DIST = SQRT (((COOR (1, IPNTR) - X) **2) +
     &         ((COOR (2, IPNTR) - Y) **2))
            IF (DIST .LT. DMIN) THEN
               DMIN = DIST
               JJ = IPOINT (IPNTR)
            ENDIF
         ENDIF
  100 CONTINUE
      RETURN
      END
