C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE CLOSEP (MP, N15, X, Y, IPOINT, COOR, LINKP, JJ)
C***********************************************************************

C  SUBROUTINE CLOSE = FINDS THE CLOSEST EXISTING POINT TO THE MOUSE

C***********************************************************************

C  SUBROUTINE CALLED BY:
C     INPUT  = INPUTS MESH DEFINITIONS FROM THE LIGHT TABLE

C***********************************************************************

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

C***********************************************************************

      DIMENSION IPOINT (MP), COOR (2, MP), LINKP (2, MP)

      LOGICAL ADDLNK

      ADDLNK = .FALSE.
      DMIN = 100000.

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
