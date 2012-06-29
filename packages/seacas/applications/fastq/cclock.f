C $Id: cclock.f,v 1.2 1991/03/21 15:44:23 gdsjaar Exp $
C $Log: cclock.f,v $
C Revision 1.2  1991/03/21 15:44:23  gdsjaar
C Changed all 3.14159... to atan2(0.0, -1.0)
C
c Revision 1.1.1.1  1990/11/30  11:04:21  gdsjaar
c FASTQ Version 2.0X
c
c Revision 1.1  90/11/30  11:04:19  gdsjaar
c Initial revision
c 
C
CC* FILE: [.QMESH]CCLOCK.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE CCLOCK (X, Y, N, CCW, ERR, INDETR)
C***********************************************************************
C
C  SUBROUTINE CCLOCK = DETERMINES IF THE PERIMETER OF A REGION IS STATED
C                      IN COUNTER-CLOCKWISE FASHION
C
C***********************************************************************
C
C  SUBROUTINE CALLED BY:
C     PERIM = GENERATES THE PERIMETER OF A REGION
C
C***********************************************************************
C
C  VARIABLES USED:
C     CCW    = .TRUE. IF THE PERIMETER IS IN COUNTER-CLOCKWISE ORDER
C     ERR    = .TRUE. IF THE ORDER COULD NOT BE DETERMINED, OR IF AN
C              ERROR OCCURS CHECKING THE ORDER
C     N      = THE NUMBER OF NODES IN THE PERIMETER (MUST BE AT LEAST 3)
C
C***********************************************************************
C
      DIMENSION X (N), Y (N)
C
      LOGICAL CCW, ERR, INDETR
C
      ERR = .TRUE.
      INDETR = .FALSE.
      PI = ATAN2(0.0, -1.0)
      TWOPI = PI+PI
C
      IF (N .LT. 3) THEN
         CALL MESAGE ('PERIMETER MUST CONTAIN MORE THAN 3 NODES')
         GOTO 110
      ENDIF
C
      SPIRO = 0.0
      IF ( (X (1) .EQ. X (N)) .AND. (Y (1) .EQ. Y (N)) ) THEN
         CALL MESAGE ('PERIMETER CONTAINS DUPLICATE NODE LOCATIONS')
         GOTO 110
      ENDIF
      AGOLD = ATAN2 (Y (1) - Y (N), X (1) - X (N))
      DO 100 I = 1, N
         NEXT = I + 1
         IF (NEXT .GT. N)NEXT = 1
         IF ( (X (NEXT) .EQ. X (I)) .AND. (Y (NEXT) .EQ. Y (I)) ) THEN
            CALL MESAGE ('PERIMETER CONTAINS DUPLICATE NODE LOCATIONS')
            GOTO 110
         ENDIF
         AGNEW = ATAN2 (Y (NEXT) - Y (I), X (NEXT) - X (I))
         DIFF = AGNEW - AGOLD
         IF (DIFF .GT. PI) DIFF = DIFF-TWOPI
         IF (DIFF .LT. -PI) DIFF = DIFF+TWOPI
         IF (ABS (ABS (DIFF) - PI) .LT. 1.0E-3) THEN
            CALL MESAGE ('PERIMETER CONTAINS SWITCHBACKS')
            GOTO 110
         ENDIF
         SPIRO = SPIRO + DIFF
         AGOLD = AGNEW
  100 CONTINUE
      CCW = .TRUE.
      IF (SPIRO .LT. 0.0) CCW = .FALSE.
      IF ( (ABS (SPIRO) .LT. PI) .OR. (ABS (SPIRO) .GT. (3.*PI))) THEN
         INDETR = .TRUE.
         RETURN
      ENDIF
      ERR = .FALSE.
      RETURN
C
C  ERROR IN THE ROUTINE
C
  110 CONTINUE
      RETURN
C
      END
