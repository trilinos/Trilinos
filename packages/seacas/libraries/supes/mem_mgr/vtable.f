C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details
      SUBROUTINE VTABLE (NEWLOC, NEWLEN, VOID, LVOID, NVOIDS, CHRCOL,
     *   ERR)
      IMPLICIT INTEGER (A-Z)

C     THIS SUBROUTINE INSERTS NEW VOIDS IN THE VOID TABLE AND
C     THEN CHECKS FOR CONTIGUOUS VOIDS WHICH ARE THEN JOINED.

C     ERROR CODES

C     ERROR VECTOR AND FLAGS.
C     THE ERROR PARAMETERS BELONG IN MDINIT ALSO.

      INCLUDE 'params.inc'

C     VFULL  = NO ROOM IN VOID TABLE
C     BDVOID = OVERLAPPING VOIDS

      DIMENSION VOID(LVOID,CHRCOL,2)

      IF (NEWLEN .GT. 0) THEN

         IF (NVOIDS .GE. LVOID) THEN
            ERR = VFULL
            RETURN
         END IF

C        FIND LOCATION FOR NEW ENTRY.

         CALL SRCHI(VOID,1,NVOIDS,NEWLOC,ERR,ROW)
         IF (ERR .NE. 0) THEN
            ERR = BDVOID
            RETURN
         END IF

C        NEW ENTRY IN TABLE.

         IF (ROW .LE. NVOIDS) THEN

C           MAKE ROOM FOR NEW ENTRY.

            CALL SHFTI (VOID, LVOID*CHRCOL, 2, ROW, NVOIDS, -1)

         END IF

         VOID(ROW,1,1) = NEWLOC
         VOID(ROW,1,2) = NEWLEN
         NVOIDS = NVOIDS + 1

      END IF

C     CHECK TABLE TO SEE IF ANY VOIDS HAVE JOINED OR ARE ZERO LENGTH.

C     NOTE THAT A STANDARD DO LOOP CANNOT BE USED BECAUSE THE UPPER
C     LIMIT OF THE LOOP CAN CHANGE INSIDE THE LOOP.

      I = 1
  100 IF (I .GE. NVOIDS) GO TO 110
         IF (VOID(I,1,1)+VOID(I,1,2) .EQ. VOID(I+1,1,1)) THEN

C           THESE TWO VOIDS SHOULD BE JOINED.

            VOID(I,1,2) = VOID(I,1,2) + VOID(I+1,1,2)
            CALL SHFTI (VOID, LVOID*CHRCOL, 2, I+2, NVOIDS, 1)
            NVOIDS = NVOIDS - 1
            GO TO 100

         ELSE IF (VOID(I,1,2) .EQ. 0) THEN

C           THIS VOID IS ZERO LENGTH.

            CALL SHFTI (VOID, LVOID*CHRCOL, 2, I+1, NVOIDS, 1)
            NVOIDS = NVOIDS - 1

         ELSE IF (VOID(I,1,1)+VOID(I,1,2) .GT. VOID(I+1,1,1)) THEN

C           OVERLAPPING VOIDS

            ERR = BDVOID
            RETURN

         END IF

         I = I + 1
         GO TO 100

  110 CONTINUE

C     CHECK LAST VOID

      IF (NVOIDS .GE. 1) THEN
         IF (VOID(NVOIDS,1,2) .EQ. 0) NVOIDS = NVOIDS - 1
      END IF

      ERR = SUCESS
      RETURN
      END
