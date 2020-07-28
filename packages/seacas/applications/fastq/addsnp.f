C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE ADDSNP(MSNAP,SNAPDX,NSNAP,INDEX,VALUE,ERR)
C***********************************************************************

C   SUBROUTINE ADDSNP = ADDS SNAP GRID LINE DEFINITIONS

C***********************************************************************

C  VARIABLES USED:
C     MSNAP  = DIMENSION OV SNAP ARRAYS
C     SNAPDX = THE SNAP GRID VALUES ARRAY (X AND Y)
C     NSNAP  = THE NUMBER OF SNAP GRID VALUES IN X AND Y
C     INDEX  = 1 FOR X VALUES, 2 FOR Y VALUES
C     VALUE  = THE GRID VALUE TO BE ADDED
C     KOUNT  = THE LOCATION OF THE SNAPDX VALUE JUST LESS THAN VALUE

C***********************************************************************

      DIMENSION SNAPDX(2,MSNAP),NSNAP(2)

      LOGICAL ERR

C  ADD THE SNAP GRID VALUE WHERE IT FITS IN NUMERICAL ORDER

      ERR=.FALSE.

      IF(NSNAP(INDEX).GT.0)THEN
         KOUNT=0
         DO 100 I=1,NSNAP(INDEX)
            IF(VALUE.LT.SNAPDX(INDEX,I))GOTO 110
            KOUNT=I

C  IF THE VALUE IS ALREADY THERE, THEN DON'T ADD IT AGAIN - JUST RETURN

            IF(VALUE.EQ.SNAPDX(INDEX,I))RETURN

  100    CONTINUE
  110    CONTINUE
         IF(NSNAP(INDEX).EQ.MSNAP)THEN
            CALL MESAGE('** NO MORE ROOM FOR ADDITIONAL GRID LINES **')
            WRITE(*,10000)MSNAP
            ERR=.TRUE.
            RETURN
         ENDIF
         NSNAP(INDEX)=NSNAP(INDEX)+1
         DO 120 I=NSNAP(INDEX),KOUNT+2,-1
            SNAPDX(INDEX,I)=SNAPDX(INDEX,I-1)
  120    CONTINUE
         SNAPDX(INDEX,KOUNT+1)=VALUE

C  JUST PUT THE FIRST VALUE WHERE IT BELONGS

      ELSE
         NSNAP(INDEX)=1
         SNAPDX(INDEX,1)=VALUE
      ENDIF

      RETURN

10000 FORMAT(' THE MAXIMUM NUMBER OF GRID LINES IS: ',I10)

      END
