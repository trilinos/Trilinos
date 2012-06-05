C $Id: bcross.f,v 1.3 2004/01/26 17:28:18 gdsjaar Exp $
C $Log: bcross.f,v $
C Revision 1.3  2004/01/26 17:28:18  gdsjaar
C Removed several unused variables from getang subroutine.
C
C Initialized a variable
C
C Revision 1.2  2004/01/22 14:25:22  gdsjaar
C Attempt to fix strange problem on x86_64 AMD Opteron system using
C Portland Group 5.1-3 compilers. The getang function would work
C correctly if compiled with no optimization and in debug mode, but
C would crash if compiled optimized. The location of the crash was not
C in a place that made any sense that something was wrong.
C
C After much trial and error, it was found that adding a 'SAVE'
C statement at the beginning of the file fixed the problem.
C
C Also cleaned out some unused parameters being passed to the function.
C
C Revision 1.1.1.1  1990/11/30 11:03:55  gdsjaar
C FASTQ Version 2.0X
C
c Revision 1.1  90/11/30  11:03:52  gdsjaar
c Initial revision
c 
C
CC* FILE: [.PAVING]BCROSS.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE BCROSS (MXND, MLN, XN, YN, ZN, LXK, KXL, NXL, LXN,
     &   LNODES, I1, I2, J1, J2, NLOOP, BOK, LLL, XMIN, XMAX, YMIN,
     &   YMAX, ZMIN, ZMAX, DEV1, KREG, ERR)
C***********************************************************************
C
C  SUBROUTINE BCROSS = CHECKS TO MAKE SURE THAT A BOUNDARY IS NOT
C                      BECOMING A PERMANENT CROSS
C
C***********************************************************************
C
      DIMENSION XN(MXND), YN(MXND), ZN(MXND)
      DIMENSION NXL(2, 3*MXND), LXN(4, MXND)
      DIMENSION LXK(4, MXND), KXL(2, 3*MXND)
      DIMENSION LNODES(MLN, MXND)
C
      LOGICAL BOK, LCROSS, ERR
C
      CHARACTER*3 DEV1
C
      BOK = .TRUE.
      ERR = .FALSE.
C
      J0 = LNODES(2, J1)
      J3 = LNODES(3, J2)
C
C  IF J0 TO I2, OR J3 TO I1 IS TO BECOME A BOUNDARY LINE,
C  THEN TEST THE CONNECTION TO SEE IF IT INTERSECTS ANY OTHER
C  BOUNDARY LINES
C
      KOUNT = 0
C
C  FIRST TEST THE J0 TO I2 LINE
C
      IF ((LXN(2, J0) .LT. 0) .AND. (LXN(2, I2) .LT. 0)) THEN
         NODE1 = I1
         NODE2 = I2
  100    CONTINUE
         NODE1 = NODE2
         NODE2 = LNODES(3, NODE2)
         IF ((LXN(2, NODE1) .LT. 0) .AND. (LXN(2, NODE2) .LT. 0)) THEN
            IF (NODE2 .EQ. J0) THEN
               CALL GETANG (MXND, MLN, XN, YN, LNODES, LXK, KXL,
     &            NXL, LXN, NODE1, J0, J1, ANGLE, ERR)
               IF (ANGLE .LT. 0) THEN
                  BOK = .FALSE.
                  GOTO 130
               ELSE
                  GOTO 110
               ENDIF
C
            ELSE
               CALL INTSCT (XN(NODE1), YN(NODE1), XN(NODE2), YN(NODE2),
     &            XN(J0), YN(J0), XN(I2), YN(I2), U, W, LCROSS)
               IF (LCROSS) THEN
                  BOK = .FALSE.
                  GOTO 130
               ENDIF
            ENDIF
         ENDIF
         KOUNT = KOUNT + 1
         IF (KOUNT .LT. NLOOP) THEN
            GOTO 100
         ELSE
            ERR = .TRUE.
            GOTO 130
         ENDIF
      ENDIF
C
  110 CONTINUE
C
C  NEXT TEST THE J3 TO I1 LINE
C
      KOUNT = 0
      IF ((LXN(2, J3) .LT. 0) .AND. (LXN(2, I1) .LT. 0)) THEN
         NODE1 = J3
         NODE2 = LNODES(3, J3)
  120    CONTINUE
         NODE1 = NODE2
         NODE2 = LNODES(3, NODE2)
         IF ((LXN(2, NODE1) .LT. 0) .AND. (LXN(2, NODE2) .LT. 0)) THEN
            IF (NODE2 .EQ. I1) THEN
               CALL GETANG (MXND, MLN, XN, YN, LNODES, LXK, KXL,
     &            NXL, LXN, NODE1, I1, I2, ANGLE, ERR)
               IF (ANGLE .LT. 0) THEN
                  BOK = .FALSE.
                  GOTO 130
               ELSE
                  GOTO 130
               ENDIF
C
            ELSE
               CALL INTSCT (XN(NODE1), YN(NODE1), XN(NODE2), YN(NODE2),
     &            XN(J3), YN(J3), XN(I1), YN(I1), U, W, LCROSS)
               IF (LCROSS) THEN
                  BOK = .FALSE.
                  GOTO 130
               ENDIF
            ENDIF
         ENDIF
         KOUNT = KOUNT + 1
         IF (KOUNT .LT. NLOOP) THEN
            GOTO 120
         ELSE
            ERR = .TRUE.
            GOTO 130
         ENDIF
      ENDIF
C
  130 CONTINUE
C
      RETURN
C
      END
