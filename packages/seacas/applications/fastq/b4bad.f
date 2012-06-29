C $Id: b4bad.f,v 1.2 1998/07/14 18:18:23 gdsjaar Exp $
C $Log: b4bad.f,v $
C Revision 1.2  1998/07/14 18:18:23  gdsjaar
C Removed unused variables, cleaned up a little.
C
C Changed BLUE labels to GREEN to help visibility on black background
C (indirectly requested by a couple users)
C
C Revision 1.1.1.1  1990/11/30 11:03:50  gdsjaar
C FASTQ Version 2.0X
C
c Revision 1.1  90/11/30  11:03:49  gdsjaar
c Initial revision
c 
C
CC* FILE: [.PAVING]B4BAD.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE B4BAD (MXND, MLN, XN, YN, LXK, KXL, NXL, LXN, LNODES,
     &   ANGLE, I1, I2, J1, J2, NLOOP, KOUNTL, BOK, ERR)
C***********************************************************************
C
C  SUBROUTINE BCROSS = CHECKS TO MAKE SURE THAT A BOUNDARY IS NOT
C                      BECOMING A PERMANENT CROSS
C
C***********************************************************************
C
      DIMENSION NXL (2, 3*MXND), LXN (4, MXND)
      DIMENSION LXK (4, MXND), KXL (2, 3*MXND)
      DIMENSION LNODES (MLN, MXND), ANGLE (MXND)
      DIMENSION NODE (4)
C
      LOGICAL BOK, ERR
C
      BOK = .TRUE.
      ERR = .FALSE.
C
C  GET THE NODES THAT FORM THE REMAINING ELEMENT
C
      IF (KOUNTL .EQ. 4) THEN
         NODE(1) = I2
         NODE(2) = LNODES (3, NODE(1))
         NODE(3) = LNODES (3, NODE(2))
         NODE(4) = LNODES (3, NODE(3))
      ELSEIF (NLOOP - KOUNTL - 2 .EQ. 4) THEN
         NODE(1) = I1
         NODE(2) = LNODES (3, J2)
         NODE(3) = LNODES (3, NODE(2))
         NODE(4) = LNODES (3, NODE(3))
      ELSE
         GOTO 110
      ENDIF
C
C  NOW CHECK ALL THE NODES TO SEE IF THEY ARE ON THE BOUNDARY
C  AND CAN BE CLASSIFIED AS CORNERS
C
      DO 100 I = 1, 4
         IF ( (LXN (2, NODE (I)) .LT. 0) .AND.
     &      (LNODES (6, NODE (I)) .GE. 3) ) THEN
            BOK = .FALSE.
            GOTO 110
         ENDIF
  100 CONTINUE
C
  110 CONTINUE
C
      RETURN
C
      END
