C $Id: invert.f,v 1.3 2004/01/26 17:28:18 gdsjaar Exp $
C $Log: invert.f,v $
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
C Revision 1.1.1.1  1990/11/30 11:10:24  gdsjaar
C FASTQ Version 2.0X
C
c Revision 1.1  90/11/30  11:10:23  gdsjaar
c Initial revision
c 
C
CC* FILE: [.PAVING]INVERT.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE INVERT (MXND, MLN, XN, YN, ZN, LXK, KXL, NXL, LXN, LLL,
     &   LNODES, XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX, DEV1, KREG, NODE,
     &   XDEL, YDEL)
C***********************************************************************
C
C  SUBROUTINE INVERT = CHECKS FOR AN INVERSION OR CROSSING OF A BOUNDARY
C                      UPON ITSELF AND CORRECTS IT WHERE NECESSARY
C
C***********************************************************************
C
      DIMENSION XN(MXND), YN(MXND), ZN(MXND)
      DIMENSION LXN(4, MXND), NXL(2, 3*MXND)
      DIMENSION LXK(4, MXND), KXL(2, 3*MXND)
      DIMENSION LNODES (MLN, MXND)
C
      LOGICAL ERR, VCROSS
C
      CHARACTER*3 DEV1
C
      ERR = .FALSE.

      XOLD = XN (NODE)
      YOLD = YN (NODE)
C
      N2 = NODE
      N3 = LNODES (3, N2)
      N4 = LNODES (3, N3)
      N1 = LNODES (2, N2)
      N0 = LNODES (2, N1)
C
C  GET THE ANGLES BEFORE MOVEMENT
C
      IF (LXN (4, N1) .EQ. 0)
     &   CALL GETANG (MXND, MLN, XN, YN, LNODES, LXK, KXL, NXL,
     &   LXN, N0, N1, N2, ANG1A, ERR)
      IF (LXN (4, N2) .EQ. 0)
     &   CALL GETANG (MXND, MLN, XN, YN, LNODES, LXK, KXL, NXL,
     &   LXN, N1, N2, N3, ANG2A, ERR)
      IF (LXN (4, N3) .EQ. 0)
     &   CALL GETANG (MXND, MLN, XN, YN, LNODES, LXK, KXL, NXL,
     &   LXN, N2, N3, N4, ANG3A, ERR)
C
C  NOW PLACE THE NODE TEMPORARILY AT THE NEW PROPOSED LOCATION
C
      XN (NODE) = XN (NODE) + XDEL
      YN (NODE) = YN (NODE) + YDEL
C
C  GET THE ANGLE BEING ADJUSTED AT THE NODE ITSELF
C
      IF ((LXN (4, N2) .EQ. 0) .AND. (ANG2A .GT. 0.)) THEN
         CALL GETANG (MXND, MLN, XN, YN, LNODES, LXK, KXL, NXL,
     &      LXN, N1, N2, N3, ANG2B, ERR)
C
C  ADJUST THE NODE LOCATION IF NECESSARY
C
         IF (ANG2B .LT. 0.) THEN
            CALL VINTER (MXND, XN, YN, N1, N3, N2, XOLD, YOLD,
     &         XNEW, YNEW, VCROSS)
            IF (VCROSS) THEN
               XN (NODE) = XNEW
               YN (NODE) = YNEW
            ENDIF
         ENDIF
      ENDIF
C
C  GET THE ANGLE BEING ADJUSTED ON THE CCW SIDE OF THIS NODE
C
      IF ((LXN (4, N1) .EQ. 0) .AND. (ANG1A .GT. 0.)) THEN
         CALL GETANG (MXND, MLN, XN, YN, LNODES, LXK, KXL, NXL,
     &      LXN, N0, N1, N2, ANG1B, ERR)
C
C  ADJUST THE NODE LOCATION IF NECESSARY
C
         IF (ANG1B .LT. 0.) THEN
            CALL VINTER (MXND, XN, YN, N1, N0, N2, XOLD, YOLD,
     &         XNEW, YNEW, VCROSS)
            IF (VCROSS) THEN
               XN (NODE) = XNEW
               YN (NODE) = YNEW
            ENDIF
         ENDIF
      ENDIF
C
C  GET THE ANGLE BEING ADJUSTED ON THE CW SIDE OF THIS NODE
C
      IF ((LXN (4, N3) .EQ. 0) .AND. (ANG3A .GT. 0.)) THEN
         CALL GETANG (MXND, MLN, XN, YN, LNODES, LXK, KXL, NXL,
     &      LXN, N2, N3, N4, ANG3B, ERR)
C
C  ADJUST THE NODE LOCATION IF NECESSARY
C
         IF (ANG3B .LT. 0.) THEN
            CALL VINTER (MXND, XN, YN, N3, N4, N2, XOLD, YOLD,
     &         XNEW, YNEW, VCROSS)
            IF (VCROSS) THEN
               XN (NODE) = XNEW
               YN (NODE) = YNEW
            ENDIF
         ENDIF
      ENDIF
C
C  RESTORE THE OLD LOCATION AND THE XDEL AND YDEL TO THE CORRECTED
C  VALUES
C
      XDEL = XN (NODE) - XOLD
      YDEL = YN (NODE) - YOLD
      XN (NODE) = XOLD
      YN (NODE) = YOLD
C
      RETURN
C
      END
