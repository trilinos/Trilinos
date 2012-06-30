C $Id: close6.f,v 1.2 1998/07/14 18:18:31 gdsjaar Exp $
C $Log: close6.f,v $
C Revision 1.2  1998/07/14 18:18:31  gdsjaar
C Removed unused variables, cleaned up a little.
C
C Changed BLUE labels to GREEN to help visibility on black background
C (indirectly requested by a couple users)
C
C Revision 1.1.1.1  1990/11/30 11:04:54  gdsjaar
C FASTQ Version 2.0X
C
c Revision 1.1  90/11/30  11:04:53  gdsjaar
c Initial revision
c 
C
CC* FILE: [.PAVING]CLOSE6.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/31/90
CC* MODIFICATION: ADDED ARGUMENTS TO CALL TO CLOSE6 TO PASS MINIMUM
CC**              ELEMENT SIZE (SIZMIN) AND GETSIZ PARAMETERS OF
CC**              EMIN AND EMAX
C
      SUBROUTINE CLOSE6 (MXND, MXCORN, MLN, NUID, XN, YN, LXK, KXL, NXL,
     &   LXN, ANGLE, BNSIZE, LNODES, NODE, NLOOP, KKKOLD, LLLOLD,
     &   NNNOLD, NAVAIL, IAVAIL, DONE, XMIN, XMAX, YMIN, YMAX, DEV1,
     &   LLL, KKK, NNN, LCORN, NCORN, GRAPH, VIDEO, SIZEIT, NOROOM,
     &   ERR, XNOLD, YNOLD, NXKOLD, LINKEG, LISTEG, BMESUR, MLINK,
     &   NPNOLD, NPEOLD, NNXK, REMESH, REXMIN, REXMAX, REYMIN, REYMAX,
     &   IDIVIS, SIZMIN, EMAX, EMIN)
C***********************************************************************
C
C  SUBROUTINE CLOSE6 = FINISHES UP A LOOP WITH ONLY 6 LINES IN THE LOOP
C
C***********************************************************************
C
      DIMENSION XN (MXND), YN (MXND), NUID (MXND)
      DIMENSION LXK (4, MXND), KXL (2, 3*MXND)
      DIMENSION NXL (2, 3*MXND), LXN (4, MXND)
      DIMENSION ANGLE (MXND), LNODES (MLN, MXND), BNSIZE (2, MXND)
      DIMENSION LCORN (MXCORN)
C
      DIMENSION XNOLD(NPNOLD), YNOLD(NPNOLD)
      DIMENSION NXKOLD(NNXK, NPEOLD)
      DIMENSION LINKEG(2, MLINK), LISTEG(4 * NPEOLD), BMESUR(NPNOLD)
C
      LOGICAL GRAPH, VIDEO, ERR, SIZEIT, NOROOM
C
      CHARACTER*3 DEV1
C
      ERR = .FALSE.
C
C  SET ALL THE LOOP NODES TO BE INTERIOR
C
      INODE = NODE
      IKOUNT = 0
  100 CONTINUE
      IKOUNT = IKOUNT + 1
      IF (IKOUNT .GT. 6) THEN
         CALL MESAGE ('** PROBLEMS IN CLOSE6 WITH TOO MANY IN LOOP **')
         ERR = .TRUE.
         GOTO 110
      ENDIF
C
      LNODES (4, INODE) = - 2
      CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES, INODE, ERR)
      IF (ERR) GOTO 110
C
      INODE = LNODES (3, INODE)
      IF (INODE .NE. NODE) GOTO 100
C
C  NOW GET THE CORRECT INTERPRETATION OF THE SHAPE
C
      CALL CNTCRN (MXND, MXCORN, MLN, LNODES, LCORN, NCORN, NLOOP,
     &   NODE, ERR)
      IF (ERR) GOTO 110
C
C  PROCESS A TRIANGLE SHAPE WITH VARIOUS INTERVAL COMBINATIONS
C
      IF (NCORN .EQ. 3) THEN
         I1 = LNODES (7, LCORN(1))
         I2 = LNODES (7, LCORN(2))
         I3 = LNODES (7, LCORN(3))
C
C  HANDLE A 4-1-1 TRIANGLE
C
         IF (MAX0 (I1, I2, I3) .EQ. 4) THEN
            IF (I1 .EQ. 4) THEN
               XNEW = ( ( (XN (LCORN (1)) + XN (LCORN (2)) ) * .5) +
     &            XN (LCORN (3)) ) * .5
               YNEW = ( ( (YN (LCORN (1)) + YN (LCORN (2)) ) * .5) +
     &            YN (LCORN (3)) ) * .5
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/31/90
CC* MODIFICATION: ADDED ARGUMENTS TO CALL TO ADD1CN TO PASS MINIMUM
CC**              ELEMENT SIZE (SIZMIN) AND GETSIZ PARAMETERS OF
CC**              EMIN AND EMAX
C
               CALL ADD1CN (MXND, MLN, XN, YN, NUID, LXK, KXL, NXL, LXN,
     &            ANGLE, BNSIZE, LNODES, NNN, KKK, LLL, NNNOLD, LLLOLD,
     &            NLOOP, XNEW, YNEW, LCORN(1), IAVAIL, NAVAIL,
     &            GRAPH, VIDEO, SIZEIT, NOROOM, ERR, XNOLD, YNOLD,
     &            NXKOLD, LINKEG, LISTEG, BMESUR, MLINK, NPNOLD,
     &            NPEOLD, NNXK, REMESH, REXMIN, REXMAX, REYMIN, REYMAX,
     &            IDIVIS, SIZMIN, EMAX, EMIN)
               IF ((NOROOM) .OR. (ERR)) GOTO 110
            ELSEIF (I2 .EQ. 4) THEN
               XNEW = ( ( (XN (LCORN (2)) + XN (LCORN (3)) ) * .5) +
     &            XN (LCORN (1)) ) * .5
               YNEW = ( ( (YN (LCORN (2)) + YN (LCORN (3)) ) * .5) +
     &            YN (LCORN (1)) ) * .5
               CALL ADD1CN (MXND, MLN, XN, YN, NUID, LXK, KXL, NXL, LXN,
     &            ANGLE, BNSIZE, LNODES, NNN, KKK, LLL, NNNOLD, LLLOLD,
     &            NLOOP, XNEW, YNEW, LCORN(2), IAVAIL, NAVAIL,
     &            GRAPH, VIDEO, SIZEIT, NOROOM, ERR, XNOLD, YNOLD,
     &            NXKOLD, LINKEG, LISTEG, BMESUR, MLINK, NPNOLD,
     &            NPEOLD, NNXK, REMESH, REXMIN, REXMAX, REYMIN, REYMAX,
     &            IDIVIS, SIZMIN, EMAX, EMIN)
               IF ((NOROOM) .OR. (ERR)) GOTO 110
            ELSE
               XNEW = ( ( (XN (LCORN (3)) + XN (LCORN (1)) ) * .5) +
     &            XN (LCORN (2)) ) * .5
               YNEW = ( ( (YN (LCORN (3)) + YN (LCORN (1)) ) * .5) +
     &            YN (LCORN (2)) ) * .5
               CALL ADD1CN (MXND, MLN, XN, YN, NUID, LXK, KXL, NXL, LXN,
     &            ANGLE, BNSIZE, LNODES, NNN, KKK, LLL, NNNOLD, LLLOLD,
     &            NLOOP, XNEW, YNEW, LCORN(3), IAVAIL, NAVAIL,
     &            GRAPH, VIDEO, SIZEIT, NOROOM, ERR, XNOLD, YNOLD,
     &            NXKOLD, LINKEG, LISTEG, BMESUR, MLINK, NPNOLD,
     &            NPEOLD, NNXK, REMESH, REXMIN, REXMAX, REYMIN, REYMAX,
     &            IDIVIS, SIZMIN, EMAX, EMIN)
               IF ((NOROOM) .OR. (ERR)) GOTO 110
            ENDIF
C
C  HANDLE A 3-2-1 TRIANGLE
C
         ELSEIF (MAX0 (I1, I2, I3) .EQ. 3) THEN
            IF (I1 .EQ. 1) THEN
               CALL ADD2EL (MXND, MLN, XN, YN, NUID, LXK, KXL, NXL, LXN,
     &            ANGLE, LNODES, NNN, KKK, LLL, NNNOLD, LLLOLD,
     &            NLOOP, LNODES (2, LCORN(1)),
     &            LNODES (3, LCORN(2)), IAVAIL, NAVAIL, GRAPH, VIDEO,
     &            NOROOM, ERR)
               IF ((NOROOM) .OR. (ERR)) GOTO 110
            ELSEIF (I2 .EQ. 1) THEN
               CALL ADD2EL (MXND, MLN, XN, YN, NUID, LXK, KXL, NXL, LXN,
     &            ANGLE, LNODES, NNN, KKK, LLL, NNNOLD, LLLOLD,
     &            NLOOP, LNODES (2, LCORN(2)),
     &            LNODES (3, LCORN(3)), IAVAIL, NAVAIL, GRAPH, VIDEO,
     &            NOROOM, ERR)
               IF ((NOROOM) .OR. (ERR)) GOTO 110
            ELSE
               CALL ADD2EL (MXND, MLN, XN, YN, NUID, LXK, KXL, NXL, LXN,
     &            ANGLE, LNODES, NNN, KKK, LLL, NNNOLD, LLLOLD,
     &            NLOOP, LNODES (2, LCORN(3)),
     &            LNODES (3, LCORN(1)), IAVAIL, NAVAIL, GRAPH, VIDEO,
     &            NOROOM, ERR)
               IF ((NOROOM) .OR. (ERR)) GOTO 110
            ENDIF
C
C  HANDLE A 2-2-2 TRIANGLE
C
         ELSE
            XNEW1 = ( ( (XN (LCORN (1)) + XN (LCORN (2)) ) * .5) +
     &         XN (LCORN (3)) ) * .5
            XNEW2 = ( ( (XN (LCORN (2)) + XN (LCORN (3)) ) * .5) +
     &         XN (LCORN (1)) ) * .5
            XNEW3 = ( ( (XN (LCORN (3)) + XN (LCORN (1)) ) * .5) +
     &         XN (LCORN (2)) ) * .5
C
            YNEW1 = ( ( (YN (LCORN (1)) + YN (LCORN (2)) ) * .5) +
     &         YN (LCORN (3)) ) * .5
            YNEW2 = ( ( (YN (LCORN (2)) + YN (LCORN (3)) ) * .5) +
     &         YN (LCORN (1)) ) * .5
            YNEW3 = ( ( (YN (LCORN (3)) + YN (LCORN (1)) ) * .5) +
     &         YN (LCORN (2)) ) * .5
C
            XNEW = (XNEW1 + XNEW2 + XNEW3) / 3.
            YNEW = (YNEW1 + YNEW2 + YNEW3) / 3.
            CALL ADD1CN (MXND, MLN, XN, YN, NUID, LXK, KXL, NXL, LXN,
     &         ANGLE, BNSIZE, LNODES, NNN, KKK, LLL, NNNOLD, LLLOLD,
     &         NLOOP, XNEW, YNEW, LCORN(1), IAVAIL, NAVAIL,
     &         GRAPH, VIDEO, SIZEIT, NOROOM, ERR, XNOLD, YNOLD, NXKOLD,
     &         LINKEG, LISTEG, BMESUR, MLINK, NPNOLD, NPEOLD, NNXK,
     &         REMESH, REXMIN, REXMAX, REYMIN, REYMAX, IDIVIS, SIZMIN,
     &         EMAX, EMIN)
            IF ((NOROOM) .OR. (ERR)) GOTO 110
         ENDIF
C
C  PROCESS A RECTANGLE SHAPE WITH VARIOUS INTERVAL COMBINATIONS
C
      ELSEIF (NCORN .EQ. 4) THEN
         I1 = LNODES (7, LCORN(1))
         I2 = LNODES (7, LCORN(2))
         I3 = LNODES (7, LCORN(3))
         I4 = LNODES (7, LCORN(4))
C
C  HANDLE A 3-1-1-1 RECTANGLE
C
         IF (MAX0 (I1, I2, I3, I4) .EQ. 3) THEN
            IF (I1 .EQ. 3) THEN
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/31/90
CC* MODIFICATION: ADDED ARGUMENTS TO CALL TO ADD2CN TO PASS MINIMUM
CC**              ELEMENT SIZE (SIZMIN) AND GETSIZ PARAMETERS OF
CC**              EMIN AND EMAX
C
               CALL ADD2CN (MXND, MLN, XN, YN, NUID, LXK, KXL, NXL, LXN,
     &            ANGLE, BNSIZE, LNODES, NNN, KKK, LLL, NNNOLD, LLLOLD,
     &            NLOOP, LCORN(1), IAVAIL, NAVAIL, GRAPH, VIDEO,
     &            SIZEIT, NOROOM, ERR, XNOLD, YNOLD, NXKOLD, LINKEG,
     &            LISTEG, BMESUR, MLINK, NPNOLD, NPEOLD, NNXK, REMESH,
     &            REXMIN, REXMAX, REYMIN, REYMAX, IDIVIS, SIZMIN, EMAX,
     &            EMIN)
               IF ((NOROOM) .OR. (ERR)) GOTO 110
            ELSEIF (I2 .EQ. 3) THEN
               CALL ADD2CN (MXND, MLN, XN, YN, NUID, LXK, KXL, NXL, LXN,
     &            ANGLE, BNSIZE, LNODES, NNN, KKK, LLL, NNNOLD, LLLOLD,
     &            NLOOP, LCORN(2), IAVAIL, NAVAIL, GRAPH, VIDEO,
     &            SIZEIT, NOROOM, ERR, XNOLD, YNOLD, NXKOLD, LINKEG,
     &            LISTEG, BMESUR, MLINK, NPNOLD, NPEOLD, NNXK, REMESH,
     &            REXMIN, REXMAX, REYMIN, REYMAX, IDIVIS, SIZMIN, EMAX,
     &            EMIN)
               IF ((NOROOM) .OR. (ERR)) GOTO 110
            ELSEIF (I3 .EQ. 3) THEN
               CALL ADD2CN (MXND, MLN, XN, YN, NUID, LXK, KXL, NXL, LXN,
     &            ANGLE, BNSIZE, LNODES, NNN, KKK, LLL, NNNOLD, LLLOLD,
     &            NLOOP, LCORN(3), IAVAIL, NAVAIL, GRAPH, VIDEO,
     &            SIZEIT, NOROOM, ERR, XNOLD, YNOLD, NXKOLD, LINKEG,
     &            LISTEG, BMESUR, MLINK, NPNOLD, NPEOLD, NNXK, REMESH,
     &            REXMIN, REXMAX, REYMIN, REYMAX, IDIVIS, SIZMIN, EMAX,
     &            EMIN)
               IF ((NOROOM) .OR. (ERR)) GOTO 110
            ELSE
               CALL ADD2CN (MXND, MLN, XN, YN, NUID, LXK, KXL, NXL, LXN,
     &            ANGLE, BNSIZE, LNODES, NNN, KKK, LLL, NNNOLD, LLLOLD,
     &            NLOOP, LCORN(4), IAVAIL, NAVAIL, GRAPH, VIDEO,
     &            SIZEIT, NOROOM, ERR, XNOLD, YNOLD, NXKOLD, LINKEG,
     &            LISTEG, BMESUR, MLINK, NPNOLD, NPEOLD, NNXK, REMESH,
     &            REXMIN, REXMAX, REYMIN, REYMAX, IDIVIS, SIZMIN, EMAX,
     &            EMIN)
               IF ((NOROOM) .OR. (ERR)) GOTO 110
            ENDIF
C
C  HANDLE A 2-2-1-1 RECTANGLE
C
         ELSEIF (MAX0 ((I1+I2), (I2+I3), (I3+I4), (I4+I1)) .EQ. 4) THEN
            IF ( (I1+I2) .EQ. 4) THEN
               XNEW = ( XN (LNODES (3, LCORN(1))) +
     &            XN (LNODES (3, LCORN(2))) +
     &            ( XN (LCORN (3)) * .5) +
     &            XN (LCORN (4)) +
     &            ( XN (LCORN (1)) * .5) ) * .25
               YNEW = ( YN (LNODES (3, LCORN(1))) +
     &            YN (LNODES (3, LCORN(2))) +
     &            ( YN (LCORN (3)) * .5) +
     &            YN (LCORN (4)) +
     &            ( YN (LCORN (1)) * .5) ) * .25
               CALL ADD1CN (MXND, MLN, XN, YN, NUID, LXK, KXL, NXL, LXN,
     &            ANGLE, BNSIZE, LNODES, NNN, KKK, LLL, NNNOLD, LLLOLD,
     &            NLOOP, XNEW, YNEW, LCORN(1), IAVAIL, NAVAIL,
     &            GRAPH, VIDEO, SIZEIT, NOROOM, ERR, XNOLD, YNOLD,
     &            NXKOLD, LINKEG, LISTEG, BMESUR, MLINK, NPNOLD,
     &            NPEOLD, NNXK, REMESH, REXMIN, REXMAX, REYMIN, REYMAX,
     &            IDIVIS, SIZMIN, EMAX, EMIN)
               IF ((NOROOM) .OR. (ERR)) GOTO 110
            ELSEIF ( (I2+I3) .EQ. 4) THEN
               XNEW = ( XN (LNODES (3, LCORN(2))) +
     &            XN (LNODES (3, LCORN(3))) +
     &            ( XN (LCORN (4)) * .5) +
     &            XN (LCORN (1)) +
     &            ( XN (LCORN (2)) * .5) ) * .25
               YNEW = ( YN (LNODES (3, LCORN(2))) +
     &            YN (LNODES (3, LCORN(3))) +
     &            ( YN (LCORN (4)) * .5) +
     &            YN (LCORN (1)) +
     &            ( YN (LCORN (2)) * .5) ) * .25
               CALL ADD1CN (MXND, MLN, XN, YN, NUID, LXK, KXL, NXL, LXN,
     &            ANGLE, BNSIZE, LNODES, NNN, KKK, LLL, NNNOLD, LLLOLD,
     &            NLOOP, XNEW, YNEW, LCORN(2), IAVAIL, NAVAIL,
     &            GRAPH, VIDEO, SIZEIT, NOROOM, ERR, XNOLD, YNOLD,
     &            NXKOLD, LINKEG, LISTEG, BMESUR, MLINK, NPNOLD,
     &            NPEOLD, NNXK, REMESH, REXMIN, REXMAX, REYMIN, REYMAX,
     &            IDIVIS, SIZMIN, EMAX, EMIN)
               IF ((NOROOM) .OR. (ERR)) GOTO 110
            ELSEIF ( (I3+I4) .EQ. 4) THEN
               XNEW = ( XN (LNODES (3, LCORN(3))) +
     &            XN (LNODES (3, LCORN(4))) +
     &            ( XN (LCORN (1)) * .5) +
     &            XN (LCORN (2)) +
     &            ( XN (LCORN (3)) * .5) ) * .25
               YNEW = ( YN (LNODES (3, LCORN(3))) +
     &            YN (LNODES (3, LCORN(4))) +
     &            ( YN (LCORN (1)) * .5) +
     &            YN (LCORN (2)) +
     &            ( YN (LCORN (3)) * .5) ) * .25
               CALL ADD1CN (MXND, MLN, XN, YN, NUID, LXK, KXL, NXL, LXN,
     &            ANGLE, BNSIZE, LNODES, NNN, KKK, LLL, NNNOLD, LLLOLD,
     &            NLOOP, XNEW, YNEW, LCORN(3), IAVAIL, NAVAIL,
     &            GRAPH, VIDEO, SIZEIT, NOROOM, ERR, XNOLD, YNOLD,
     &            NXKOLD, LINKEG, LISTEG, BMESUR, MLINK, NPNOLD,
     &            NPEOLD, NNXK, REMESH, REXMIN, REXMAX, REYMIN, REYMAX,
     &            IDIVIS, SIZMIN, EMAX, EMIN)
               IF ((NOROOM) .OR. (ERR)) GOTO 110
            ELSE
               XNEW = ( XN (LNODES (3, LCORN(4))) +
     &            XN (LNODES (3, LCORN(1))) +
     &            ( XN (LCORN (2)) * .5) +
     &            XN (LCORN (3)) +
     &            ( XN (LCORN (4)) * .5) ) * .25
               YNEW = ( YN (LNODES (3, LCORN(4))) +
     &            YN (LNODES (3, LCORN(1))) +
     &            ( YN (LCORN (2)) * .5) +
     &            YN (LCORN (3)) +
     &            ( YN (LCORN (4)) * .5) ) * .25
               CALL ADD1CN (MXND, MLN, XN, YN, NUID, LXK, KXL, NXL, LXN,
     &            ANGLE, BNSIZE, LNODES, NNN, KKK, LLL, NNNOLD, LLLOLD,
     &            NLOOP, XNEW, YNEW, LCORN(4), IAVAIL, NAVAIL,
     &            GRAPH, VIDEO, SIZEIT, NOROOM, ERR, XNOLD, YNOLD,
     &            NXKOLD, LINKEG, LISTEG, BMESUR, MLINK, NPNOLD,
     &            NPEOLD, NNXK, REMESH, REXMIN, REXMAX, REYMIN, REYMAX,
     &            IDIVIS, SIZMIN, EMAX, EMIN)
               IF ((NOROOM) .OR. (ERR)) GOTO 110
            ENDIF
C
C  HANDLE A 2-1-2-1 RECTANGLE
C
         ELSE
            IF (I1 .EQ. 2) THEN
               CALL ADD2EL (MXND, MLN, XN, YN, NUID, LXK, KXL, NXL, LXN,
     &            ANGLE, LNODES, NNN, KKK, LLL, NNNOLD, LLLOLD,
     &            NLOOP, LNODES (3, LCORN(1)),
     &            LNODES (3, LCORN(3)), IAVAIL, NAVAIL,
     &            GRAPH, VIDEO, NOROOM, ERR)
               IF ((NOROOM) .OR. (ERR)) GOTO 110
            ELSE
               CALL ADD2EL (MXND, MLN, XN, YN, NUID, LXK, KXL, NXL, LXN,
     &            ANGLE, LNODES, NNN, KKK, LLL, NNNOLD, LLLOLD,
     &            NLOOP, LNODES (3, LCORN(2)),
     &            LNODES (3, LCORN(4)), IAVAIL, NAVAIL,
     &            GRAPH, VIDEO, NOROOM, ERR)
               IF ((NOROOM) .OR. (ERR)) GOTO 110
            ENDIF
         ENDIF
C
C  PROCESS A SEMICIRCLE SHAPE WITH VARIOUS INTERVAL COMBINATIONS
C
      ELSEIF (NCORN .EQ. 2) THEN
         I1 = LNODES (7, LCORN(1))
         I2 = LNODES (7, LCORN(2))
C
C  HANDLE A 5-1 SEMICIRCLE
C
         IF (MAX0 (I1, I2) .EQ. 5) THEN
            IF (I1 .EQ. 1) THEN
               CALL ADD2EL (MXND, MLN, XN, YN, NUID, LXK, KXL, NXL, LXN,
     &            ANGLE, LNODES, NNN, KKK, LLL, NNNOLD, LLLOLD,
     &            NLOOP, LNODES (2, LCORN(1)),
     &            LNODES (3, LCORN(2)), IAVAIL, NAVAIL,
     &            GRAPH, VIDEO, NOROOM, ERR)
               IF ((NOROOM) .OR. (ERR)) GOTO 110
            ELSE
               CALL ADD2EL (MXND, MLN, XN, YN, NUID, LXK, KXL, NXL, LXN,
     &            ANGLE, LNODES, NNN, KKK, LLL, NNNOLD, LLLOLD,
     &            NLOOP, LNODES (2, LCORN(2)),
     &            LNODES (3, LCORN(1)), IAVAIL, NAVAIL,
     &            GRAPH, VIDEO, NOROOM, ERR)
               IF ((NOROOM) .OR. (ERR)) GOTO 110
            ENDIF
C
C  HANDLE A 4-2 SEMICIRCLE
C
         ELSEIF (MAX0 (I1, I2) .EQ. 4) THEN
            IF (I1 .EQ. 2) THEN
               CALL ADD2EL (MXND, MLN, XN, YN, NUID, LXK, KXL, NXL, LXN,
     &            ANGLE, LNODES, NNN, KKK, LLL, NNNOLD, LLLOLD,
     &            NLOOP, LNODES (3, LCORN(1)),
     &            LNODES (3, LNODES (3, LCORN(2))), IAVAIL, NAVAIL,
     &            GRAPH, VIDEO, NOROOM, ERR)
               IF ((NOROOM) .OR. (ERR)) GOTO 110
            ELSE
               CALL ADD2EL (MXND, MLN, XN, YN, NUID, LXK, KXL, NXL, LXN,
     &            ANGLE, LNODES, NNN, KKK, LLL, NNNOLD, LLLOLD,
     &            NLOOP, LNODES (3, LCORN(2)),
     &            LNODES (3, LNODES (3, LCORN(1))), IAVAIL, NAVAIL,
     &            GRAPH, VIDEO, NOROOM, ERR)
               IF ((NOROOM) .OR. (ERR)) GOTO 110
            ENDIF
C
C  HANDLE A 3-3 SEMICIRCLE
C
         ELSE
            CALL ADD2EL (MXND, MLN, XN, YN, NUID, LXK, KXL, NXL, LXN,
     &         ANGLE, LNODES, NNN, KKK, LLL, NNNOLD, LLLOLD,
     &         NLOOP, LNODES (3, LCORN(1)),
     &         LNODES (3, LCORN(2)), IAVAIL, NAVAIL,
     &         GRAPH, VIDEO, NOROOM, ERR)
            IF ((NOROOM) .OR. (ERR)) GOTO 110
C            CALL ADD2CN (MXND, MLN, XN, YN, NUID, LXK, KXL, NXL, LXN,
C     &         ANGLE, BNSIZE, LNODES, NNN, KKK, LLL, NNNOLD, LLLOLD,
C     &         NLOOP, LCORN(1), IAVAIL, NAVAIL, GRAPH, VIDEO, SIZEIT,
C     &            NOROOM, ERR, XNOLD, YNOLD, NXKOLD, LINKEG, LISTEG,
C     &            BMESUR, MLINK, NPNOLD, NPEOLD, NNXK, REMESH, REXMIN,
C     &            REXMAX, REYMIN, REYMAX, IDIVIS, SIZMIN, EMAX, EMIN)
C            IF ((NOROOM) .OR. (ERR)) GOTO 110
         ENDIF
C
C  PROCESS A TEAR-DROP SHAPE (1 CORNER)
C
      ELSEIF (NCORN .EQ. 1) THEN
         CALL ADD2EL (MXND, MLN, XN, YN, NUID, LXK, KXL, NXL, LXN,
     &      ANGLE, LNODES, NNN, KKK, LLL, NNNOLD, LLLOLD,
     &      NLOOP, LNODES (3, LCORN(1)),
     &      LNODES (2, LNODES (2, LCORN(1))), IAVAIL, NAVAIL,
     &      GRAPH, VIDEO, NOROOM, ERR)
         IF ((NOROOM) .OR. (ERR)) GOTO 110
C
C  PROCESS A TRUE CIRCLE (OR ANYTHING ELSE FOR THAT MATTER)
C
      ELSE
         I1 = NODE
         I2 = LNODES (3, NODE)
         I3 = LNODES (3, LNODES (3, NODE))
         I4 = LNODES (3, LNODES (3, LNODES (3, NODE)))
         I5 = LNODES (2, LNODES (2, NODE))
         I6 = LNODES (2, NODE)
C
         DIST1 = SQRT ( ((XN (I1) - XN (I4)) ** 2) +
     &      ((YN (I1) - YN (I4)) ** 2) )
         DIST2 = SQRT ( ((XN (I2) - XN (I5)) ** 2) +
     &      ((YN (I2) - YN (I5)) ** 2) )
         DIST3 = SQRT ( ((XN (I3) - XN (I6)) ** 2) +
     &      ((YN (I3) - YN (I6)) ** 2) )
C
         IF ( (DIST1 .LE. DIST2) .AND. (DIST1 .LE. DIST3) ) THEN
            CALL ADD2EL (MXND, MLN, XN, YN, NUID, LXK, KXL, NXL, LXN,
     &         ANGLE, LNODES, NNN, KKK, LLL, NNNOLD, LLLOLD,
     &         NLOOP, I1, I4, IAVAIL, NAVAIL, GRAPH, VIDEO,
     &         NOROOM, ERR)
            IF ((NOROOM) .OR. (ERR)) GOTO 110
         ELSEIF ( (DIST2 .LE. DIST3) .AND. (DIST2 .LE. DIST1) ) THEN
            CALL ADD2EL (MXND, MLN, XN, YN, NUID, LXK, KXL, NXL, LXN,
     &         ANGLE, LNODES, NNN, KKK, LLL, NNNOLD, LLLOLD,
     &         NLOOP, I2, I5, IAVAIL, NAVAIL, GRAPH, VIDEO,
     &         NOROOM, ERR)
            IF ((NOROOM) .OR. (ERR)) GOTO 110
         ELSE
            CALL ADD2EL (MXND, MLN, XN, YN, NUID, LXK, KXL, NXL, LXN,
     &         ANGLE, LNODES, NNN, KKK, LLL, NNNOLD, LLLOLD,
     &         NLOOP, I3, I6, IAVAIL, NAVAIL, GRAPH, VIDEO,
     &         NOROOM, ERR)
            IF ((NOROOM) .OR. (ERR)) GOTO 110
         ENDIF
      ENDIF
C
  110 CONTINUE
C
      RETURN
C
      END
