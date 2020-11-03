C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE CLOSE6 (MXND, MXCORN, MLN, NUID, XN, YN, LXK, KXL, NXL,
     &   LXN, ANGLE, BNSIZE, LNODES, NODE, NLOOP, KKKOLD, LLLOLD,
     &   NNNOLD, NAVAIL, IAVAIL, DONE, XMIN, XMAX, YMIN, YMAX, DEV1,
     &   LLL, KKK, NNN, LCORN, NCORN, GRAPH, VIDEO, SIZEIT, NOROOM,
     &   ERR, XNOLD, YNOLD, NXKOLD, LINKEG, LISTEG, BMESUR, MLINK,
     &   NPNOLD, NPEOLD, NNXK, REMESH, REXMIN, REXMAX, REYMIN, REYMAX,
     &   IDIVIS, SIZMIN, EMAX, EMIN)
C***********************************************************************

C  SUBROUTINE CLOSE6 = FINISHES UP A LOOP WITH ONLY 6 LINES IN THE LOOP

C***********************************************************************

      DIMENSION XN (MXND), YN (MXND), NUID (MXND)
      DIMENSION LXK (4, MXND), KXL (2, 3*MXND)
      DIMENSION NXL (2, 3*MXND), LXN (4, MXND)
      DIMENSION ANGLE (MXND), LNODES (MLN, MXND), BNSIZE (2, MXND)
      DIMENSION LCORN (MXCORN)

      DIMENSION XNOLD(NPNOLD), YNOLD(NPNOLD)
      DIMENSION NXKOLD(NNXK, NPEOLD)
      DIMENSION LINKEG(2, MLINK), LISTEG(4 * NPEOLD), BMESUR(NPNOLD)

      LOGICAL GRAPH, VIDEO, ERR, SIZEIT, NOROOM

      CHARACTER*3 DEV1

      ERR = .FALSE.

C  SET ALL THE LOOP NODES TO BE INTERIOR

      INODE = NODE
      IKOUNT = 0
  100 CONTINUE
      IKOUNT = IKOUNT + 1
      IF (IKOUNT .GT. 6) THEN
         CALL MESAGE ('** PROBLEMS IN CLOSE6 WITH TOO MANY IN LOOP **')
         ERR = .TRUE.
         GOTO 110
      ENDIF

      LNODES (4, INODE) = - 2
      CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES, INODE, ERR)
      IF (ERR) GOTO 110

      INODE = LNODES (3, INODE)
      IF (INODE .NE. NODE) GOTO 100

C  NOW GET THE CORRECT INTERPRETATION OF THE SHAPE

      CALL CNTCRN (MXND, MXCORN, MLN, LNODES, LCORN, NCORN, NLOOP,
     &   NODE, ERR)
      IF (ERR) GOTO 110

C  PROCESS A TRIANGLE SHAPE WITH VARIOUS INTERVAL COMBINATIONS

      IF (NCORN .EQ. 3) THEN
         I1 = LNODES (7, LCORN(1))
         I2 = LNODES (7, LCORN(2))
         I3 = LNODES (7, LCORN(3))

C  HANDLE A 4-1-1 TRIANGLE

         IF (MAX0 (I1, I2, I3) .EQ. 4) THEN
            IF (I1 .EQ. 4) THEN
               XNEW = ( ( (XN (LCORN (1)) + XN (LCORN (2)) ) * .5) +
     &            XN (LCORN (3)) ) * .5
               YNEW = ( ( (YN (LCORN (1)) + YN (LCORN (2)) ) * .5) +
     &            YN (LCORN (3)) ) * .5

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

C  HANDLE A 3-2-1 TRIANGLE

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

C  HANDLE A 2-2-2 TRIANGLE

         ELSE
            XNEW1 = ( ( (XN (LCORN (1)) + XN (LCORN (2)) ) * .5) +
     &         XN (LCORN (3)) ) * .5
            XNEW2 = ( ( (XN (LCORN (2)) + XN (LCORN (3)) ) * .5) +
     &         XN (LCORN (1)) ) * .5
            XNEW3 = ( ( (XN (LCORN (3)) + XN (LCORN (1)) ) * .5) +
     &         XN (LCORN (2)) ) * .5

            YNEW1 = ( ( (YN (LCORN (1)) + YN (LCORN (2)) ) * .5) +
     &         YN (LCORN (3)) ) * .5
            YNEW2 = ( ( (YN (LCORN (2)) + YN (LCORN (3)) ) * .5) +
     &         YN (LCORN (1)) ) * .5
            YNEW3 = ( ( (YN (LCORN (3)) + YN (LCORN (1)) ) * .5) +
     &         YN (LCORN (2)) ) * .5

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

C  PROCESS A RECTANGLE SHAPE WITH VARIOUS INTERVAL COMBINATIONS

      ELSEIF (NCORN .EQ. 4) THEN
         I1 = LNODES (7, LCORN(1))
         I2 = LNODES (7, LCORN(2))
         I3 = LNODES (7, LCORN(3))
         I4 = LNODES (7, LCORN(4))

C  HANDLE A 3-1-1-1 RECTANGLE

         IF (MAX0 (I1, I2, I3, I4) .EQ. 3) THEN
            IF (I1 .EQ. 3) THEN

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

C  HANDLE A 2-2-1-1 RECTANGLE

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

C  HANDLE A 2-1-2-1 RECTANGLE

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

C  PROCESS A SEMICIRCLE SHAPE WITH VARIOUS INTERVAL COMBINATIONS

      ELSEIF (NCORN .EQ. 2) THEN
         I1 = LNODES (7, LCORN(1))
         I2 = LNODES (7, LCORN(2))

C  HANDLE A 5-1 SEMICIRCLE

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

C  HANDLE A 4-2 SEMICIRCLE

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

C  HANDLE A 3-3 SEMICIRCLE

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

C  PROCESS A TEAR-DROP SHAPE (1 CORNER)

      ELSEIF (NCORN .EQ. 1) THEN
         CALL ADD2EL (MXND, MLN, XN, YN, NUID, LXK, KXL, NXL, LXN,
     &      ANGLE, LNODES, NNN, KKK, LLL, NNNOLD, LLLOLD,
     &      NLOOP, LNODES (3, LCORN(1)),
     &      LNODES (2, LNODES (2, LCORN(1))), IAVAIL, NAVAIL,
     &      GRAPH, VIDEO, NOROOM, ERR)
         IF ((NOROOM) .OR. (ERR)) GOTO 110

C  PROCESS A TRUE CIRCLE (OR ANYTHING ELSE FOR THAT MATTER)

      ELSE
         I1 = NODE
         I2 = LNODES (3, NODE)
         I3 = LNODES (3, LNODES (3, NODE))
         I4 = LNODES (3, LNODES (3, LNODES (3, NODE)))
         I5 = LNODES (2, LNODES (2, NODE))
         I6 = LNODES (2, NODE)

         DIST1 = SQRT ( ((XN (I1) - XN (I4)) ** 2) +
     &      ((YN (I1) - YN (I4)) ** 2) )
         DIST2 = SQRT ( ((XN (I2) - XN (I5)) ** 2) +
     &      ((YN (I2) - YN (I5)) ** 2) )
         DIST3 = SQRT ( ((XN (I3) - XN (I6)) ** 2) +
     &      ((YN (I3) - YN (I6)) ** 2) )

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

  110 CONTINUE

      RETURN

      END
