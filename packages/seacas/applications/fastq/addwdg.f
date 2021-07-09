C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE ADDWDG (MXND, MLN, NUID, XN, YN, ZN, LXK, KXL, NXL,
     &   LXN, LNODES, ANGLE, BNSIZE, NLOOP, IAVAIL, NAVAIL, LLL, KKK,
     &   NNN, LLLOLD, NNNOLD, TANG, KANG, NSTART, NEND, XMIN, XMAX,
     &   YMIN, YMAX, ZMIN, ZMAX, GRAPH, VIDEO, DEV1, KREG, NOROOM, ERR)
C***********************************************************************

C  SUBROUTINE ADDWDG = ADDS WEDGES IN A ROW

C***********************************************************************

C  ADD WEDGES BASED ON THE TOTAL TURNED ANGLE:
C      FOR TURNING ANGLES LESS THAN 135 DEGREES - 1 WEDGE
C      FOR TURNING ANGLES BETWEEN 135 AND 225 DEGREES - TRY 2 WEDGES
C      FOR TURNING ANGLES BETWEEN 225 AND 315 DEGREES - TRY 3 WEDGES
C      FOR TURNING ANGLES GREATER THAN 315 DEGREES - TRY 4 WEDGES

C***********************************************************************

      DIMENSION XN (MXND), YN (MXND), ZN (MXND), NUID (MXND)
      DIMENSION LXK (4, MXND), KXL (2, 3*MXND)
      DIMENSION NXL (2, 3*MXND), LXN (4, MXND)
      DIMENSION LNODES (MLN, MXND), BNSIZE (2, MXND), ANGLE (MXND)
      DIMENSION INODE (4)

      LOGICAL GRAPH, ERR, MAXSIZ, VIDEO, NOROOM, PWEDGE

      CHARACTER*3 DEV1

      MAXSIZ = .TRUE.
      ERR = .FALSE.
      PWEDGE = .FALSE.

      IF (TANG .LT. 2.3561945) THEN
         NWANT = 1
      ELSEIF (TANG .LT. 3.9269908) THEN
         NWANT = 2
      ELSEIF (TANG .LT. 5.4977871) THEN
         NWANT = 3
      ELSE
         NWANT = 4
      ENDIF

      CALL NSPLIT (MXND, MLN, LNODES, ANGLE, NSTART, KANG, INODE,
     &   NNODE, NWANT, MAXSIZ)

      DO 100 I = 1, NNODE
         CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES,
     &      LNODES (2, INODE(I)), ERR)
         IF (ERR) GOTO 110
         CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES,
     &      LNODES (2, LNODES (2, INODE(I))), ERR)
         IF (ERR) GOTO 110
         CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES,
     &      LNODES (3, INODE(I)), ERR)
         IF (ERR) GOTO 110
         CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES,
     &      LNODES (3, LNODES (3, INODE(I))), ERR)
         IF (ERR) GOTO 110
         CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES,
     &      INODE(I), ERR)
         IF (ERR) GOTO 110
         CALL WEDGE (MXND, MLN, NUID, LXK, KXL, NXL, LXN, XN, YN,
     &      LNODES, BNSIZE, IAVAIL, NAVAIL, LLL, KKK, NNN, LLLOLD,
     &      NNNOLD, INODE(I), IDUM, NLOOP, PWEDGE, GRAPH, VIDEO,
     &      NOROOM, ERR)
         IF ((NOROOM) .OR. (ERR)) GOTO 110
         IF (VIDEO) THEN
            CALL RPLOTL (MXND, XN, YN, ZN, NXL, XMIN, XMAX, YMIN,
     &         YMAX, ZMIN, ZMAX, LLL, DEV1, KREG)
            CALL SNAPIT (2)
         ENDIF
  100 CONTINUE

  110 CONTINUE
      RETURN

      END
