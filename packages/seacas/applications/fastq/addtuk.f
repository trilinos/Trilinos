C $Id: addtuk.f,v 1.1 1990/11/30 11:03:10 gdsjaar Exp $
C $Log: addtuk.f,v $
C Revision 1.1  1990/11/30 11:03:10  gdsjaar
C Initial revision
C
C
CC* FILE: [.PAVING]ADDTUK.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE ADDTUK (MXND, MLN, NUID, XN, YN, ZN, LXK, KXL, NXL,
     &   LXN, LNODES, ANGLE, NLOOP, IAVAIL, NAVAIL, LLL, KKK, NNN, TANG,
     &   KANG, NSTART, NEND, NODE, XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX,
     &   GRAPH, VIDEO, DEV1, NOROOM, ERR)
C***********************************************************************
C
C  SUBROUTINE ADDTUK = ADDS TUCKS IN A ROW
C
C***********************************************************************
C
C  ADD TUCKS BASED ON THE TOTAL TURNED ANGLE:
C      FOR TURNING ANGLES LESS THAN 135 DEGREES - 1 TUCK
C      FOR TURNING ANGLES BETWEEN 135 AND 225 DEGREES - TRY 2 TUCKS
C      FOR TURNING ANGLES BETWEEN 225 AND 315 DEGREES - TRY 3 TUCKS
C      FOR TURNING ANGLES GREATER THAN 315 DEGREES - TRY 4 TUCKS
C
C***********************************************************************
C
      DIMENSION XN (MXND), YN (MXND), ZN (MXND), NUID (MXND)
      DIMENSION LXK (4, MXND), KXL (2, 3*MXND)
      DIMENSION NXL (2, 3*MXND), LXN (4, MXND)
      DIMENSION LNODES (MLN, MXND), ANGLE (MXND)
      DIMENSION INODE (4)
C
      LOGICAL GRAPH, ERR, MAXSIZ, VIDEO, NOROOM
C
      CHARACTER*3 DEV1
C
      ERR = .FALSE.
      MAXSIZ = .FALSE.
C
      IF (TANG .LT. 2.3561945) THEN
         NWANT = 1
      ELSEIF (TANG .LT. 3.9269908) THEN
         IF (KANG .GT. 2) THEN
            NWANT = 2
         ELSE
            NWANT = 1
         ENDIF
      ELSEIF (TANG .LT. 5.4977871) THEN
         IF (KANG .GT. 4) THEN
            NWANT = 3
         ELSEIF (KANG .GT. 2) THEN
            NWANT = 2
         ELSE
            NWANT = 1
         ENDIF
      ELSE
         IF (KANG .GT. 6) THEN
            NWANT = 4
         ELSEIF (KANG .GT. 4) THEN
            NWANT = 3
         ELSEIF (KANG .GT. 2) THEN
            NWANT = 2
         ELSE
            NWANT = 1
         ENDIF
      ENDIF
C
      CALL NSPLIT (MXND, MLN, LNODES, ANGLE, NSTART, KANG, INODE,
     &   NNODE, NWANT, MAXSIZ)
C
      DO 100 I = 1, NNODE
         IF (LXN (1, I) .GT. 0) THEN
C
C  MARK THE SMOOTHING
C
            CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES,
     &         LNODES (2, INODE(I)), ERR)
            IF (ERR) GOTO 110
            CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES,
     &         LNODES (2, LNODES (2, INODE(I))), ERR)
            IF (ERR) GOTO 110
            CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES,
     &         LNODES (3, INODE(I)), ERR)
            IF (ERR) GOTO 110
            CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES,
     &         LNODES (3, LNODES (3, INODE(I))), ERR)
            IF (ERR) GOTO 110
            CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES,
     &         LNODES (3, LNODES (3, LNODES (3, INODE(I)))), ERR)
            IF (ERR) GOTO 110
C
C  MAKE SURE THAT THE NSTART, NEND, AND NODE GET UPDATED IF THEY ARE TO
C  BE DELETED
C
            IF ( (INODE(I) .EQ. NSTART) .OR.
     &         (LNODES (3, INODE(I)) .EQ. NSTART)) THEN
               NSTART = LNODES (3, LNODES (3, INODE(I)))
            ENDIF
            IF ( (INODE(I) .EQ. NEND) .OR.
     &         (LNODES (3, INODE(I)) .EQ. NEND) ) THEN
               NEND = LNODES (3, LNODES (3, INODE(I)))
            ENDIF
            IF ( (INODE(I) .EQ. NODE) .OR.
     &         (LNODES (3, INODE(I)) .EQ. NODE) ) THEN
               NODE = LNODES (3, LNODES (3, INODE(I)))
            ENDIF
C
C  TAKE THE TUCK
C
            CALL TUCK (MXND, MLN, NUID, XN, YN, LXK, KXL, NXL, LXN,
     &         LNODES, IAVAIL, NAVAIL, LLL, KKK, NNN, INODE (I), NLOOP,
     &         GRAPH, NOROOM, ERR)
            IF ((NOROOM) .OR. (ERR)) GOTO 110
            IF (VIDEO) THEN
               CALL RPLOTL (MXND, XN, YN, ZN, NXL, XMIN, XMAX, YMIN,
     &            YMAX, ZMIN, ZMAX, LLL, DEV1, KREG)
               CALL SNAPIT (1)
            ENDIF
         ENDIF
  100 CONTINUE
C
  110 CONTINUE
C
      RETURN
C
      END
