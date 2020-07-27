C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE ADJTRI (MXND, MLN, LNODES, XN, YN, ZN, NUID, LXK, KXL,
     &   NXL, LXN, NNN, NAVAIL, IAVAIL, NODE, KELEM, ANG, TOLER1,
     &   TOLER2, N1, N2, N3, KREG, XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX,
     &   KKK, LLL, DEV1, DONE, CHECK, GRAPH, VIDEO, NOROOM, ERR, KKKADD)
C***********************************************************************

C  SUBROUTINE ADJTRI = ADJUSTS A TRIANGULAR SHAPED ELEMENT WHERE
C                      POSSIBLE

C***********************************************************************

C  SUBROUTINE CALLED BY TRIDEL

C***********************************************************************

C  THERE ARE THREE POSSIBILITIES FOR CHANGE:
C     1) ANYTHING OVER TOLER1 GETS THE CORRESPONDING ELEMENT
C        DELETED
C     2) ANYTHING OVER TOLER2 AND HOOKED TO ANOTHER 3-LINE NODE GETS
C        THE CORRESPONDING ELEMENT DELETED
C     3) AN ELONGATED ELEMENT OVER 150 DEGREES GETS A 3 ELEMENT
C        REPLACEMENT FOR THE TWO ELEMENTS THERE

C***********************************************************************

      DIMENSION LXK(4, MXND), NXL(2, 3*MXND), KXL(2, 3*MXND)
      DIMENSION LXN(4, MXND), XN(MXND), YN(MXND), ZN(MXND), NUID(MXND)
      DIMENSION LNODES (MLN, MXND)

      CHARACTER*3 DEV1

      LOGICAL NOROOM, ERR, DONE, GRAPH, CHECK, VIDEO

C  CHECK FOR CASE 1

      IF (ANG .GT. TOLER1) THEN
         IF (GRAPH) THEN
            CALL RPLOTL (MXND, XN, YN, ZN, NXL, XMIN, XMAX, YMIN,
     &         YMAX, ZMIN, ZMAX, LLL, DEV1, KREG)
            CALL LCOLOR ('PINK ')
            CALL D2NODE (MXND, XN, YN, NODE, N1)
            CALL D2NODE (MXND, XN, YN, NODE, N2)
            CALL LCOLOR ('WHITE')
            CALL SFLUSH
         ENDIF
  100    CONTINUE
         CALL DELEM (MXND, XN, YN, NUID, LXK, KXL, NXL, LXN,
     &      NNN, NAVAIL, IAVAIL, NODE, KELEM, NODE1, NODE3, DONE,
     &      CHECK, NOROOM, ERR)
         IF ((NOROOM) .OR. (ERR)) GOTO 120
         IF (LXN (3, NODE1) .LE. 0) THEN
            NODE = NODE1
            KELEM = KXL (1, LXN (1, NODE))
            CHECK = .FALSE.
            GOTO 100
         ELSEIF (LXN (3, NODE3) .LE. 0) THEN
            NODE = NODE3
            KELEM = KXL (1, LXN (1, NODE))
            CHECK = .FALSE.
            GOTO 100
         ENDIF
         CHECK = .TRUE.
         IF ((ERR) .OR. (DONE)) GOTO 120
      ENDIF

C  CHECK FOR CASE 2

      IF ( (ANG .GT. TOLER2) .AND. (LXN (4, N3) .EQ. 0) .AND.
     &   (LXN (3, N3) .GT. 0) .AND. (LXN (2, N3) .GT. 0)) THEN
         IF (GRAPH) THEN
            CALL RPLOTL (MXND, XN, YN, ZN, NXL, XMIN, XMAX, YMIN,
     &         YMAX, ZMIN, ZMAX, LLL, DEV1, KREG)
            CALL LCOLOR ('PINK ')
            CALL D2NODE (MXND, XN, YN, NODE, N1)
            CALL D2NODE (MXND, XN, YN, NODE, N2)
            CALL LCOLOR ('WHITE')
            CALL SFLUSH
         ENDIF
  110    CONTINUE
         CALL DELEM (MXND, XN, YN, NUID, LXK, KXL, NXL, LXN,
     &      NNN, NAVAIL, IAVAIL, NODE, KELEM, NODE1, NODE3, DONE,
     &      CHECK, NOROOM, ERR)
         IF ((NOROOM) .OR. (ERR)) GOTO 120
         IF (LXN (3, NODE1) .LE. 0) THEN
            NODE = NODE1
            KELEM = KXL (1, LXN (1, NODE))
            CHECK = .FALSE.
            GOTO 110
         ELSEIF (LXN (3, NODE3) .LE. 0) THEN
            NODE = NODE3
            KELEM = KXL (1, LXN (1, NODE))
            CHECK = .FALSE.
            GOTO 110
         ENDIF
         CHECK = .TRUE.
         IF ((ERR) .OR. (DONE)) GOTO 120
      ENDIF

C  CHECK FOR CASE 3

      CALL LONGEL (MXND, MLN, LNODES, XN, YN, NUID, LXK, KXL, NXL,
     &   LXN, NNN, NAVAIL, IAVAIL, NODE, KELEM, ANG, TOLER2,
     &   N1, N2, KREG, XMIN, XMAX, YMIN, YMAX, KKK, LLL, DONE, GRAPH,
     &   VIDEO, NOROOM, ERR, KKKADD)

  120 CONTINUE
      RETURN

      END
