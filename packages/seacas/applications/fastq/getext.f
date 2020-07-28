C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE GETEXT (MP, ML, MS, MR, N, IPOINT, COOR, ILINE, LTYPE,
     &   LCON, NLPS, IFLINE, ILLIST, NSPR, IFSIDE, ISLIST, LINKP,
     &   LINKL, LINKS, LINKR, REXTRM, BXMIN, BXMAX, BYMIN, BYMAX)
C***********************************************************************

C  SUBROUTINE GETEXT = GETS THE REGION AND BODY EXTREMES

C***********************************************************************

      DIMENSION IPOINT (MP), COOR (2, MP)
      DIMENSION ILINE (ML), LTYPE (ML), LCON (3, ML)
      DIMENSION NLPS (MS), IFLINE (MS), ILLIST (MS * 3)
      DIMENSION NSPR (MR), IFSIDE (MR), ISLIST (MR * 4)
      DIMENSION LINKP (2, MP), LINKL (2, ML), LINKS (2, MS)
      DIMENSION LINKR (2, MR)
      DIMENSION REXTRM (4, MR), N (29)

      LOGICAL FOUND, GETMAX, ADDLNK
      LOGICAL NUMPLT, TEST

C  GET THE POINTS EXTREMES

      ADDLNK = .FALSE.
      GETMAX = .TRUE.
      FOUND = .FALSE.
      DO 100 I = 1, N (1)
         CALL LTSORT (MP, LINKP, IABS (IPOINT (I)), II, ADDLNK)
         IF (II .GT. 0) THEN
            IF (FOUND) THEN
               BXMAX = AMAX1 (COOR (1, II), BXMAX)
               BYMAX = AMAX1 (COOR (2, II), BYMAX)
               BXMIN = AMIN1 (COOR (1, II), BXMIN)
               BYMIN = AMIN1 (COOR (2, II), BYMIN)
            ELSE
               BXMAX = COOR (1, II)
               BXMIN = COOR (1, II)
               BYMAX = COOR (2, II)
               BYMIN = COOR (2, II)
               FOUND = .TRUE.
            ENDIF
         ENDIF
  100 CONTINUE

C  GET ALL THE LINES EXTREMES

      IF (FOUND) THEN
         DO 110 I = 1, N (2)
            CALL LTSORT (ML, LINKL, IABS (ILINE (I)), II, ADDLNK)
            IF (II .GT. 0) THEN
               CALL DLINE (MP, ML, COOR, LINKP, ILINE (II),
     &            LTYPE (II), LCON (1, II), LCON (2, II), LCON (3, II),
     &            NUMPLT, X1, Y1, TEST, GETMAX, BXMIN, BXMAX, BYMIN,
     &            BYMAX)
            ENDIF
  110    CONTINUE
      ELSE
         BXMIN = 0.
         BXMAX = 20.
         BYMIN = 0.
         BYMAX = 15.
      ENDIF

C  CALCULATE THE EXTREMES FOR EACH REGION

      DO 120 I = 1, N (22)
         CALL LTSORT (MR, LINKR, I, II, ADDLNK)
         IF (II .GT. 0) THEN
            CALL REGEXT (MP, ML, MS, MR, N, II, COOR, ILINE, LTYPE,
     &         LCON, NLPS, IFLINE, ILLIST, NSPR, IFSIDE, ISLIST, LINKP,
     &         LINKL, LINKS, LINKR, XMIN, XMAX, YMIN, YMAX)
            REXTRM (1, II) = XMIN
            REXTRM (2, II) = XMAX
            REXTRM (3, II) = YMIN
            REXTRM (4, II) = YMAX
         ENDIF
  120 CONTINUE

      RETURN

      END
