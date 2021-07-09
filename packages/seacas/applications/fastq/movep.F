C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE MOVEP (MP, ML, MS, MR, N, COOR, LTYPE, LCON, NLPS,
     &   IFLINE, ILLIST, NSPR, IFSIDE, ISLIST, LINKP, LINKL, LINKS,
     &   LINKR, IMOVE, ALPHA, X, Y)
C***********************************************************************

C  SUBROUTINE MOVEP = MOVES A POINT, AND REDRAWS ALL LINES DEPENDENT ON
C                     THAT POINT

C***********************************************************************

      DIMENSION COOR (2, MP), LTYPE (ML), LCON (3, ML)
      DIMENSION NLPS (MS), IFLINE (MS)
      DIMENSION ILLIST (MS * 3)
      DIMENSION NSPR (MR), IFSIDE (MR), ISLIST (MR * 4)
      DIMENSION LINKP (2, MP), LINKL (2, ML), LINKS (2, MS)
      DIMENSION LINKR (2, MR)
      DIMENSION N (29)

      LOGICAL NUMPLT, ALPHA, ADDLNK, TEST, GETMAX

      KNUM = 0
      ADDLNK = .FALSE.
      TEST = .FALSE.
      GETMAX = .FALSE.
      IF (.NOT.ALPHA)CALL PLTSTD (1, 0.)

C  ERASE ANY LINES USING THAT POINT

      DO 120 I = 1, N (19)
         CALL LTSORT (ML, LINKL, I, II, ADDLNK)
         IF (II .GT. 0) THEN
            DO 100 J = 1, 3
               IF (IMOVE .EQ. LCON (J, II)) THEN
                  LT = LTYPE (II)
                  IP1 = LCON (1, II)
                  IP2 = LCON (2, II)
                  IP3 = LCON (3, II)
                  CALL LTSORT (MP, LINKP, IP1, IPNTR1, ADDLNK)
                  CALL LTSORT (MP, LINKP, IP2, IPNTR2, ADDLNK)
                  IF (IP3 .NE. 0) THEN
                     CALL LTSORT (MP, LINKP, IABS (IP3), IPNTR3, ADDLNK)
                  ELSE
                     IPNTR3 = 0
                  ENDIF
                  IF ( (IPNTR1 .GT. 0) .AND. (IPNTR2 .GT. 0) .AND.
     &               ( (LT .EQ. 1) .OR. (IPNTR3 .GT. 0))) THEN
                     IF (.NOT.ALPHA) THEN
                        CALL DLINE (MP, ML, COOR, LINKP, KNUM, LT, IP1,
     &                     IP2, IP3, NUMPLT, DUM1, DUM2, TEST, GETMAX,
     &                     DUM3, DUM4, DUM5, DUM6)
                        CALL PLTFLU
                     ENDIF
                  ENDIF
                  GOTO 110
               ENDIF
  100       CONTINUE
  110       CONTINUE
         ENDIF
  120 CONTINUE

C  ERASE THE POINT

      ADDLNK = .FALSE.
      CALL LTSORT (MP, LINKP, IMOVE, II, ADDLNK)
      IF (.NOT.ALPHA) THEN
#if NeedsDoubleEscape
         CALL MPD2SY (1, COOR (1, II), COOR (2, II), '\\CX')
         CALL MPD2SY (1, COOR (1, II), COOR (2, II), '\\CSQ')
#else
         CALL MPD2SY (1, COOR (1, II), COOR (2, II), '\CX')
         CALL MPD2SY (1, COOR (1, II), COOR (2, II), '\CSQ')
#endif
      ENDIF

C  REDEFINE THE POINT

      COOR (1, II) = X
      COOR (2, II) = Y

      IF (.NOT.ALPHA)CALL PLTSTD (1, 7.)

C  REDRAW ANY LINES USING THAT POINT

      DO 150 I = 1, N (19)
         CALL LTSORT (ML, LINKL, I, II, ADDLNK)
         IF (II .GT. 0) THEN
            DO 130 J = 1, 3
               IF (IMOVE .EQ. LCON (J, II)) THEN
                  LT = LTYPE (II)
                  IP1 = LCON (1, II)
                  IP2 = LCON (2, II)
                  IP3 = LCON (3, II)
                  CALL LTSORT (MP, LINKP, IP1, IPNTR1, ADDLNK)
                  CALL LTSORT (MP, LINKP, IP2, IPNTR2, ADDLNK)
                  IF (IP3 .NE. 0) THEN
                     CALL LTSORT (MP, LINKP, IABS (IP3), IPNTR3, ADDLNK)
                  ELSE
                     IPNTR3 = 0
                  ENDIF
                  IF ( (IPNTR1 .GT. 0) .AND. (IPNTR2 .GT. 0) .AND.
     &               ( (LT .EQ. 1) .OR. (IPNTR3 .GT. 0))) THEN
                     IF (.NOT.ALPHA) THEN
                        CALL DLINE (MP, ML, COOR, LINKP, KNUM, LT, IP1,
     &                     IP2, IP3, NUMPLT, DUM1, DUM2, TEST, GETMAX,
     &                     DUM3, DUM4, DUM5, DUM6)
                        CALL PLTFLU
                     ENDIF
                  ENDIF
                  GOTO 140
               ENDIF
  130       CONTINUE
  140       CONTINUE
         ENDIF
  150 CONTINUE

C  REDRAW THE POINT

      IF (.NOT.ALPHA)CALL PLTSTD (1, 3.)

      ADDLNK = .FALSE.
      CALL LTSORT (MP, LINKP, IMOVE, II, ADDLNK)
      IF (.NOT.ALPHA) THEN
#if NeedsDoubleEscape
         CALL MPD2SY (1, COOR (1, II), COOR (2, II), '\\CX')
#else
         CALL MPD2SY (1, COOR (1, II), COOR (2, II), '\CX')
#endif
         CALL PLTSTD (1, 7.)
         CALL PLTFLU
      ENDIF

      RETURN
      END
