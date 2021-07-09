C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE ERASEL (MP, ML, MS, MR, N, COOR, ILINE, LTYPE, LCON,
     &   NLPS, IFLINE, ILLIST, NSPR, IFSIDE, ISLIST, LINKP, LINKL,
     &   LINKS, LINKR, IBAD, ALPHA)
C***********************************************************************

C  SUBROUTINE ERASEL = ERASES A LINE, AND ALL DATA DEPENDENT ON THAT
C                      LINE (I.E. SIDES OR REGIONS)

C***********************************************************************

      DIMENSION COOR(2, MP), ILINE(ML), LTYPE(ML), LCON(3, ML), NLPS(MS)
      DIMENSION NSPR(MR), IFLINE(MS), ILLIST(MS*3), IFSIDE(MR)
      DIMENSION ISLIST(MR*4), LINKP(2, MP), LINKL(2, ML), LINKS(2, MS)
      DIMENSION LINKR(2, MR), N(29)

      LOGICAL NUMPLT, ALPHA, ADDLNK, TEST, GETMAX

      IZ = 0
      KNUM = 0
      ADDLNK = .FALSE.
      TEST = .FALSE.
      GETMAX = .FALSE.

C  ERASE ANY REGIONS USING THAT LINE

      DO 120 I = 1, N(22)
         CALL LTSORT(MR, LINKR, I, II, ADDLNK)
         IF (II .GT. 0) THEN
            DO 110 J = IFSIDE(II), IFSIDE(II) + NSPR(II) - 1
               CALL LTSORT(MS, LINKS, ISLIST(J), IPNTR, ADDLNK)
               IF ((ISLIST(J) .GT. 0) .AND. (IPNTR .GT. 0)) THEN
                  JJ = ISLIST(J)
                  DO 100 K = IFLINE(JJ), IFLINE(JJ) + NLPS(JJ) - 1
                     KK = ILLIST(K)
                     IF (KK .EQ. IBAD) THEN
                        ADDLNK = .TRUE.
                        CALL LTSORT(MR, LINKR, I, IZ, ADDLNK)
                        ADDLNK = .FALSE.
                     END IF
  100             CONTINUE
               ELSE IF (ISLIST(J) .LT. 0) THEN
                  JJ = IABS(ISLIST(J))
                  IF (JJ .EQ. IBAD) THEN
                     ADDLNK = .TRUE.
                     CALL LTSORT(MR, LINKR, I, IZ, ADDLNK)
                     ADDLNK = .FALSE.
                  END IF
               END IF
  110       CONTINUE

            CALL LTSORT(MR, LINKR, I, IPNTR, ADDLNK)
            IF ((.NOT.ALPHA) .AND. (IPNTR .LE. 0)) THEN
               CALL REGEXT (MP, ML, MS, MR, N, II, COOR, ILINE,
     &            LTYPE, LCON, NLPS, IFLINE, ILLIST, NSPR, IFSIDE,
     &            ISLIST, LINKP, LINKL, LINKS, LINKR, XMIN, XMAX,
     &            YMIN, YMAX)
               XMID = (XMAX + XMIN)/2.
               YMID = (YMAX + YMIN)/2.
#if NeedsDoubleEscape
               CALL MPD2SY(1, XMID, YMID, '\\CDI')
#else
               CALL MPD2SY(1, XMID, YMID, '\CDI')
#endif
               CALL PLTFLU
            END IF
         END IF
  120 CONTINUE

C  ERASE ANY SIDES USING THAT LINE

      DO 150 I = 1, N(20)
         CALL LTSORT(MS, LINKS, I, II, ADDLNK)
         IF (II .GT. 0) THEN
            DO 130 J = IFLINE(II), NLPS(II) + IFLINE(II) - 1
               JJ = ILLIST(J)
               IF (JJ .EQ. IBAD) THEN
                  ADDLNK = .TRUE.
                  CALL LTSORT(MS, LINKS, I, IZ, ADDLNK)
                  ADDLNK = .FALSE.
                  GO TO 140
               END IF
  130       CONTINUE
  140       CONTINUE
         END IF
  150 CONTINUE

C  DELETE THE LINE

      CALL LTSORT(ML, LINKL, IBAD, II, ADDLNK)
      IF (II .GT. 0) THEN
         LT = LTYPE(II)
         IP1 = LCON(1, II)
         IP2 = LCON(2, II)
         IP3 = LCON(3, II)
         CALL LTSORT(MP, LINKP, IP1, IPNTR1, ADDLNK)
         CALL LTSORT(MP, LINKP, IP2, IPNTR2, ADDLNK)
         IF (IP3 .NE. 0) THEN
            CALL LTSORT(MP, LINKP, IABS(IP3), IPNTR3, ADDLNK)
         ELSE
            IPNTR3 = 0
         END IF
         IF ((IPNTR1 .GT. 0) .AND. (IPNTR2 .GT. 0) .AND.
     &      ((LT .EQ. 1) .OR. (IPNTR3 .GT. 0))) THEN
            IF (.NOT.ALPHA) THEN
               CALL DLINE(MP, ML, COOR, LINKP, KNUM, LT, IP1, IP2,
     &            IP3, NUMPLT, DUM1, DUM2, TEST, GETMAX, DUM3, DUM4,
     &            DUM5, DUM6)
               CALL PLTFLU
            END IF
         END IF
         ADDLNK = .TRUE.
         CALL LTSORT(ML, LINKL, IBAD, IZ, ADDLNK)
         ADDLNK = .FALSE.
      END IF

      RETURN

      END
