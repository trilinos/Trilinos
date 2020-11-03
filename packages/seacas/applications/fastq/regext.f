C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE REGEXT (MP, ML, MS, MR, N, II, COOR, ILINE, LTYPE,
     &   LCON, NLPS, IFLINE, ILLIST, NSPR, IFSIDE, ISLIST, LINKP,
     &   LINKL, LINKS, LINKR, XMIN, XMAX, YMIN, YMAX)
C***********************************************************************

C  SUBROUTINE REGEXT = GETS THE REGION EXTREMES

C***********************************************************************

      DIMENSION COOR (2, MP)
      DIMENSION ILINE (ML), LTYPE (ML), LCON (3, ML)
      DIMENSION NLPS (MS), IFLINE (MS), ILLIST (MS * 3)
      DIMENSION NSPR (MR), IFSIDE (MR), ISLIST (MR * 4)
      DIMENSION LINKP (2, MP), LINKL (2, ML), LINKS (2, MS)
      DIMENSION LINKR (2, MR)
      DIMENSION N (29)

      LOGICAL FOUND, GETMAX, ADDLNK
      LOGICAL NUMPLT, TEST

      ADDLNK = .FALSE.
      GETMAX = .TRUE.
      FOUND = .FALSE.

      DO 110 J = IFSIDE (II), IFSIDE (II) + NSPR (II) - 1

C  GET SIDE EXTREMES

         IF ( ISLIST(J) .GT. 0) THEN
            CALL LTSORT (MS, LINKS, ISLIST (J), IPNTR, ADDLNK)
            IF (IPNTR .GT. 0) THEN
               DO 100 K = IFLINE (IPNTR), IFLINE (IPNTR) +
     &              NLPS (IPNTR) - 1
                 CALL LTSORT (ML, LINKL, ILLIST (K), KK, ADDLNK)
                 IF (KK .GT. 0) THEN
                    IF (.NOT.FOUND) THEN
                       CALL LTSORT (MP, LINKP, IABS (LCON (1, KK)),
     &                      IPNT, ADDLNK)
                       IF (IPNT .GT. 0) THEN
                          XMAX = COOR (1, IPNT)
                          XMIN = COOR (1, IPNT)
                          YMAX = COOR (2, IPNT)
                          YMIN = COOR (2, IPNT)
                          FOUND = .TRUE.
                       ENDIF
                    ENDIF
                    IF (FOUND) THEN
                       CALL DLINE (MP, ML, COOR, LINKP, ILINE (KK),
     &                      LTYPE (KK), LCON (1, KK), LCON (2, KK),
     &                      LCON (3, KK), NUMPLT, X1, Y1, TEST, GETMAX,
     &                      XMIN, XMAX, YMIN, YMAX)
                    ENDIF
                 ENDIF
 100           CONTINUE
            END IF

C  GET LINE EXTREMES

         ELSEIF (ISLIST (J) .LT. 0) THEN
            JJ = IABS (ISLIST (J))
            CALL LTSORT (ML, LINKL, JJ, KK, ADDLNK)
            IF (KK .GT. 0) THEN
               IF (.NOT.FOUND) THEN
                  CALL LTSORT (MP, LINKP, IABS (LCON (1, KK)), IPNT,
     &               ADDLNK)
                  IF (IPNT .GT. 0) THEN
                     XMAX = COOR (1, IPNT)
                     XMIN = COOR (1, IPNT)
                     YMAX = COOR (2, IPNT)
                     YMIN = COOR (2, IPNT)
                     FOUND = .TRUE.
                  ENDIF
               ENDIF
               IF (FOUND) THEN
                  CALL DLINE (MP, ML, COOR, LINKP, ILINE (KK),
     &               LTYPE (KK), LCON (1, KK), LCON (2, KK),
     &               LCON (3, KK), NUMPLT, X1, Y1, TEST, GETMAX, XMIN,
     &               XMAX, YMIN, YMAX)
               ENDIF
            ENDIF
         ENDIF
  110 CONTINUE

      RETURN

      END
