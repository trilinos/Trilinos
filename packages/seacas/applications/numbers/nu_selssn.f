C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE SELSSN (SELECT, NUMNP, NLIST, LIST,
     *   IDSS, NNSS, IPNSS, LTNSS, NUMSS, NUMSEL)
C=======================================================================
C    IDSS  (NUMSS) SIDE SET IDS
C    NNSS  (NUMSS) SIDE SET NODE    COUNTS
C    IPNSS (NUMSS) SIDE SET NODE    POINTERS
C    LTNSS (LSSNL) SIDE SET NODE    LIST

      LOGICAL SELECT(*)
      INTEGER LIST(*), IDSS(*), NNSS(*), IPNSS(*), LTNSS(*)
      CHARACTER*80 STRA

      CALL INILOG (NUMNP, .FALSE., SELECT)

      DO 20 II = 1, NLIST
         IFLG = LOCINT (LIST(II), NUMSS, IDSS)
         IF (IFLG .EQ. 0) THEN
            WRITE (STRA, 50) LIST(II)
            CALL SQZSTR (STRA, LSTR)
            CALL PRTERR ('ERROR', STRA(:LSTR))
         ELSE
            IBEG = IPNSS(IFLG)
            IEND = IBEG + NNSS(IFLG) - 1
            DO 10 I=IBEG, IEND
               SELECT(LTNSS(I)) = .TRUE.
   10       CONTINUE
         END IF
   20 CONTINUE

      NUMSEL = NUMEQL (.TRUE., NUMNP, SELECT)
      WRITE (STRA, 30) NUMSEL
      CALL SQZSTR(STRA, LSTR)
      WRITE (*, 40) STRA(:LSTR)
   30 FORMAT (I10,' nodes selected')
   40 FORMAT (/5X,A)

   50 FORMAT (' Set Flag ',I10,' not found. ')
      RETURN
      END
