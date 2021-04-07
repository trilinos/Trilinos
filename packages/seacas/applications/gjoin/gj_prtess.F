C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE PRTESS (ISTAT, NESS1, NESS2, IDESS, NEESS, ISCR)
C=======================================================================

C   --*** PRTESS *** (GJOIN) Print summary for element side sets
C   --   Written by Amy Gilkey - revised 02/04/88
C   --
C   --PRTESS prints a summary for the element side sets, including their
C   --status.
C   --
C   --Parameters:
C   --   ISTAT - IN - the status of each element side set:
C   --      0 = same
C   --      - = delete
C   --      n = combine with block n
C   --   NESS1 - IN - the number of element side sets in the first database
C   --   NESS2 - IN - the number of element side sets in the second database
C   --   IDESS - IN - the element side set ID for each set
C   --   NEESS - IN - the number of elements for each set
C   --   ISCR - SCRATCH - size = NESS1+NESS2

      INTEGER ISTAT(*)
      INTEGER IDESS(*)
      INTEGER NEESS(*)
      INTEGER ISCR(*)

      CHARACTER*2 STR2
      CHARACTER*20 STRA, STRB
      CHARACTER CH

      NUMESS = NESS1 + NESS2

      IF (NUMESS .LE. 0) RETURN

      WRITE (*, *)
      WRITE (*, *) 'Element side sets:'

      WRITE (STRA, 10000) NUMESS
10000  FORMAT ('(#', I5, ')')
      CALL PCKSTR (1, STRA)
      LSTRA = LENSTR (STRA)

      DO 110 IESS = 1, NUMESS
         IF (ISTAT(IESS) .EQ. 0) THEN
            WRITE (STRA, 10000) IESS
            CALL PCKSTR (1, STRA)
            CH = ' '
            IF (IESS .GT. NESS1) CH = '*'
            WRITE (*, 10010) CH, IDESS(IESS), STRA(:LSTRA),
     &         NEESS(IESS)
         ELSE IF (ISTAT(IESS) .GT. 0) THEN
            IF (ISTAT(IESS) .EQ. IESS) THEN
               CALL GETALL (IESS, NUMESS, ISTAT, NINSET, ISCR)
               WRITE (STRB, '(I5)') IDESS(IESS)
               CALL SQZSTR (STRB, LSTRB)
               DO 100 ISET = 1, NINSET
                  IF (ISET .EQ. 1) THEN
#if NeedsDoubleEscape
                     STR2 = '\\  '
#else
                     STR2 = '\ '
#endif
                  ELSE IF (ISET .EQ. NINSET) THEN
                     STR2 = '/ '
                  ELSE
                     STR2 = ' |'
                  END IF
                  I = ISCR(ISET)
                  WRITE (STRA, 10000) I
                  CALL PCKSTR (1, STRA)
                  CH = ' '
                  IF (I .GT. NESS1) CH = '*'
                  WRITE (*, 10010) CH, IDESS(I), STRA(:LSTRA),
     &               NEESS(I),
     &               STR2, 'combined into ID ', STRB(:LSTRB)
  100          CONTINUE
            ELSE
               CONTINUE
            END IF
         ELSE IF (ISTAT(IESS) .LT. 0) THEN
            WRITE (STRA, 10000) IESS
            CALL PCKSTR (1, STRA)
            CH = ' '
            IF (IESS .GT. NESS1) CH = '*'
            WRITE (*, 10010) CH, IDESS(IESS), STRA(:LSTRA),
     &         NEESS(IESS), '<deleted>'
         END IF
  110 CONTINUE

      RETURN
10010  FORMAT (2X, A1, 1X, 'Set', I6, 1X, A, ':',
     &   I9, ' elements and sides', :, 3X, A, :, 3X, 5A)
      END
