C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE PRTNPS (ISTAT, NNPS1, NNPS2, IDNPS, NNNPS, ISCR)
C=======================================================================

C   --*** PRTNPS *** (GJOIN) Print summary for nodal point sets
C   --   Written by Amy Gilkey - revised 02/04/88
C   --
C   --PRTNPS prints a summary for the nodal point sets, including their
C   --status.
C   --
C   --Parameters:
C   --   ISTAT - IN - the status of each nodal point set:
C   --      0 = same
C   --      - = delete
C   --      n = combine with block n
C   --   NNPS1 - IN - the number of nodal point sets in the first database
C   --   NNPS2 - IN - the number of nodal point sets in the second database
C   --   IDNPS - IN - the nodal point set ID for each set
C   --   NNNPS - IN - the number of nodes for each set
C   --   ISCR - SCRATCH - size = NNPS1+NNPS2

      INTEGER ISTAT(*)
      INTEGER IDNPS(*)
      INTEGER NNNPS(*)
      INTEGER ISCR(*)

      CHARACTER*2 STR2
      CHARACTER*20 STRA, STRB
      CHARACTER CH

      NUMNPS = NNPS1 + NNPS2

      IF (NUMNPS .LE. 0) RETURN

      WRITE (*, *)
      WRITE (*, *) 'Nodal point sets:'

      WRITE (STRA, 10000) NUMNPS
10000  FORMAT ('(#', I5, ')')
      CALL PCKSTR (1, STRA)
      LSTRA = LENSTR (STRA)

      DO 110 INPS = 1, NUMNPS
         IF (ISTAT(INPS) .EQ. 0) THEN
            WRITE (STRA, 10000) INPS
            CALL PCKSTR (1, STRA)
            CH = ' '
            IF (INPS .GT. NNPS1) CH = '*'
            WRITE (*, 10010) CH, IDNPS(INPS), STRA(:LSTRA),
     &         NNNPS(INPS)
         ELSE IF (ISTAT(INPS) .GT. 0) THEN
            IF (ISTAT(INPS) .EQ. INPS) THEN
               CALL GETALL (INPS, NUMNPS, ISTAT, NINSET, ISCR)
               WRITE (STRB, '(I5)') IDNPS(INPS)
               CALL SQZSTR (STRB, LSTRB)
               DO 100 ISET = 1, NINSET
                  IF (ISET .EQ. 1) THEN
#if NeedsDoubleEscape
                     STR2 = '\\ '
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
                  IF (I .GT. NNPS1) CH = '*'
                  WRITE (*, 10010) CH, IDNPS(I), STRA(:LSTRA),
     &               NNNPS(I),
     &               STR2, 'combined into ID ', STRB(:LSTRB)
  100          CONTINUE
            ELSE
               CONTINUE
            END IF
         ELSE IF (ISTAT(INPS) .LT. 0) THEN
            WRITE (STRA, 10000) INPS
            CALL PCKSTR (1, STRA)
            CH = ' '
            IF (INPS .GT. NNPS1) CH = '*'
            WRITE (*, 10010) CH, IDNPS(INPS), STRA(:LSTRA),
     &         NNNPS(INPS), '<deleted>'
         END IF
  110 CONTINUE

      RETURN
10010  FORMAT (2X, A1, 1X, 'Set', I6, 1X, A, ':',
     &   I9, ' nodes', :, 3X, A, :, 3X, 5A)
      END
