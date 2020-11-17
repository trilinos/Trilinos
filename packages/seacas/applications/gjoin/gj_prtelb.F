C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE PRTELB (ISTAT, NELBL1, NELBL2,
     &   IDELB, NUMELB, NUMLNK, NUMATR, ISCR)
C=======================================================================

C   --*** PRTELB *** (GJOIN) Print summary for element blocks
C   --   Written by Amy Gilkey - revised 02/04/88
C   --
C   --PRTELB prints a summary of the element blocks, including their
C   --status.
C   --
C   --Parameters:
C   --   ISTAT - IN - the status of each element block:
C   --      0 = same
C   --      - = delete
C   --      n = combine with block n
C   --   NELBL1 - IN - the number of element blocks in the first database
C   --   NELBL2 - IN - the number of element blocks in the second database
C   --   IDELB - IN - the element block IDs for each block
C   --   NUMELB - IN - the number of elements in each block
C   --   NUMLNK - IN - the number of nodes per element in each block
C   --   NUMATR - IN - the number of attributes in each block
C   --   ISCR - SCRATCH - size = NELBL1+NELBL2

      INTEGER ISTAT(*)
      INTEGER IDELB(*)
      INTEGER NUMELB(*)
      INTEGER NUMLNK(*)
      INTEGER NUMATR(*)
      INTEGER ISCR(*)

      CHARACTER*2 STR2
      CHARACTER*20 STRA, STRB
      CHARACTER CH

      NELBLK = NELBL1 + NELBL2

      IF (NELBLK .LE. 0) RETURN

      WRITE (*, *)
      WRITE (*, *) 'Element blocks:'

      WRITE (STRA, 10000) NELBLK
10000  FORMAT ('(#', I5, ')')
      CALL PCKSTR (1, STRA)
      LSTRA = LENSTR (STRA)

      DO 110 IELB = 1, NELBLK
         IF (ISTAT(IELB) .EQ. 0) THEN
            WRITE (STRA, 10000) IELB
            CALL PCKSTR (1, STRA)
            CH = ' '
            IF (IELB .GT. NELBL1) CH = '*'
            WRITE (*, 10010) CH, IDELB(IELB), STRA(:LSTRA),
     &         NUMELB(IELB), NUMLNK(IELB)
         ELSE IF (ISTAT(IELB) .GT. 0) THEN
            IF (ISTAT(IELB) .EQ. IELB) THEN
               CALL GETALL (IELB, NELBLK, ISTAT, NINSET, ISCR)
               WRITE (STRB, '(I10)') IDELB(IELB)
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
                  IF (I .GT. NELBL1) CH = '*'
                  WRITE (*, 10010) CH, IDELB(I), STRA(:LSTRA),
     &               NUMELB(I), NUMLNK(I),
     &               STR2, 'combined into ID ', STRB(:LSTRB)
  100          CONTINUE
            ELSE
               CONTINUE
            END IF
         ELSE IF (ISTAT(IELB) .LT. 0) THEN
            WRITE (STRA, 10000) IELB
            CALL PCKSTR (1, STRA)
            CH = ' '
            IF (IELB .GT. NELBL1) CH = '*'
            WRITE (*, 10010) CH, IDELB(IELB), STRA(:LSTRA),
     &         NUMELB(IELB), NUMLNK(IELB), '<deleted>'
         END IF
  110 CONTINUE

      RETURN
10010  FORMAT (2X, A1, 1X, 'Block', I10, 1X, A, ':',
     &   I9, ' elements', I3, '-node', :, 1X, A, :, 1X, 10A)
      END
