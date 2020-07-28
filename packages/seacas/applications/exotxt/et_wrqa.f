C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE WRQA (NTXT, NQAREC, QAREC, NINFO, INFO, QAINFO)
C=======================================================================

C   --*** WRQA *** (TXTEXO) Write QA and information records
C   --   Written by Amy Gilkey - revised 09/30/87
C   --
C   --WRQA writes the QA records and the information records.
C   --
C   --Parameters:
C   --   NTXT   - IN - the text file
C   --   NQAREC - IN - the number of QA records
C   --   QAREC  - IN - the QA records containing:
C   --                 (1) - the analysis code name
C   --                 (2) - the analysis code QA descriptor
C   --                 (3) - the analysis date
C   --                 (4) - the analysis time
C   --   NINFO  - IN - the number of information records
C   --   INFO   - IN - the information records
C   --
C   --Database must be positioned at start of QA records upon entry;
C   --upon exit positioned at end of information records.

      CHARACTER*32 QAREC(4,*)
      CHARACTER*32 QAINFO(6)
      CHARACTER*80 INFO(*)

      WRITE (NTXT, '(A)') '! QA Records'
      WRITE (NTXT, 10010) NQAREC+1, '! QA records'

      DO 100 IQA = 1, NQAREC
         WRITE (NTXT, '(A)') (QAREC(I,IQA), I=1,4)
  100 CONTINUE

C ... Add record for this code
      write (ntxt, '(A)') QAINFO(1), QAINFO(3), QAINFO(5)(:8),
     $     QAINFO(6)(:8)
      WRITE (NTXT, '(A)') '! Information Records'
      WRITE (NTXT, 10010) NINFO, '! information records'

      IF (NINFO .GT. 0) THEN
         DO 110 I = 1, NINFO
C ... Some codes are embedding newlines in the info records which screws up the text output
C     Filter them (and other strange characters) out here...
           do 105 j=1, 80
             if (ichar(info(i)(j:j)) .lt.  32 .or.
     *           ichar(info(i)(j:j)) .gt. 126) then
               info(i)(j:j) = ' '
             end if
 105         continue
            WRITE (NTXT, 10000) INFO(I)
  110    CONTINUE
      END IF

      RETURN
10000  FORMAT (99 (A, :, 1X))
10010  FORMAT (1I10, 6X, A)
      END
