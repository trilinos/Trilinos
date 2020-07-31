C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE DBOQA (NDB, NQAREC, QAREC, NINFO, INFO)
C=======================================================================

C   --*** DBOQA *** (EXOLIB) Write QA and information records
C   --   Written by Amy Gilkey - revised 02/08/88
C   --
C   --DBOQA writes the QA records and the information records.
C   --
C   --Parameters:
C   --   NDB - IN - the database number
C   --   NQAREC - IN - the number of QA records (must be at least one);
C   --      written only if >= 0
C   --   QAREC - IN - the QA records containing:
C   --      (1) - the analysis code name
C   --      (2) - the analysis code QA descriptor
C   --      (3) - the analysis date
C   --      (4) - the analysis time
C   --   NINFO - IN - the number of information records; written only if >= 0
C   --   INFO - IN - the information records
C   --
C   --Database must be positioned at start of QA records upon entry;
C   --upon exit positioned at end of information records.

      INTEGER NDB
      INTEGER NQAREC
      CHARACTER*8 QAREC(4,*)
      INTEGER NINFO
      CHARACTER*80 INFO(*)

      IF (NQAREC .LT. 0) GOTO 120

      WRITE (NDB) NQAREC

      DO 100 IQA = 1, NQAREC
         WRITE (NDB) (QAREC(I,IQA), I=1,4)
  100 CONTINUE

      IF (NINFO .LT. 0) GOTO 120

      WRITE (NDB) NINFO

      DO 110 I = 1, NINFO
         WRITE (NDB) INFO(I)
  110 CONTINUE

  120 CONTINUE
      RETURN
      END
