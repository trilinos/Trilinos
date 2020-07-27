C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details
C=======================================================================
      SUBROUTINE DBOQA (NDB, QAINFO, NQAREC, QAREC, NINFO, INFO)
C=======================================================================

C   --*** DBOQA *** (EXOLIB) Write QA and information records
C   --   Written by Amy Gilkey - revised 02/08/88
C   --   Modified by for EXODUSIIV2- 9/11/95
C   --
C   --Parameters:
C   --   NDB    - IN - the database number
C   --   QAINFO - IN - the QA record for the current run of algebra
C   --   NQAREC - IN - the number of QA records written only if >= 0
C   --   QAREC  - IN - the QA records containing:
C   --                 (1) - the analysis code name
C   --                 (2) - the analysis code QA descriptor
C   --                 (3) - the analysis date
C   --                 (4) - the analysis time
C   --   NINFO  - IN - the number of information records; written only if >= 0
C   --   INFO   - IN - the information records

      include 'exodusII.inc'

      INTEGER NDB
      CHARACTER*(MXSTLN)  QAINFO(6)
      CHARACTER*(MXSTLN)  QAREC(4,*)
      CHARACTER*(MXLNLN)  INFO(*)

      nqarec = nqarec + 1
      QAREC(1,nqarec) = QAINFO(1)
      QAREC(2,nqarec) = QAINFO(3)
      QAREC(3,nqarec) = QAINFO(5)
      QAREC(4,nqarec) = QAINFO(6)

C     Write QA records to the file ndbout
C     There is always at least one QA record
      call expqa(ndb, nqarec, qarec, ierr)

C     Write information records to the file ndbout
      IF (NINFO .GT. 0) THEN
         call expinf(ndb, ninfo, info, ierr)
      END IF

      RETURN
      END
