C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details
C=======================================================================
      SUBROUTINE DBIQA (NDB, OPTION, NQAREC, QAREC, NINFO, INFO)
C=======================================================================

C   --*** DBIQA *** (EXOLIB) Read QA and information records
C   --   Written by Amy Gilkey - revised 02/08/88
C   --   Modified for ExodusIIV2 database format - 9/10/95
C   --
C   --Parameters:
C   --   NDB    - IN  - the database number
C   --   OPTION - IN  - ' ' to not store, '*' to store all, else store options:
C   --                  'Q' to store QA records
C   --                  'I' to store information records
C   --   NQAREC - IN  - the number of QA records; <0 if end-of-file
C   --   QAREC  - OUT - the QA records containing: (if OPTION)
C   --                  (1) - the analysis code name
C   --                  (2) - the analysis code QA descriptor
C   --                  (3) - the analysis date
C   --                  (4) - the analysis time
C   --   NINFO  - IN  - the number of information records; <0 if end-of-file
C   --   INFO   - OUT - the information records (if OPTION)

      include 'exodusII.inc'

      INTEGER NDB
      CHARACTER*(*) OPTION
      INTEGER NQAREC
      CHARACTER*(MXSTLN) QAREC(4,*)
      INTEGER NINFO
      CHARACTER*(MXLNLN) INFO(*)
      INTEGER IERR
      LOGICAL ALL

      ALL  = (OPTION .EQ. '*')

C     Read in Qa Records
      IF ((ALL .OR. (INDEX(OPTION, 'Q') .GT. 0))
     &                    .AND. (NQAREC .GT. 0)) THEN
         CALL EXGQA(NDB, QAREC, IERR)
      END IF

C     Read in Information Records
      IF ((ALL .OR. (INDEX(OPTION, 'I') .GT. 0))
     &                     .AND. (NINFO .GT. 0)) THEN
         CALL EXGINF(NDB, INFO, IERR)
      END IF

      RETURN
      END
