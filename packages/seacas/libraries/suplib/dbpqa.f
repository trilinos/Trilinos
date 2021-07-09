C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details
C=======================================================================
      SUBROUTINE DBPQA (OPTION, NQAREC, QAREC, NINFO, INFO)
C=======================================================================
C   --*** DBPQA *** (EXOLIB) Print QA and information records
C   --   Written by Amy Gilkey - revised 02/02/88
C   --
C   --DBPQA displays the QA records and the information records.
C   --
C   --Parameters:
C   --   OPTION - IN - '*' to print all, else print options:
C   --                 'Q' to print QA records
C   --                 'I' to print information records
C   --   NQAREC - IN - the number of QA records (if OPTION)
C   --   QAREC  - IN - the QA records containing: (if OPTION)
C   --                 (1) - the analysis code name
C   --                 (2) - the analysis code QA descriptor
C   --                 (3) - the analysis date
C   --                 (4) - the analysis time
C   --   NINFO - IN - the number of information records (if OPTION)
C   --   INFO  - IN - the information records (if OPTION)

      CHARACTER*(*) OPTION
      INTEGER NQAREC
      CHARACTER*(32) QAREC(4,*)
      INTEGER NINFO
      CHARACTER*(80) INFO(*)

      LOGICAL ALL

      ALL = OPTION .EQ. '*'

      IF (ALL .OR. (INDEX (OPTION, 'Q') .GT. 0)) THEN
         IF (NQAREC .GT. 0) THEN
            WRITE (*, *)
            DO 100 IQA = 1, NQAREC
               call sqzstr(qarec(1, iqa), icode)
               call sqzstr(qarec(2, iqa), ivers)
               call sqzstr(qarec(3, iqa), idate)
               call sqzstr(qarec(4, iqa), itime)
               WRITE (*, 10000)qarec(1,iqa)(:icode),
     &         qarec(2, iqa)(:max(ivers,11)), qarec(3, iqa)(:idate),
     &         qarec(4, iqa)(:itime)
10000           FORMAT (1X, 'Code:  ', A, '  version  ', A,
     &                  '  on  ', A, '  at  ', A)
  100       CONTINUE
         END IF
      END IF

      IF (ALL .OR. (INDEX (OPTION, 'I') .GT. 0)) THEN
         IF (NINFO .GT. 0) THEN
            WRITE (*, *)
            DO 110 I = 1, NINFO
               WRITE (*, '(1x,a)') INFO(I)(:LENSTR(INFO(I)))
  110       CONTINUE
         END IF
      END IF

      RETURN
      END
