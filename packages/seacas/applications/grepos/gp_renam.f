C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE RENAM (TYPE, NAMLST, NUMNM, OLD, NEW)
C=======================================================================

      CHARACTER*(*) TYPE
      CHARACTER*(*) NAMLST(NUMNM), OLD, NEW
      CHARACTER*1024 STRING

C ... Determine location of NAME to be changed

      IMAT = LOCSTR (OLD, NUMNM, NAMLST)
      IF (IMAT .EQ. 0) THEN
         WRITE (STRING, 90) OLD, TYPE
   90    FORMAT (A,' does not exist in ',A,' variable list.')
         CALL SQZSTR (STRING, LSTR)
         CALL PRTERR ('ERROR', STRING(:LSTR))
         RETURN
      ELSE
         NAMLST(IMAT) = NEW
         WRITE (STRING, 100) TYPE, OLD, NEW
  100    FORMAT (A,' variable name ',A,' changed to ',A)
         CALL SQZSTR (STRING, LSTR)
         CALL PRTERR ('CMDSPEC', STRING(:LSTR))
      END IF

      RETURN
      END
