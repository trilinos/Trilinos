C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE NEWID (TYPE, IDLST, NUMID, IDNEW, IDOLD)
C=======================================================================

      CHARACTER*(*) TYPE
      DIMENSION IDLST(*)
      CHARACTER*80 STRING
      CHARACTER*8 STRA

      IF (TYPE(:1) .EQ. 'M') THEN
         STRA = 'Material'
      ELSE IF (TYPE(:1) .EQ. 'S') THEN
         STRA = 'Sideset'
      ELSE IF (TYPE(:1) .EQ. 'N') THEN
         STRA = 'Nodeset'
      ELSE
         CALL PRTERR ('PROGRAM', 'unrecognized id type in NEWID')
         RETURN
      END IF

C ... Check for existence of new id in old list - not allowed

      IF (NUMID .LE. 0) RETURN

      IMAT = LOCINT (IDNEW, NUMID, IDLST)

      IF (IMAT .NE. 0) THEN
         CALL PRTERR ('ERROR', 'Cannot change to an existing ID')
         RETURN
      END IF

C ... Determine location of ID to be changed

      IMAT = LOCINT (IDOLD, NUMID, IDLST)
      IF (IMAT .EQ. 0) THEN
         WRITE (STRING, 90) STRA, IDOLD
   90    FORMAT (A,1X,I5,' does not exist')
         CALL SQZSTR (STRING, LSTR)
         CALL PRTERR ('ERROR', STRING(:LSTR))
         RETURN
      ELSE
         IDLST(IMAT) = IDNEW
         WRITE (STRING, 100) STRA, IDOLD, STRA, IDNEW
  100    FORMAT (A,1X,I10,' changed to ',A,1X,I10)
         CALL SQZSTR (STRING, LSTR)
         CALL PRTERR ('CMDSPEC', STRING(:LSTR))
      END IF

      RETURN
      END
