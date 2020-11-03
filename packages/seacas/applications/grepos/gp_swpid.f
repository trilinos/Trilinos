C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE SWPID (TYPE, IDLST, NUMID, ID, DELOK)
C=======================================================================

      CHARACTER*(*) TYPE
      DIMENSION     IDLST(*)
      CHARACTER*80  STRING
      CHARACTER*8   STRA
      LOGICAL       DELOK

      IF (TYPE(:1) .EQ. 'M') THEN
         STRA = 'Material'
      ELSE IF (TYPE(:1) .EQ. 'S') THEN
         STRA = 'Sideset'
      ELSE IF (TYPE(:1) .EQ. 'N') THEN
         STRA = 'Nodeset'
      ELSE
         CALL PRTERR ('PROGRAM', 'unrecognized id type in SWPID')
         RETURN
      END IF

C ... Determine location of ID to be changed

      IMAT = LOCINT (ID, NUMID, IDLST)
      IF (IMAT .EQ. 0) THEN
         WRITE (STRING, 90) STRA, ID
   90    FORMAT (A,1X,I11,' does not exist')
         CALL SQZSTR (STRING, LSTR)
         CALL PRTERR ('ERROR', STRING(:LSTR))
         DELOK = .FALSE.
      ELSE
         WRITE (STRING, 100) STRA, ID
  100    FORMAT (A,1X,I11,' will be swapped')
         IDLST(IMAT) = -IDLST(IMAT)
         CALL SQZSTR (STRING, LSTR)
         CALL PRTERR ('CMDSPEC', STRING(:LSTR))
         DELOK = .TRUE.
      END IF

      RETURN
      END
