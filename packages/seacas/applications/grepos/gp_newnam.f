C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE NEWNAM (TYPE, IDLST, NAMLST, NUMID, ID, NEWSTR)
C=======================================================================

      CHARACTER*(*) TYPE
      DIMENSION IDLST(*)
      CHARACTER*(*) NAMLST(NUMID), NEWSTR
      CHARACTER*1024 STRING
      CHARACTER*16 STRA, STRB

      STRB = 'name'
      IF (TYPE(:1) .EQ. 'B') THEN
         STRA = 'Element Block'
      ELSE IF (TYPE(:1) .EQ. 'b') THEN
         STRA = 'Material Block'
         STRB = 'type'
      ELSE IF (TYPE(:1) .EQ. 'C') THEN
         STRA = 'Coordinate'
      ELSE IF (TYPE(:1) .EQ. 'S') THEN
         STRA = 'Sideset'
      ELSE IF (TYPE(:1) .EQ. 'N') THEN
         STRA = 'Nodeset'
      ELSE
         CALL PRTERR ('PROGRAM', 'unrecognized id type in NEWNAM')
         RETURN
      END IF

      IF (NUMID .LE. 0) RETURN

C ... Determine location of ID to be changed

      IMAT = LOCINT (ID, NUMID, IDLST)
      IF (IMAT .EQ. 0) THEN
         WRITE (STRING, 90) STRA, ID
   90    FORMAT (A,1X,I5,' does not exist')
         CALL SQZSTR (STRING, LSTR)
         CALL PRTERR ('ERROR', STRING(:LSTR))
         RETURN
      ELSE
         NAMLST(IMAT) = NEWSTR
         WRITE (STRING, 100) STRA, ID, STRB, NEWSTR
  100    FORMAT (A,1X,I9,' ',A,' changed to ',A)
         CALL SQZSTR (STRING, LSTR)
         CALL PRTERR ('CMDSPEC', STRING(:LSTR))
      END IF

      RETURN
      END
