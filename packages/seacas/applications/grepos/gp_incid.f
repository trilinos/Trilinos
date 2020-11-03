C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE INCID (TYPE, IDLST, NUMID, IDINC)
C=======================================================================

      CHARACTER*(*) TYPE
      INTEGER       IDLST(*)
      CHARACTER*80  STRING
      CHARACTER*8   STRA

C ... This routine increments all IDS by IDINC.

      IF (TYPE(:1) .EQ. 'M') THEN
         STRA = 'Material'
      ELSE IF (TYPE(:1) .EQ. 'S') THEN
         STRA = 'Sideset'
      ELSE IF (TYPE(:1) .EQ. 'N') THEN
         STRA = 'Nodeset'
      ELSE
         CALL PRTERR ('PROGRAM', 'unrecognized id type in INCID')
         RETURN
      END IF

C ... Just increment all ids by IDINC.  Note that we assume IDs valid
C     at entrance to routine, therefore, list will be valid at exit.

      IF (NUMID .LE. 0) RETURN
      DO 20 ID = 1, NUMID
         WRITE (STRING, 10) STRA, IDLST(ID), STRA, IDLST(ID)+IDINC
   10    FORMAT (A,1X,I10,' changed to ',A,1X,I10)
         CALL SQZSTR (STRING, LSTR)
         CALL PRTERR ('CMDSPEC', STRING(:LSTR))
         IDLST(ID) = IDLST(ID) + IDINC
   20 CONTINUE

      RETURN
      END
