C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE COMID (TYPE, IDLST, ISTAT, NUMID, IDFINAL, ID)
C=======================================================================

C   --   ISTAT - IN/OUT - the status of each item
C   --      0 = same
C   --      - = delete
C   --      n = combine with entity n

      CHARACTER*(*) TYPE
      DIMENSION     IDLST(*)
      INTEGER       ISTAT(*)
      CHARACTER*80  STRING
      CHARACTER*8   STRA

      IF (TYPE(:1) .EQ. 'M') THEN
         STRA = 'Material'
      ELSE IF (TYPE(:1) .EQ. 'S') THEN
         STRA = 'Sideset'
      ELSE IF (TYPE(:1) .EQ. 'N') THEN
         STRA = 'Nodeset'
      ELSE
         CALL PRTERR ('PROGRAM', 'unrecognized id type in COMID')
         RETURN
      END IF

C ... Check that IDFINAL and ID are different
      if (idfinal .eq. id) then
         write (string, 80) type(:lenstr(type)), idfinal
 80      format('Cannot combine ',A,1x,I11,' with itself. Ignoring')
         call sqzstr(string, lstr)
         call prterr ('CMDWARN', string(:lstr))
         return
      end if

C ... Determine location of FINAL ID (entity combined to)
      IMATF = LOCINT (IDFINAL, NUMID, IDLST)
      IF (IMATF .EQ. 0) THEN
         WRITE (STRING, 90) STRA, IDFINAL
         CALL SQZSTR (STRING, LSTR)
         CALL PRTERR ('ERROR', STRING(:LSTR))
         RETURN
      END IF

C ... Check that the entity with id IDFINAL is active
      if (istat(imatf) .ne. 0 .and. istat(imatf) .ne. imatf) then
         WRITE (STRING, 85) STRA, IDFINAL
 85      FORMAT (A,1X,I5,' has been modified. This is not allowed')
         CALL SQZSTR (STRING, LSTR)
         CALL PRTERR ('ERROR', STRING(:LSTR))
         RETURN
      END IF

C ... Determine location of ID to be changed
      IMAT  = LOCINT (ID, NUMID, IDLST)
      IF (IMAT .EQ. 0) THEN
         WRITE (STRING, 90) STRA, ID
   90    FORMAT (A,1X,I5,' does not exist')
         CALL SQZSTR (STRING, LSTR)
         CALL PRTERR ('ERROR', STRING(:LSTR))
      ELSE
         WRITE (STRING, 100) STRA, ID, STRA, IDFINAL
  100    FORMAT (A,1X,I11,' combined with ',A,1x,I11)
         ISTAT(IMAT)  = IMATF
         ISTAT(IMATF) = IMATF
         CALL SQZSTR (STRING, LSTR)
         CALL PRTERR ('CMDSPEC', STRING(:LSTR))
      END IF

      RETURN
      END
