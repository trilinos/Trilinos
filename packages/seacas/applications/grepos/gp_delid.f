C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE DELID (TYPE, IDLST, ISTAT, NUMID, ID, DODEL)
C=======================================================================

C   --   ISTAT - IN/OUT - the status of each item
C   --      0 = same
C   --      - = delete
C   --      n = combine with entity n

      CHARACTER*(*) TYPE
      INTEGER  IDLST(*)
      INTEGER  ISTAT(*)
      LOGICAL  DODEL
      CHARACTER*80  STRING
      CHARACTER*8   STRA

      IF (TYPE(:1) .EQ. 'M') THEN
         STRA = 'Material'
      ELSE IF (TYPE(:1) .EQ. 'S') THEN
         STRA = 'Sideset'
      ELSE IF (TYPE(:1) .EQ. 'N') THEN
         STRA = 'Nodeset'
      ELSE
         CALL PRTERR ('PROGRAM', 'unrecognized id type in DELID')
         RETURN
      END IF

C ... Determine location of ID to be changed

      IMAT = LOCINT (ID, NUMID, IDLST)
      IF (IMAT .EQ. 0) THEN
         WRITE (STRING, 90) STRA, ID
   90    FORMAT (A,1X,I11,' does not exist')
         CALL SQZSTR (STRING, LSTR)
         CALL PRTERR ('ERROR', STRING(:LSTR))
      ELSE
        if (dodel) then
         WRITE (STRING, 100) STRA, ID, ' deleted'
       else
         WRITE (STRING, 100) STRA, ID, ' undeleted'
       end if
  100    FORMAT (A,1X,I11,A)
         if (idlst(imat) .eq. 0) then
           istat(imat) = -9999
         else
           if (dodel) then
             ISTAT(IMAT) = -IDLST(IMAT)
           else
             ISTAT(IMAT) = 0
           end if
         end if
         CALL SQZSTR (STRING, LSTR)
         CALL PRTERR ('CMDSPEC', STRING(:LSTR))
      END IF

      RETURN
      END
