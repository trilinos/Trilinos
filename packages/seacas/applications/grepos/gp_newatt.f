C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE NEWATT (IDLST, ID, NUM, IDATT, NUMATR, ATTNAM, NEWNAM)
C=======================================================================

      include 'gp_namlen.blk'
      INTEGER IDLST(*)
      INTEGER NUMATR(*)
      CHARACTER*(maxnam) ATTNAM(*), NEWNAM
      CHARACTER*1024 STRING
      CHARACTER*16 STRA

      STRA = 'Element Block'

      IF (NUM .LE. 0) RETURN

C ... Determine location of ID to be changed

      IMAT = LOCINT (ID, NUM, IDLST)
      IF (IMAT .EQ. 0) THEN
        WRITE (STRING, 90) STRA, ID
 90     FORMAT (A,1X,I5,' does not exist')
        CALL SQZSTR (STRING, LSTR)
        CALL PRTERR ('ERROR', STRING(:LSTR))
        RETURN
      END IF

C ... Determine beginning of attribute names for this block
      IOFF = 0
      DO I=1, IMAT-1
        IOFF = IOFF + NUMATR(I)
      END DO

      ATTNAM(IOFF + IDATT) = NEWNAM
      ATTNAM(IOFF + IDATT) = NEWNAM
      WRITE (STRING, 100) IDATT, STRA, ID, NEWNAM
 100  FORMAT ('Name of Attribute ', I9, ' on ', A,1X,I9,
     *  ' changed to ',A)
         CALL SQZSTR (STRING, LSTR)
         CALL PRTERR ('CMDSPEC', STRING(:LSTR))

      RETURN
      END
