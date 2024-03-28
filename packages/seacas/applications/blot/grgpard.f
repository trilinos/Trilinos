C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE GRGPARD (PARTYP, INDEV, IPARMS, IPSTR)
C=======================================================================
C     .. Same as GRGPAR, but only handles 'DEVICE' and 'SOFTCHAR'for PARTYP
C     .. this is to avoid issues with passing logical or integer for `iparms`

      include 'grpcom.blk'

      CHARACTER*(*) PARTYP
      INTEGER INDEV
      LOGICAL IPARMS
      CHARACTER*(*) IPSTR

      IF ((INDEV .NE. 1) .AND. (INDEV .NE. 2)) THEN
         IDEV = ICURDV
      ELSE
         IDEV = INDEV
      END IF

      IF (PARTYP .EQ. 'DEVICE') THEN
         IF (DEVOK(IDEV)) THEN
            CALL CPYLOG (1, .TRUE., IPARMS)
            IPSTR = DEVNAM(IDEV)
         ELSE
            CALL CPYLOG (1, .FALSE., IPARMS)
            IPSTR = ' '
         END IF
      ELSE IF (PARTYP .EQ. 'SOFTCHAR') THEN
         CALL CPYLOG (1, SOFTCH(IDEV), IPARMS)
         IF (SOFTCH(IDEV)) THEN
            IPSTR = 'software'
         ELSE
            IPSTR = 'hardware'
         END IF

      END IF

      RETURN
      END
