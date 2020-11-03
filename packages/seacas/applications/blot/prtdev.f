C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE PRTDEV (STRING, DEVNAM, INCDEV)
C=======================================================================

C   --*** PRTDEV *** (BLOT) Internal to SHODEV
C   --   Written by Amy Gilkey - revised 06/01/87
C   --
C   --PRTDEV prints a device-dependent parameter string, including the
C   --device name if there is more than one device.
C   --
C   --Parameters:
C   --   STRING - IN - the parameter string with '$' where ' on DEV' belongs
C   --   DEVNAM - IN - the device name
C   --   INCDEV - IN - true iff device name is to be included

C   --Routines Called:
C   --   LENSTR - (STRLIB) Find string length

      CHARACTER*(*) STRING
      CHARACTER*(*) DEVNAM
      LOGICAL INCDEV

      LOGICAL DEVOK

      L = LENSTR (STRING)
      I = INDEX (STRING, '$')
      IF (I .EQ. 0) I = L+1

      DEVOK = (DEVNAM .NE. '$$$')

      IF (I .GE. L) THEN
         IF (.NOT. (INCDEV .AND. DEVOK)) THEN
            WRITE (*, 10000) STRING(1:I-1)
         ELSE
            WRITE (*, 10000) STRING(1:I-1), ' on ', DEVNAM
         END IF
      ELSE
         IF (.NOT. (INCDEV .AND. DEVOK)) THEN
            WRITE (*, 10000) STRING(1:I-1), STRING(I+1:L)
         ELSE
            WRITE (*, 10000) STRING(1:I-1), ' on ', DEVNAM,
     &         STRING(I+1:L)
         END IF
      END IF

      RETURN
10000  FORMAT (1X, 6A)
      END
