C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE SHODEV (SHOTYP)
C=======================================================================

C   --*** SHODEV *** (BLOT) Display device dependent options
C   --   Written by Amy Gilkey - revised 04/26/88
C   --
C   --SHODEV displays the device dependent options requested by the
C   --option type:
C   --   SOFTCHAR - software versus hardware characters
C   --   FONT - font to use
C   --   COLOR - number of standard colors to use and maximum on device
C   --   SPECTRUM - number of spectrum colors and maximum on device
C   --   SNAP - number of frames to snap
C   --   AUTO - automatic versus user-directed plotting
C   --
C   --Parameters:
C   --   SHOTYP - IN - the show option (see above)

C   --Routines Called:
C   --   GRGPAR - (GRPLIB) Get graphics device parameters
C   --   INTSTR - (STRLIB) Convert integers to strings
C   --   LENSTR - (STRLIB) Find string length
C   --   SQZSTR - (STRLIB) Delete extra blanks from string

      CHARACTER*(*) SHOTYP

      LOGICAL ONEDEV, TWODEV
      CHARACTER*3 DEVNAM(2)
      INTEGER IPARMS(3)
      CHARACTER*80 STRING
      CHARACTER*20 STR20

      CALL GRGPARD ('DEVICE', 1, ONEDEV, DEVNAM(1))
      CALL GRGPARD ('DEVICE', 2, TWODEV, DEVNAM(2))
      ISDEV = 1
      IF (.NOT. ONEDEV) ISDEV = 2
      IEDEV = 2
      IF (.NOT. TWODEV) IEDEV = 1
      TWODEV = (ONEDEV .AND. TWODEV)

      IF (SHOTYP .EQ. 'SOFTCHAR') THEN
         DO 100 I = ISDEV, IEDEV
            CALL GRGPAR (SHOTYP, I, IPARMS, STR20)
            L = LENSTR (STR20)
            STRING = 'Plot in ' // STR20(:L) // ' characters$'
            CALL PRTDEV (STRING, DEVNAM(I), TWODEV)
  100    CONTINUE

      ELSE IF (SHOTYP .EQ. 'FONT') THEN
         DO 110 I = ISDEV, IEDEV
            CALL GRGPAR (SHOTYP, I, IPARMS, STR20)
            L = LENSTR(STR20)
            STRING = 'Plot in ' // STR20(:L) // ' font$'
            CALL PRTDEV (STRING, DEVNAM(I), TWODEV)
  110    CONTINUE

      ELSE IF ((SHOTYP .EQ. 'COLOR') .OR. (SHOTYP .EQ. 'SPECTRUM')) THEN
         DO 120 I = ISDEV, IEDEV
            CALL GRGPAR ('COLOR', I, IPARMS, STR20)
            WRITE (STR20, 10000, IOSTAT=IDUM) IPARMS(1), IPARMS(2)
10000        FORMAT (I6, ' of ', I6)
            CALL SQZSTR (STR20, L)
            IF (IPARMS(2) .EQ. 0) THEN
               STRING = 'Number of colors to use$: single color device'
               CALL PRTDEV (STRING, DEVNAM(I), TWODEV)
            ELSE
               STRING = 'Number of standard colors to use$: '
     &            // STR20(:L)
               CALL PRTDEV (STRING, DEVNAM(I), TWODEV)
            END IF

            CALL GRGPAR ('SPECTRUM', I, IPARMS, STR20)
            IF (IPARMS(3) .GT. 0) THEN
               WRITE (STR20, 10000, IOSTAT=IDUM) IPARMS(1), IPARMS(2)
               CALL SQZSTR (STR20, L)
               IF (IPARMS(3) .EQ. 1) THEN
                  STRING = 'Number of spectrum colors to use$: '
     &               // STR20(:L)
                  CALL PRTDEV (STRING, DEVNAM(I), TWODEV)
               ELSE
                  STRING = '***$: ' // STR20(:L)
                  CALL PRTDEV (STRING, DEVNAM(I), TWODEV)
               END IF
            END IF
  120    CONTINUE

      ELSE IF (SHOTYP .EQ. 'SNAP') THEN
         I = ISDEV
         CALL GRGPAR (SHOTYP, I, IPARMS, STR20)
         IF (IPARMS(1) .LT. 0) THEN
            STRING =  'Cannot snap frames$'
            CALL PRTDEV (STRING, DEVNAM(I), .TRUE.)
         ELSE IF (IPARMS(1) .LE. 0) THEN
            STRING =  'Do not snap frames$'
            CALL PRTDEV (STRING, DEVNAM(I), TWODEV)
         ELSE
            CALL INTSTR (1, 0, IPARMS(1), STR20, L)
            STRING = 'Snap ' // STR20(:L) // ' frames for each plot$'
            CALL PRTDEV (STRING, DEVNAM(I), TWODEV)
         END IF

      ELSE IF (SHOTYP .EQ. 'AUTO') THEN
         I = ISDEV
         CALL GRGPAR (SHOTYP, I, IPARMS, STR20)
         L = LENSTR (STR20)
         STRING = 'Plot in ' // STR20(:L) // ' mode$'
         CALL PRTDEV (STRING, DEVNAM(I), TWODEV)
      END IF

      RETURN
      END
