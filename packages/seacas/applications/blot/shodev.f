C Copyright(C) 2009-2017 National Technology & Engineering Solutions of
C Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C Redistribution and use in source and binary forms, with or without
C modification, are permitted provided that the following conditions are
C met:
C
C     * Redistributions of source code must retain the above copyright
C       notice, this list of conditions and the following disclaimer.
C
C     * Redistributions in binary form must reproduce the above
C       copyright notice, this list of conditions and the following
C       disclaimer in the documentation and/or other materials provided
C       with the distribution.
C     * Neither the name of NTESS nor the names of its
C       contributors may be used to endorse or promote products derived
C       from this software without specific prior written permission.
C
C THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
C "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
C LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
C A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
C OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
C SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
C LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
C DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
C THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
C (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
C OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

C $Log: shodev.f,v $
C Revision 1.2  2009/03/25 12:36:48  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.1  1994/04/07 20:12:10  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:57:20  gdsjaar
c Added RCS Id and Log to all files
c
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

      CALL GRGPAR ('DEVICE', 1, ONEDEV, DEVNAM(1))
      CALL GRGPAR ('DEVICE', 2, TWODEV, DEVNAM(2))
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
