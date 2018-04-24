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

C $Log: grsnap.f,v $
C Revision 1.3  2009/03/25 12:36:44  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.2  1998/06/12 15:53:15  gdsjaar
C 1. Problem with TIMES array. Blot accesses a dummy timestep even if
C there were no timesteps on the database. Array wasn't allocated, so
C writing off into never-never land.
C
C 2. Inconsistency among some common blocks. Some places weren't using
C the include but had the definition hardwired in. Removed those.
C
C 3. Added 'EXTERNAL BLKDAT' to all routines that used data values set
C in BLKDAT
C
C 4. Cleanup of some A vs. IA argument passing.
C
C Revision 1.1  1994/04/07 20:02:46  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
CRevision 1.2  1990/12/14  08:51:54  gdsjaar
CAdded RCS Id and Log to all files
C
C=======================================================================
      SUBROUTINE GRSNAP (CMD, INDEV)
C=======================================================================

C   --*** GRSNAP *** (GRPLIB) Perform movie snap operations
C   --   Written by Amy Gilkey, revised 01/25/88
C   --
C   --GRSNAP sets up a device for snapping multiple frames for movies.
C   --This code is very device-dependent.
C   --
C   --At the present time, the only device that may be snapped is a camera
C   --(attached to a terminal) or a Dicomed.
C   --
C   --An attached camera device must be pre-initialized and NSNAP must
C   --be set non-negative (zero ok).
C   --
C   --Parameters:
C   --   CMD - IN - the snap operation:
C   --      QUERY = initialize device code, nsnap < 0 if cannot snap frames
C   --      INIT  = initialize device to snap nsnap frames
C   --      ABORT = abort plot
C   --      START = start plot
C   --      STOP  = end plot and snap frames
C   --      EXIT  = deactivate device
C   --   INDEV - IN - the device; 0 for current device (only QUERY allowed
C   --      for non-current device)
C   --
C   --Common Variables:
C   --   Uses ICURDV, DEVCOD, NSNAP of /GRPCOM/

C   --Routines Called:
C   --   VDESCP - (VDI) Send escape sequence to device

      COMMON /GRPCOC/ DEVNAM(2), DEVCOD(2)
      CHARACTER*3 DEVNAM
      CHARACTER*8 DEVCOD
      COMMON /GRPCOM/ ICURDV, ISHARD, DEVOK(2), TALKOK(2),
     &   NSNAP(2), IFONT(2), SOFTCH(2), AUTOPL(2),
     &   MAXCOL(2), NUMCOL(0:1,2), MAPALT(2), MAPUSE(2)
      LOGICAL ISHARD, DEVOK, TALKOK, SOFTCH, AUTOPL

      CHARACTER*(*) CMD
      INTEGER INDEV

      INTEGER RBUF(2)

      CHARACTER*8 IDSNAP(2)
      SAVE IDSNAP
C      --IDSNAP - an artificial device code related to the device code

      DATA IDSNAP / ' ', ' ' /

      IF ((INDEV .NE. 1) .AND. (INDEV .NE. 2)) THEN
         IDEV = ICURDV
      ELSE
         IDEV = INDEV
      END IF

C   --Initialize device code
      IF (IDSNAP(IDEV) .EQ. ' ') THEN
         IDSNAP(IDEV) = 'NONE'
         IF (NSNAP(IDEV) .GE. 0) THEN
            IDSNAP(IDEV) = 'CAMERA'
         ELSE IF (DEVCOD(IDEV) .EQ. 'DICOMED') THEN
            IDSNAP(IDEV) = DEVCOD(IDEV)
         END IF

C      --Initialize the device
         IF (IDSNAP(IDEV) .EQ. 'CAMERA') THEN
C         --Initialized before the routine was called
            CONTINUE
         ELSE IF (IDSNAP(IDEV) .EQ. 'DICOMED') THEN
            CONTINUE
         END IF
      END IF

C   --Set nsnap = -999 if cannot snap device
      IF (IDSNAP(IDEV) .EQ. 'NONE') NSNAP(IDEV) = -999

      IF (CMD .NE. 'QUERY') THEN

C      --Skip if not current device
         IF (IDEV .NE. ICURDV) GOTO 110

C      --Skip if snaps are not requested on device
         IF (NSNAP(ICURDV) .LE. 0) GOTO 110

         CALL PLTFLU

         IF (CMD .EQ. 'INIT') THEN

C         --Initialize number of snaps

            IF (IDSNAP(IDEV) .EQ. 'CAMERA') THEN
               CONTINUE

            ELSE IF (IDSNAP(IDEV) .EQ. 'DICOMED') THEN
               CONTINUE
            END IF

         ELSE IF (CMD .EQ. 'ABORT') THEN

C         --Abort plot

            IF (IDSNAP(IDEV) .EQ. 'CAMERA') THEN
               CONTINUE

            ELSE IF (IDSNAP(IDEV) .EQ. 'DICOMED') THEN
C            --Close segment and delete all segments
               CALL VDESCP (201, 0, 0)
               RBUF(1) = 0
               CALL VDESCP (203, 1, RBUF)
            END IF

         ELSE IF (CMD .EQ. 'START') THEN

C         --Start plot

            IF (IDSNAP(IDEV) .EQ. 'CAMERA') THEN
               CONTINUE

            ELSE IF (IDSNAP(IDEV) .EQ. 'DICOMED') THEN
C            --Open segment
               RBUF(1) = 1
               RBUF(2) = 1
               CALL VDESCP (200, 2, RBUF)
            END IF

         ELSE IF (CMD .EQ. 'STOP') THEN

C         --End plot and snap frames

            IF (IDSNAP(IDEV) .EQ. 'CAMERA') THEN
C            --Snap n frames
               CONTINUE

            ELSE IF (IDSNAP(IDEV) .EQ. 'DICOMED') THEN
C            --Close segment
               CALL VDESCP (201, 0, 0)

C            --Snap n-1 frames (plot segment with newpage)
               DO 100 I = 1, NSNAP(ICURDV)-1
                  CALL VDNWPG
  100          CONTINUE

C            --Delete all segments
               RBUF(1) = 0
               CALL VDESCP (203, 1, RBUF)
            END IF

         ELSE IF (CMD .EQ. 'EXIT') THEN

C         --Deactivate device

            IF (IDSNAP(IDEV) .EQ. 'CAMERA') THEN
               CONTINUE

            ELSE IF (IDSNAP(IDEV) .EQ. 'DICOMED') THEN
               CONTINUE
            END IF

         ELSE
            GOTO 110
         END IF
      END IF

  110 CONTINUE
      RETURN
      END
