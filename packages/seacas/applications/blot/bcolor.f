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

C $Log: bcolor.f,v $
C Revision 1.3  2009/03/25 12:36:42  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.2  1994/06/13 18:29:00  gdsjaar
C Fixed background and foreground color setting. (I think)
C
c Revision 1.1  1994/04/07  19:54:46  gdsjaar
c Initial checkin of ACCESS/graphics/blotII2
c
c Revision 1.2  1990/12/14  08:47:36  gdsjaar
c Added RCS Id and Log to all files
c
C============================================================================
      SUBROUTINE BCOLOR (INIT, INLINE, IFLD, INTYP, IFIELD,
     &   CFIELD, BLKCOL)
C============================================================================

C   --*** BCOLOR ***  (BLOT) Process BLKCOL command
C   --   Written by John Glick - 11/22/88
C   --
C   --Parameters:
C   --   INIT  - IN - .TRUE. iff the purpose of the call is to
C   --            initialize the BLKCOL array and not to modify it.
C   --   INLINE - IN/OUT - the parsed input lines for the log file
C   --   IFLD, INTYP, IFIELD, CFIELD, - IN/OUT - the free-field reader
C   --          index and charcter field.
C   --   BLKCOL - IN/OUT - the user selected colors of the element blocks.
C   --                    BLKCOL(0) = 1 if the user defined material
C   --                                colors should be used in mesh plots.
C   --                              = -1 if program selected colors should
C   --                                be used.
C   --                    BLKCOL(i) = the user selected color of element
C   --                               block i:
C   --                                  -2 - no color selected by user.
C   --                                  -1 - black
C   --                                   0 - white
C   --                                   1 - red
C   --                                   2 - green
C   --                                   3 - yellow
C   --                                   4 - blue
C   --                                   5 - cyan
C   --                                   6 - magenta
      LOGICAL INIT
      CHARACTER*(*) INLINE(*)
      INTEGER IFLD, INTYP(*), IFIELD(*)
      CHARACTER*(*) CFIELD(*)
      INTEGER BLKCOL(0:NELBLK)

      include 'params.blk'
      include 'dbnums.blk'
      include 'bcolr.blk'

      LOGICAL FFMATC, FFEXST
      INTEGER IRNG(3)
      LOGICAL NUMSPC
      INTEGER IDCOL, LOCSTR

      CHARACTER*(MXSTLN) COLA, COLF
      include 'cmap-lst.blk'

C *****************************************************************

      BCOLCH = .TRUE.

      IF (INIT) THEN

         BLKCOL(0) = -1
         DO 100 I = 1, NELBLK
            BLKCOL(I) = -2
  100    CONTINUE

      ELSE

C        If there are no fields on the BLKCOL command line,
C        toggle the ON/OFF flag.
         IF (INTYP(IFLD) .LE. -1) THEN
            BLKCOL(0) = -BLKCOL(0)

C        Check for ON flag.
         ELSE IF (FFMATC (IFLD, INTYP, CFIELD, 'ON', 2)) THEN
            BLKCOL(0) = 1
            CALL FFADDC ('ON', INLINE(1))

C        Check for OFF flag.
         ELSE IF (FFMATC (IFLD, INTYP, CFIELD, 'OFF', 2)) THEN
            BLKCOL(0) = -1
            CALL FFADDC ('OFF', INLINE(1))

C        Check for RESET flag.
         ELSE IF (FFMATC (IFLD, INTYP, CFIELD, 'RESET', 2)) THEN
            CALL FFADDC ('RESET', INLINE(1))
            BLKCOL(0) = -1
            DO 110 I = 1, NELBLK
               BLKCOL(I) = -2
  110       CONTINUE

         ELSE

  120       CONTINUE
            IF (FFEXST (IFLD, INTYP)) THEN

C                 Get numeric range identifying blocks whose colors
C                 will be specified.
               NUMSPC = .FALSE.
  130          CONTINUE

               CALL FFNRNG (IFLD, INTYP, CFIELD, IFIELD,
     &            'block id range', NELBLK, IRNG, *150, *170)
               CALL FFADDV (IRNG, INLINE)
               NUMSPC = .TRUE.

C                    Flag these blocks for color assignments.
               DO 140 I = IRNG(1), IRNG(2), IRNG(3)
                  IF (BLKCOL(I) .LT. 7) BLKCOL(I) = BLKCOL(I) + 10
  140          CONTINUE

               GO TO 130

  150          CONTINUE

C                 Get color to assign to the specified blocks.

C                    Check that next field has characters in it.

               IF (INTYP(IFLD) .GE. 0) THEN
                  COLA = CFIELD(IFLD)
                  IFLD = IFLD + 1
                  CALL ABRSTR (COLF, COLA, COLLST)
                  IF (COLF .EQ. ' ') THEN
                     WRITE (*, 10000) COLA
10000                 FORMAT (1X, A, ' not a valid color name.',
     &                  'rest of command not processed.')
                     GO TO 170
                  ELSE
                     IDCOL = LOCSTR (COLF, NCOLOR, COLLST)
                     IF (.NOT. NUMSPC) THEN
                        WRITE (*,10010) COLLST(IDCOL)
10010                    FORMAT (1X, 'No block ids were specified',
     &                     'for the color ', A)
                     ENDIF
                     CALL FFADDC (COLF, INLINE)
                     IF (IDCOL .GT. 0) THEN
                        IDCOL = IDCOL - 2
                        DO 160 I = 1, NELBLK
                           IF (BLKCOL(I) .GE. 7) BLKCOL(I) = IDCOL
  160                   CONTINUE
                     ELSE
                        WRITE (*, 10000) COLA
                        GO TO 170
                     ENDIF
                  ENDIF
               ELSE
                  WRITE (*, 10020)
10020              FORMAT (1X, 'No color specified following block',
     &               'id specifications')
                  GO TO 170
               ENDIF

               GO TO 120
            ENDIF
            BLKCOL(0) = 1

         ENDIF

      ENDIF
      GO TO 190

  170 CONTINUE
      BLKCOL(0) = -1
      DO 180 I = 1, NELBLK
         IF (BLKCOL(I) .GE. 7) BLKCOL(I) = BLKCOL(I) - 10
  180 CONTINUE

  190 CONTINUE
      RETURN
      END
