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

C=======================================================================
      SUBROUTINE CMDELB (VERB, INLINE, IFLD, INTYP, CFIELD, IFIELD,
     &                   IDELB, IELBST, NEWELB, *)
C=======================================================================

C   --*** CMDELB *** (MESH) Process active element block commands
C   --   Written by Amy Gilkey - revised 05/18/88
C   --
C   --Parameters:
C   --   VERB    - I/O - the verbs for the SHOW command
C   --   INLINE  - I/O - the parsed input line for the log file
C   --   IFLD    - I/O - the free-field reader index and fields
C   --   INTYP   - I/O - the free-field reader index and fields
C   --   CFIELD  - I/O - the free-field reader index and fields
C   --   IFIELD  - I/O - the free-field reader index and fields
C   --   IDELB   - IN  - the element block IDs
C   --   IELBST  - I/O - the element block status:
C   --                   -1 = OFF, 0 = ON, but not selected, 1 = selected
C   --   NEWELB  - I/O - the new element blocks flag, set 0 in ROUTER only
C   --                   0 = no new element blocks
C   --                   1 = new selected element blocks
C   --                   2 = new displayed element blocks
C   --                       (implies new selected blocks)
C   --
C   --Common Variables:
C   --   Uses IS3DIM of /D3NUMS/
C   --   Uses NELBLK of /DBNUMS/

      include 'params.blk'
      include 'dbnums.blk'
      include 'd3nums.blk'

      CHARACTER*(*) VERB
      CHARACTER*(*) INLINE
      INTEGER     INTYP(*)
      CHARACTER*(*) CFIELD(*)
      INTEGER     IFIELD(*)
      INTEGER IDELB(NELBLK)
      INTEGER IELBST(NELBLK)
      INTEGER NEWELB

      CHARACTER*(MXSTLN) WORD
      CHARACTER OPTION
      LOGICAL FFEXST, FFNUMB, FFMATC

      IF (VERB .EQ. 'VISIBLE') THEN
         CALL FFADDC (VERB, INLINE)
         IF (.NOT. IS3DIM) THEN
            CALL PRTERR ('CMDERR', 'Command allowed in 3D only')
            GOTO 190
         END IF

         NEWELB = 2

         IF (.NOT. FFEXST (IFLD, INTYP)) THEN

C         --Select all element blocks if no parameters
            DO 100 IELB = 1, NELBLK
               IF (IELBST(IELB) .LT. 0) IELBST(IELB) = 1
  100       CONTINUE

         ELSE IF (FFMATC (IFLD, INTYP, CFIELD, 'OFF', 3)) THEN

C         --Turn off all element blocks if OFF
            CALL FFADDC ('OFF', INLINE)
            DO 110 IELB = 1, NELBLK
               IF (IELBST(IELB) .GE. 0) IELBST(IELB) = -1
  110       CONTINUE

         ELSE

C         --Strip off ADD or DELETE option
            IF (FFMATC (IFLD, INTYP, CFIELD, 'ADD', 3)) THEN
               CALL FFADDC ('ADD', INLINE)
               OPTION = '+'
            ELSE IF (FFMATC (IFLD, INTYP, CFIELD, 'DELETE', 3)) THEN
               CALL FFADDC ('DELETE', INLINE)
               OPTION = '-'
            ELSE
               IF (.NOT. FFNUMB (IFLD, INTYP)) THEN
                  CALL PRTERR ('CMDERR',
     &               'Expected "OFF" or "ADD" or "DELETE"'
     &               // ' or element block ID')
                  GOTO 190
               END IF
C            --De-select all element blocks so only listed blocks are selected
               CALL INIINT (NELBLK, -1, IELBST)
               OPTION = '+'
            END IF

  120       CONTINUE
            IF (FFEXST (IFLD, INTYP)) THEN
               CALL FFINTG (IFLD, INTYP, IFIELD,
     &            'element block ID', 0, ID, *130)
               IELB = LOCINT (ID, NELBLK, IDELB)
               IF (IELB .GT. 0) THEN
                  CALL FFADDI (ID, INLINE)
                  IF (OPTION .EQ. '+') THEN
                     IF (IELBST(IELB) .LT. 0) IELBST(IELB) = 1
                  ELSE
                     IELBST(IELB) = -1
                  END IF
               ELSE
                  CALL INTSTR (1, 0, ID, WORD, LSTR)
                  CALL PRTERR ('CMDWARN', 'Element block ID '
     &               // WORD(:LSTR) // ' does not exist, ignored')
               END IF
  130          CONTINUE
               GOTO 120
            END IF

            CALL CNTELB (IELBST, NELBLK, NUMON, NUMSEL)
            IF (NUMON .EQ. 0) CALL FFADDC ('OFF', INLINE)
         END IF

      ELSE IF ((VERB .EQ. 'BLOCKS') .OR. (VERB .EQ. 'MATERIAL')) THEN
         CALL FFADDC (VERB, INLINE)
         NEWELB = MAX (NEWELB, 1)

         IF (.NOT. FFEXST (IFLD, INTYP)) THEN

C         --Select all element blocks if no parameters
            DO 140 IELB = 1, NELBLK
               IF (IELBST(IELB) .GE. 0) IELBST(IELB) = 1
  140       CONTINUE

         ELSE IF (FFMATC (IFLD, INTYP, CFIELD, 'OFF', 3)) THEN

C         --De-select all element blocks if OFF
            CALL FFADDC ('OFF', INLINE)
            DO 150 IELB = 1, NELBLK
               IF (IELBST(IELB) .GE. 0) IELBST(IELB) = 0
  150       CONTINUE

         ELSE

C         --Strip off ADD or DELETE option
            IF (FFMATC (IFLD, INTYP, CFIELD, 'ADD', 3)) THEN
               CALL FFADDC ('ADD', INLINE)
               OPTION = '+'
            ELSE IF (FFMATC (IFLD, INTYP, CFIELD, 'DELETE', 3)) THEN
               CALL FFADDC ('DELETE', INLINE)
               OPTION = '-'
            ELSE
               IF (.NOT. FFNUMB (IFLD, INTYP)) THEN
                  CALL PRTERR ('CMDERR',
     &               'Expected "OFF" or "ADD" or "DELETE"'
     &               // ' or element block ID')
                  GOTO 190
               END IF
C            --De-select all element blocks so only listed blocks are selected
               DO 160 IELB = 1, NELBLK
                  IF (IELBST(IELB) .GE. 0) IELBST(IELB) = 0
  160          CONTINUE
               OPTION = '+'
            END IF

  170       CONTINUE
            IF (FFEXST (IFLD, INTYP)) THEN
               CALL FFINTG (IFLD, INTYP, IFIELD,
     &            'element block ID', 0, ID, *180)
               IELB = LOCINT (ID, NELBLK, IDELB)
               IF (IELB .GT. 0) THEN
                  IF (IELBST(IELB) .GE. 0) THEN
                     CALL FFADDI (ID, INLINE)
                     IF (OPTION .EQ. '+') THEN
                        IELBST(IELB) = 1
                     ELSE
                        IELBST(IELB) = 0
                     END IF
                  ELSE
                     CALL INTSTR (1, 0, ID, WORD, LSTR)
                     CALL PRTERR ('CMDWARN', 'Element block '
     &                  // WORD(:LSTR) // ' is not displayed, ignored')
                  END IF
               ELSE
                  CALL INTSTR (1, 0, ID, WORD, LSTR)
                  CALL PRTERR ('CMDWARN', 'Element block ID '
     &               // WORD(:LSTR) // ' does not exist, ignored')
               END IF
  180          CONTINUE
               GOTO 170
            END IF

            CALL CNTELB (IELBST, NELBLK, NUMON, NUMSEL)
            IF (NUMSEL .EQ. 0) CALL FFADDC ('OFF', INLINE)
         END IF
      END IF

      RETURN

  190 CONTINUE
      RETURN 1
      END
