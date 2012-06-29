C    Copyright(C) 2008 Sandia Corporation.  Under the terms of Contract
C    DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
C    certain rights in this software
C    
C    Redistribution and use in source and binary forms, with or without
C    modification, are permitted provided that the following conditions are
C    met:
C    
C    * Redistributions of source code must retain the above copyright
C       notice, this list of conditions and the following disclaimer.
C              
C    * Redistributions in binary form must reproduce the above
C      copyright notice, this list of conditions and the following
C      disclaimer in the documentation and/or other materials provided
C      with the distribution.
C                            
C    * Neither the name of Sandia Corporation nor the names of its
C      contributors may be used to endorse or promote products derived
C      from this software without specific prior written permission.
C                                                    
C    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
C    "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
C    LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
C    A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
C    OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
C    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
C    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
C    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
C    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
C    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
C    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
C    

C=======================================================================
      SUBROUTINE APARSE (LINE, MAXFLD,NUMFLD, ITYPE, CFIELD, RFIELD)
C=======================================================================

C   --*** APARSE *** (ALGEBRA) Parse line into fields
C   --   Written by Amy Gilkey - revised 02/23/88
C   --
C   --APARSE converts a line to uppercase and parses the line into individual
C   --fields.  A field is one of the following:
C   --   a number (starts with 0-9 or ., rest = 0-9 or . or E or E+ or E-)
C   --   a name (starts with A-Z, rest = A-Z or $ or :)
C   --   a non-number or non-name character; ** is converted to ^
C   --
C   --Parameters:
C   --   LINE - IN - the input line to be parsed, all uppercase
C   --   MAXFLD - IN - the maximum number of parsed fields
C   --   NUMFLD - OUT - the number of parsed fields
C   --   ITYPE - OUT - the field types:
C   --      -1= no value
C   --      0 = CFIELD only valid, name or invalid number string
C   --      1 = CFIELD and RFIELD only valid
C   --      3 = CFIELD only valid, character string
C   --   CFIELD - OUT - the alphanumeric fields
C   --   RFIELD - OUT - the numeric fields

      CHARACTER*(*) LINE
      INTEGER       ITYPE(*)
      CHARACTER*(*) CFIELD(*)
      REAL          RFIELD(*)

      CHARACTER*32 CONVER
      CHARACTER CH

C   --Initialize fields

      CALL INIINT (MAXFLD, -1, ITYPE)
      CALL INISTR (MAXFLD, ' ', CFIELD)
      CALL INIREA (MAXFLD, 0.0, RFIELD)
      NUMFLD = 0

C   --Change line to upper case

      CALL EXUPCS (LINE)

C   --Repeat until done with input characters

      IEND = LENSTR (LINE)

      NCOL = 1
  100 CONTINUE
      IF (NCOL .LE. IEND) THEN

C      --Skip blank characters

  110    CONTINUE
         IF (LINE(NCOL:NCOL) .EQ. ' ') THEN
            NCOL = NCOL + 1
            GOTO 110
         END IF

         CH = LINE(NCOL:NCOL)

         NUMFLD = NUMFLD + 1

         IF (((CH .GE. '0') .AND. (CH .LE. '9'))
     &      .OR. (CH .EQ. '.')) THEN

C         --Get number string

            ILEFT = NCOL
            IRIGHT = 0
  120       CONTINUE
            IF (IRIGHT .EQ. 0) THEN
               NCOL = NCOL + 1
               CH = LINE(NCOL:NCOL)
               IF (CH .EQ. 'E') THEN
                  CH = LINE(NCOL+1:NCOL+1)
                  IF ((CH .EQ. '+') .OR. (CH .EQ. '-')) NCOL = NCOL + 1
               ELSE IF (CH .EQ. '.') THEN
                  CONTINUE
               ELSE IF ((CH .GE. '0') .AND. (CH .LE. '9')) THEN
                  CONTINUE
               ELSE
                  IRIGHT = NCOL - 1
               END IF
               GOTO 120
            END IF

C         --Convert number string and store

            IF (NUMFLD .LE. MAXFLD) THEN
               CFIELD(NUMFLD) = LINE(ILEFT:IRIGHT)
               ITYPE(NUMFLD) = 0

               CONVER = ' '
               IJUST = 32 - (IRIGHT - ILEFT + 1) + 1
               CONVER(IJUST:32) = LINE(ILEFT:IRIGHT)

               READ (CONVER, '(E32.0)', IOSTAT=ITRANS) RNUM
               IF (ITRANS .EQ. 0) THEN
                  RFIELD(NUMFLD) = RNUM
                  ITYPE(NUMFLD) = 1
               END IF
            END IF

         ELSE IF ((CH .GE. 'A') .AND. (CH .LE. 'Z')) THEN

C         --Get word string and store

            ILEFT = NCOL
            IRIGHT = 0
  130       CONTINUE
            IF (IRIGHT .EQ. 0) THEN
               NCOL = NCOL + 1
               CH = LINE(NCOL:NCOL)
               IF (((CH .GE. 'A') .AND. (CH .LE. 'Z'))
     &            .OR. ((CH .GE. '0') .AND. (CH .LE. '9'))
     &            .OR. (CH .EQ. '$') .OR. (CH .EQ. ':')
     &            .OR. (CH .EQ. '_') .OR. (CH .EQ. '.')) THEN
                  CONTINUE
               ELSE
                  IRIGHT = NCOL - 1
               END IF
               GOTO 130
            END IF

            IF (NUMFLD .LE. MAXFLD) THEN
               CFIELD(NUMFLD) = LINE(ILEFT:IRIGHT)
               ITYPE(NUMFLD) = 0
            END IF

         ELSE

C         --Get single character (or **) and store
C         --Includes { , ( ) = }

            NCOL = NCOL + 1
            IF ((CH .EQ. '*') .AND. (LINE(NCOL:NCOL) .EQ. '*')) THEN
               CH = '^'
               NCOL = NCOL + 1
            END IF

            IF (NUMFLD .LE. MAXFLD) THEN
               CFIELD(NUMFLD) = CH
               ITYPE(NUMFLD) = 3
            END IF
         END IF

         GOTO 100
      END IF

      RETURN
      END
