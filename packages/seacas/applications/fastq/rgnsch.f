C    Copyright (c) 2014, Sandia Corporation.
C    Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
C    the U.S. Government retains certain rights in this software.
C    
C    Redistribution and use in source and binary forms, with or without
C    modification, are permitted provided that the following conditions are
C    met:
C    
C        * Redistributions of source code must retain the above copyright
C          notice, this list of conditions and the following disclaimer.
C    
C        * Redistributions in binary form must reproduce the above
C          copyright notice, this list of conditions and the following
C          disclaimer in the documentation and/or other materials provided
C          with the distribution.
C    
C        * Neither the name of Sandia Corporation nor the names of its
C          contributors may be used to endorse or promote products derived
C          from this software without specific prior written permission.
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

C $Id: rgnsch.f,v 1.1 1990/11/30 11:14:59 gdsjaar Exp $
C $Log: rgnsch.f,v $
C Revision 1.1  1990/11/30 11:14:59  gdsjaar
C Initial revision
C
C
CC* FILE: [.QMESH]RGNSCH.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 8/1/90
CC* MODIFICATION: FORCED SCHEME TO BE "X" WHEN REMESHING FOR
CC*               ERROR ESTIMATION - THIS INCLUDED CHANGING THE
CC*               CALL TO ADD THE LOGICAL REMESH.
C
      SUBROUTINE RGNSCH (MCOM, ICOM, JCOM, CIN, IIN, RIN, KIN, STEP,
     &   IREGN, IPNTR, N24, MSC, SCHEME, DEFSCH, SCHSTR, LENSCH, NPER,
     &   PENTAG, TRIANG, TRNSIT, HALFC, FILL, ICODE, REMESH)
C***********************************************************************
C
C     RGNSCH - GET A REGION'S SCHEME
C
C***********************************************************************
C
      DIMENSION CIN(MCOM), IIN(MCOM), RIN(MCOM), KIN(MCOM)
      DIMENSION SCHEME(MSC)
C
      CHARACTER*72 SCHEME, DEFSCH, SCHSTR, CIN
C
      LOGICAL STEP, PENTAG, TRIANG, TRNSIT, HALFC, FILL, IANS, REMESH
C
      DATA IEXIT, IOVER, IQUIT /1, 2, 3/
C
      ICODE = 0
C
C  CHECK FOR REMESHING
C
      IF (REMESH) THEN
         SCHSTR = 'X'
      ELSE
C
C  GET THE INITIAL SCHEME
C
         IF ((ABS(IREGN) .LE. N24) .AND. (IPNTR .GT. 0)) THEN
            SCHSTR = SCHEME(IPNTR)
         ELSE
            SCHSTR = DEFSCH
         END IF
      END IF
      CALL STRCUT (SCHSTR)
      CALL STRLNG (SCHSTR, LENSCH)
C
C  STEP PROCESSING
C
      IF (STEP) THEN
         WRITE (*, 10000) SCHSTR(1:LENSCH)
         CALL INTRUP ('USE CURRENT SCHEME TO BEGIN PROCESSING', IANS,
     &      MCOM, ICOM, JCOM, CIN, IIN, RIN, KIN)
C
C  CHANGE THE SCHEME
C
         IF (.NOT.IANS) THEN
  100       CONTINUE
            IF (ICOM .LE. JCOM) THEN
               SCHSTR = CIN(ICOM)
               ICOM = ICOM + 1
               IANS = .TRUE.
            ELSE
               CALL INQSTR ('ENTER PROCESSING SCHEME: ', SCHSTR)
            END IF
            CALL STRCUT (SCHSTR)
            CALL STRLNG (SCHSTR, LENSCH)
C
C  HELP FOR SCHEME
C
            IF ((SCHSTR(1:1) .EQ. 'H') .OR.
     &         (SCHSTR(1:1) .EQ. 'h')) THEN
               CALL MESAGE (' ')
               CALL HELP_FQ (13)
               CALL MESAGE (' ')
               GO TO 100
            END IF
         END IF
C
C  BLANK SCHEME
C
         IF ((LENSCH .LE. 0) .OR. (SCHSTR(1:1) .EQ. ' ')) THEN
            CALL MESAGE ('NO INITIAL SCHEME INPUT')
            CALL MESAGE ('FORCED RECTANGLE PROCESSING USED')
            SCHSTR = ' '
            LENSCH = 1
            GO TO 120
         END IF
      END IF
C
C  DETERMINE MESHING SCHEME
C
      PENTAG = .FALSE.
      TRIANG = .FALSE.
      TRNSIT = .FALSE.
      FILL = .FALSE.
      DO 110 J = 1, LENSCH
C
C  SEE IF A PENTAGON REGION HAS BEEN FLAGGED
C
         IF ((SCHSTR(J:J) .EQ. 'U') .OR. (SCHSTR(J:J) .EQ. 'u')) THEN
            IF (NPER .GE. 10) THEN
               PENTAG = .TRUE.
               CALL MESAGE
     &            ('PENTAGON PRIMITIVE REGION PROCESSING USED')
            ELSE
               CALL MESAGE ('PENTAGON REGION GENERATION NOT')
               CALL MESAGE ('POSSIBLE WITH NO. IN PERIMETER < 10')
               CALL MESAGE ('REGULAR PROCESSING WILL BE ATTEMPTED')
            END IF
            GO TO 120
C
C  SEE IF A TRANSITION REGION HAS BEEN FLAGGED
C
         ELSE IF ((SCHSTR(J:J) .EQ. 'B') .OR.
     &      (SCHSTR(J:J) .EQ. 'b')) THEN
            IF (NPER .GE. 8) THEN
               TRNSIT = .TRUE.
               HALFC = .FALSE.
               CALL MESAGE
     &            ('TRANSITION PRIMITIVE REGION PROCESSING USED')
            ELSE
               CALL MESAGE ('TRANSITION REGION GENERATION NOT')
               CALL MESAGE ('POSSIBLE WITH NO. IN PERIMETER < 8')
               CALL MESAGE ('REGULAR PROCESSING WILL BE ATTEMPTED')
            END IF
            GO TO 120
C
C  SEE IF A SEMI-CIRCLE REGION HAS BEEN FLAGGED
C
         ELSE IF ((SCHSTR(J:J) .EQ. 'C') .OR.
     &      (SCHSTR(J:J) .EQ. 'c')) THEN
            IF (NPER .GE. 8) THEN
               TRNSIT = .TRUE.
               HALFC = .TRUE.
               CALL MESAGE
     &            ('SEMICIRCLE PRIMITIVE REGION PROCESSING USED')
            ELSE
               CALL MESAGE
     &            ('TRANSITION/SEMICIRCLE REGION GENERATION NOT')
               CALL MESAGE ('POSSIBLE WITH NO. IN PERIMETER < 8')
               CALL MESAGE ('REGULAR PROCESSING WILL BE ATTEMPTED')
            END IF
            GO TO 120
C
C  SEE IF A TRIANGULAR REGION HAS BEEN FLAGGED
C
         ELSE IF ((SCHSTR(J:J) .EQ. 'T') .OR.
     &      (SCHSTR(J:J) .EQ. 't')) THEN
            IF (NPER .GE. 6) THEN
               TRIANG = .TRUE.
               CALL MESAGE
     &            ('TRIANGLE PRIMITIVE REGION PROCESSING USED')
            ELSE
               CALL MESAGE ('TRIANGULAR REGION GENERATION NOT')
               CALL MESAGE ('POSSIBLE WITH NO. IN PERIMETER < 6')
               CALL MESAGE ('REGULAR PROCESSING WILL BE ATTEMPTED')
            END IF
            GO TO 120
C
C  SEE IF A FILL REGION HAS BEEN FLAGGED
C
         ELSE IF ((SCHSTR(J:J) .EQ. 'X') .OR.
     &      (SCHSTR(J:J) .EQ. 'x')) THEN
            FILL = .TRUE.
            CALL MESAGE ('PAVING TECHNIQUE INITIALLY USED')
            GO TO 120
C
C  SEE IF A REGULAR RECTANGULAR REGION HAS BEEN FLAGGED
C
         ELSE IF ((SCHSTR(J:J) .EQ. 'M') .OR.
     &      (SCHSTR(J:J) .EQ. 'm')) THEN
            GO TO 120
C
C  OTHER POSSIBILITIES
C
         ELSE IF ((SCHSTR(J:J) .EQ. 'E') .OR.
     &      (SCHSTR(J:J) .EQ. 'e')) THEN
            ICODE = IEXIT
            GO TO 120
         ELSE IF ((SCHSTR(J:J) .EQ. 'O') .OR.
     &      (SCHSTR(J:J) .EQ. 'o')) THEN
            ICODE = IOVER
            GO TO 120
         ELSE IF ((SCHSTR(J:J) .EQ. 'Q') .OR.
     &      (SCHSTR(J:J) .EQ. 'q')) THEN
            ICODE = IQUIT
            GO TO 120
         END IF
  110 CONTINUE
  120 CONTINUE
C
      RETURN
C
10000 FORMAT ('0INITIAL MESH DEFINED USING THIS SCHEME:', /, 5X, A)
      END
