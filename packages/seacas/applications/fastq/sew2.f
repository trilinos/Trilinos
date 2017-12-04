C    Copyright(C) 2014-2017 National Technology & Engineering Solutions of
C    Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
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
C    * Neither the name of NTESS nor the names of its
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

C $Id: sew2.f,v 1.1 1990/11/30 11:15:32 gdsjaar Exp $
C $Log: sew2.f,v $
C Revision 1.1  1990/11/30 11:15:32  gdsjaar
C Initial revision
C
C
CC* FILE: [.PAVING]SEW2.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE SEW2 (MXND, MLN, NUID, LXK, KXL, NXL, LXN, LNODES,
     &   IAVAIL, NAVAIL, LLL, KKK, NNN, I1, I2, J1, J2, NOROOM, ERR)
C***********************************************************************
C
C  SUBROUTINE SEW2 = COLLAPSES A LOOP INTO TWO POSSIBLE LOOPS
C
C***********************************************************************
C
      DIMENSION NUID (MXND)
      DIMENSION LXK (4, MXND), KXL (2, 3*MXND)
      DIMENSION NXL (2, 3*MXND), LXN (4, MXND)
      DIMENSION LNODES (MLN, MXND)
      DIMENSION L1LIST(20)
C
      LOGICAL ERR, NOROOM
C
      ERR = .FALSE.
C
C  GET THE APPROPRIATE LINES AND NODES TO BE DELETED
C
      IF ((LXN (2, J1) .LT. 0) .AND. (LXN (2, J2) .LT. 0)) THEN
         LSTAY = LNODES (5, J1)
         LGONE = LNODES (5, I1)
         NGONE1 = I1
         NGONE2 = I2
         NSTAY1 = J2
         NSTAY2 = J1
C
      ELSEIF (LXN (2, J1) .LT. 0) THEN
         LSTAY = LNODES (5, I1)
         LGONE = LNODES (5, J1)
         NGONE1 = J2
         NGONE2 = I2
         NSTAY1 = I1
         NSTAY2 = J1
C
      ELSEIF (LXN (2, J2) .LT. 0) THEN
         LSTAY = LNODES (5, I1)
         LGONE = LNODES (5, J1)
         NGONE1 = I1
         NGONE2 = J1
         NSTAY1 = J2
         NSTAY2 = I2
C
      ELSE
         LSTAY = LNODES (5, I1)
         LGONE = LNODES (5, J1)
         NGONE1 = J2
         NGONE2 = J1
         NSTAY1 = I1
         NSTAY2 = I2
      ENDIF
C
      KOLD = KXL (1, LGONE)
      KNEW = KXL (1, LSTAY)
C
C  DELETE THE OLD LINE AND REDO LINK ARRAYS
C
      IF (KNEW .EQ. 0) THEN
         KXL (1, LSTAY) = KOLD
         KXL (2, LSTAY) = 0
      ELSE
         KXL (1, LSTAY) = KNEW
         KXL (2, LSTAY) = KOLD
      ENDIF
C
      KXL (1, LGONE) = 0
      KXL (2, LGONE) = 0
      NXL (1, LGONE) = 0
      NXL (2, LGONE) = 0
C
C  FIX THE LINES PER ELEMENT ARRAY FOR THE ONE ELEMENT CHANGING
C
      IF (KOLD .GT. 0) THEN
         DO 100 II = 1, 4
            IF (LXK (II, KOLD) .EQ. LGONE) THEN
               LXK (II, KOLD) = LSTAY
               GOTO 110
            ENDIF
  100    CONTINUE
C
         CALL MESAGE ('** PROBLEMS IN SEW2 FIXING THE CHANGING'//
     &      'ELEMENT **')
         ERR = .TRUE.
         GOTO 180
C
  110    CONTINUE
      ENDIF
C
C  RECONNECT ALL LINES CONNECTING TO NGONE2 TO NSTAY2
C
      CALL GETLXN (MXND, LXN, NGONE2, L1LIST, NL, ERR)
      IF (ERR) THEN
         CALL MESAGE ('** PROBLEMS IN SEW2 FINDING LXN FOR NGONE2 **')
         GOTO 180
      ENDIF
      DO 120 II = 1, NL
         LL = L1LIST (II)
         IF (NXL (1, LL) .EQ. NGONE2) THEN
            NXL (1, LL) = NSTAY2
         ELSEIF (NXL (2, LL) .EQ. NGONE2) THEN
            NXL (2, LL) = NSTAY2
         ENDIF
  120 CONTINUE
C
C  FIX LXN ARRAY
C  UNHOOK LGONE FROM NGONE2 OR NSTAY2 AS NEEDED
C
      IF (LGONE .EQ. LNODES (5, J1)) THEN
         CALL DELLXN (MXND, LXN, NUID, NAVAIL, IAVAIL, J1,
     &      LGONE, NNN, ERR, NOROOM)
         IF ((NOROOM) .OR. (ERR)) THEN
            CALL MESAGE ('** PROBLEMS IN SEW2 DELETING NGONE2 LINES **')
            GOTO 180
         ENDIF
C
      ELSE
         CALL DELLXN (MXND, LXN, NUID, NAVAIL, IAVAIL, I2,
     &      LGONE, NNN, ERR, NOROOM)
         IF ((NOROOM) .OR. (ERR)) THEN
            CALL MESAGE ('** PROBLEMS IN SEW2 DELETING NGONE2 LINES **')
            GOTO 180
         ENDIF
C
      ENDIF
C
C  ADD ALL LINES STILL HOOKED TO NGONE2 TO THE LIST OF LINES FOR NSTAY2
C
      DO 130 II = 1, NL
         LL = L1LIST (II)
         IF (LL .NE. LGONE) THEN
            CALL ADDLXN (MXND, LXN, NUID, NAVAIL, IAVAIL,
     &         NSTAY2, LL, NNN, ERR, NOROOM)
            IF ((NOROOM) .OR. (ERR)) THEN
               CALL MESAGE ('** PROBLEMS IN SEW2 ADDING NSTAY2 '//
     &            'LINES **')
               GOTO 180
            ENDIF
         ENDIF
  130 CONTINUE
C
C  DELETE NGONE2 (UNHOOK EVERYTHING FROM IT)
C
      DO 140 II = 1, 3
         LXN (II, NGONE2) = 0
  140 CONTINUE
      LXN (4, NGONE2) = IAVAIL
      IAVAIL = NGONE2
      NAVAIL = NAVAIL+1
      NUID (NGONE2) = 0
C
C  RECONNECT ALL LINES CONNECTING TO NGONE1 TO NSTAY1
C
      CALL GETLXN (MXND, LXN, NGONE1, L1LIST, NL, ERR)
      IF (ERR) THEN
         CALL MESAGE ('** PROBLEMS IN SEW2 GETTING NGONE1 LINES **')
         GOTO 180
      ENDIF
      DO 150 II = 1, NL
         LL = L1LIST (II)
         IF (NXL (1, LL) .EQ. NGONE1) THEN
            NXL (1, LL) = NSTAY1
         ELSEIF (NXL (2, LL) .EQ. NGONE1) THEN
            NXL (2, LL) = NSTAY1
         ENDIF
  150 CONTINUE
C
C  FIX LXN ARRAY
C  UNHOOK LGONE FROM NGONE1 OR NSTAY1 AS APPROPRIATE
C
      IF (LGONE .EQ. LNODES (5, I1)) THEN
         CALL DELLXN (MXND, LXN, NUID, NAVAIL, IAVAIL, I1,
     &      LGONE, NNN, ERR, NOROOM)
         IF ((NOROOM) .OR. (ERR)) THEN
            CALL MESAGE ('** PROBLEMS IN SEW2 DELETING NGONE1 LINES **')
            GOTO 180
         ENDIF
C
      ELSE
         CALL DELLXN (MXND, LXN, NUID, NAVAIL, IAVAIL, J2,
     &      LGONE, NNN, ERR, NOROOM)
         IF ((NOROOM) .OR. (ERR)) THEN
            CALL MESAGE ('** PROBLEMS IN SEW2 DELETING NGONE1 LINES **')
            GOTO 180
         ENDIF
C
      ENDIF
C
C  ADD ALL LINES STILL HOOKED TO NGONE1 TO THE LIST OF LINES FOR NSTAY1
C
      DO 160 II = 1, NL
         LL = L1LIST (II)
         IF (LL .NE. LGONE) THEN
            CALL ADDLXN (MXND, LXN, NUID, NAVAIL, IAVAIL,
     &         NSTAY1, LL, NNN, ERR, NOROOM)
            IF ((NOROOM) .OR. (ERR)) THEN
               CALL MESAGE ('** PROBLEMS IN SEW2 ADDING NSTAY1'//
     &            ' LINES **')
               GOTO 180
            ENDIF
         ENDIF
  160 CONTINUE
C
C  DELETE NGONE1 (UNHOOK EVERYTHING FROM IT)
C
      DO 170 II = 1, 3
         LXN (II, NGONE1) = 0
  170 CONTINUE
      LXN (4, NGONE1) = IAVAIL
      IAVAIL = NGONE1
      NAVAIL = NAVAIL+1
      NUID (NGONE1) = 0
C
C  NOW FIX THE LNODES ARRAY FOR BOTH OF THE LOOPS
C
      IF (NGONE2 .EQ. J1) THEN
         LNODES (2, NSTAY2) = LNODES (2, NGONE2)
         LNODES (3, LNODES (2, NGONE2)) = NSTAY2
      ELSE
         LNODES (2, LNODES (3, NGONE2)) = NSTAY2
         LNODES (3, NSTAY2) = LNODES (3, NGONE2)
         LNODES (5, NSTAY2) = LNODES (5, NGONE2)
      ENDIF
      IF (NGONE1 .EQ. J2) THEN
         LNODES (2, LNODES (3, NGONE1)) = NSTAY1
         LNODES (3, NSTAY1) = LNODES (3, NGONE1)
         LNODES (5, NSTAY1) = LNODES (5, NGONE1)
      ELSE
         LNODES (2, NSTAY1) = LNODES (2, NGONE1)
         LNODES (3, LNODES (2, NGONE1)) = NSTAY1
      ENDIF
C
      I1 = NSTAY1
      I2 = NSTAY2
      J1 = NGONE1
      J2 = NGONE2
      CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES,
     &   I1, ERR)
      IF (ERR) GOTO 180
      CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES,
     &   LNODES (3, I1), ERR)
      IF (ERR) GOTO 180
      CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES,
     &   LNODES (3, LNODES (3, I1)), ERR)
      IF (ERR) GOTO 180
      CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES,
     &   LNODES (2, I1), ERR)
      IF (ERR) GOTO 180
      CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES,
     &   LNODES (2, LNODES (2, I1)), ERR)
      IF (ERR) GOTO 180
      CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES,
     &   I2, ERR)
      IF (ERR) GOTO 180
      CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES,
     &   LNODES (3, I2), ERR)
      IF (ERR) GOTO 180
      CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES,
     &   LNODES (3, LNODES (3, I2)), ERR)
      IF (ERR) GOTO 180
      CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES,
     &   LNODES (2, I2), ERR)
      IF (ERR) GOTO 180
      CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES,
     &   LNODES (2, LNODES (2, I2)), ERR)
      IF (ERR) GOTO 180
C
  180 CONTINUE
C
      RETURN
C
      END
