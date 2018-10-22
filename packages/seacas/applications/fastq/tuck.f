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

C $Id: tuck.f,v 1.1 1990/11/30 11:17:21 gdsjaar Exp $
C $Log: tuck.f,v $
C Revision 1.1  1990/11/30 11:17:21  gdsjaar
C Initial revision
C
C
CC* FILE: [.PAVING]TUCK.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE TUCK (MXND, MLN, NUID, XN, YN, LXK, KXL, NXL, LXN,
     &   LNODES, IAVAIL, NAVAIL, LLL, KKK, NNN, N1, NLOOP, GRAPH,
     &   NOROOM, ERR)
C***********************************************************************
C
C  SUBROUTINE TUCK = COLLAPSES TWO SIDE LINES INTO A ROW END NODE.
C                      THIS IS REFERRED TO AS A TUCK.
C
C***********************************************************************
C
      DIMENSION XN (MXND), YN (MXND), NUID (MXND)
      DIMENSION LXK (4, MXND), KXL (2, 3*MXND)
      DIMENSION NXL (2, 3*MXND), LXN (4, MXND)
      DIMENSION LNODES (MLN, MXND)
      DIMENSION L1LIST(20)
C
      LOGICAL GRAPH, ERR, NOROOM
C
      ERR = .FALSE.
C
C  CHECK TO MAKE SURE THAT THE NODE STILL EXISTS
C
      IF (LXN (1, N1) .LE. 0) THEN
         ERR = .TRUE.
         CALL MESAGE ('** PROBLEMS IN TUCK - N1 DOES NOT EXIST **')
         GOTO 290
      ENDIF
C
C  GET ALL THE DEFINITIONS IN ORDER
C
      N0 = LNODES (2, N1)
      N2 = LNODES (3, N1)
      L1 = LNODES (5, N0)
      L2 = LNODES (5, N1)
      KOLD = KXL (1, L1)
      KL2 = KXL (1, L2)
C
C  FIND L5 AND NC2
C
      DO 100 I = 1, 4
         LTEST = LXK (I, KOLD)
         IF (LTEST .NE. L1) THEN
            IF (NXL (1, LTEST) .EQ. N1) THEN
               L5 = LTEST
               NC2 = NXL (2, LTEST)
               GOTO 110
            ELSEIF (NXL (2, LTEST) .EQ. N1) THEN
               L5 = LTEST
               NC2 = NXL (1, LTEST)
               GOTO 110
            ENDIF
         ENDIF
  100 CONTINUE
      CALL MESAGE ('** PROBLEMS IN TUCK FINDING L5 AND NC2 **')
      ERR = .TRUE.
      GOTO 290
  110 CONTINUE
C
C  FIND L4 AND NC1
C
      DO 120 I = 1, 4
         LTEST = LXK (I, KOLD)
         IF ( (LTEST .NE. L1) .AND. (LTEST .NE. L5) ) THEN
            IF (NXL (1, LTEST) .EQ. N0) THEN
               L4 = LTEST
               NC1 = NXL (2, LTEST)
               GOTO 130
            ELSEIF (NXL (2, LTEST) .EQ. N0) THEN
               L4 = LTEST
               NC1 = NXL (1, LTEST)
               GOTO 130
            ENDIF
         ENDIF
  120 CONTINUE
      CALL MESAGE ('** PROBLEMS IN TUCK FINDING L4 AND NC1 **')
      ERR = .TRUE.
      GOTO 290
  130 CONTINUE
C
C  FIND L3
C
      DO 140 I = 1, 4
         LTEST = LXK (I, KOLD)
         IF ( (LTEST .NE. L1) .AND. (LTEST .NE. L5) .AND.
     &      (LTEST .NE. L4) ) THEN
            L3 = LTEST
            GOTO 150
         ENDIF
  140 CONTINUE
      CALL MESAGE ('** PROBLEMS IN TUCK FINDING L3 **')
      ERR = .TRUE.
      GOTO 290
  150 CONTINUE
C
C  FIND THE ELEMENT KL5
C
      IF (KXL (1, L5) .EQ. KOLD) THEN
         KL5 = KXL (2, L5)
      ELSEIF (KXL (2, L5) .EQ. KOLD) THEN
         KL5 = KXL (1, L5)
      ELSE
         CALL MESAGE ('** PROBLEMS IN TUCK FINDING KL5 **')
         ERR = .TRUE.
         GOTO 290
      ENDIF
C
C  NOW THAT ALL THE NECESSARY VARAIBLES HAVE BEEN DEFINED,
C  START BY DELETING LINE L1, L2, AND L5
C
      KXL (1, L1) = 0
      KXL (2, L1) = 0
      NXL (1, L1) = 0
      NXL (2, L1) = 0
      KXL (1, L2) = 0
      KXL (2, L2) = 0
      NXL (1, L2) = 0
      NXL (2, L2) = 0
      KXL (1, L5) = 0
      KXL (2, L5) = 0
      NXL (1, L5) = 0
      NXL (2, L5) = 0
C
C  NOW FIX THE KXL ARRAY FOR LINE L3 HAVING KL5 INSTEAD OF KOLD
C
      IF (KXL (1, L3) .EQ. KOLD) THEN
         KXL (1, L3) = KL5
      ELSEIF (KXL (2, L3) .EQ. KOLD) THEN
         KXL (2, L3) = KL5
      ELSE
         CALL MESAGE ('** PROBLEMS IN TUCK REPLACING KOLD FOR L3 **')
         ERR = .TRUE.
         GOTO 290
      ENDIF
C
C  NOW FIX THE KXL ARRAY FOR LINE L3 HAVING KL5 INSTEAD OF KOLD
C
      IF (KXL (1, L4) .EQ. KOLD) THEN
         KXL (1, L4) = KL2
      ELSEIF (KXL (2, L4) .EQ. KOLD) THEN
         KXL (2, L4) = KL2
      ELSE
         CALL MESAGE ('** PROBLEMS IN TUCK REPLACING KOLD FOR L4 **')
         ERR = .TRUE.
         GOTO 290
      ENDIF
C
C  FIX THE LINES PER ELEMENT ARRAY FOR ELEMENT KL5 TO REFLECT
C  L3 TAKING L5'S PLACE
C
      DO 160 I = 1, 4
         IF (LXK (I, KL5) .EQ. L5) THEN
            LXK (I, KL5) = L3
            GOTO 170
         ENDIF
  160 CONTINUE
      CALL MESAGE ('** PROBLEMS IN TUCK FINDING L5 IN KL5 **')
      ERR = .TRUE.
      GOTO 290
  170 CONTINUE
C
C  FIX THE LINES PER ELEMENT ARRAY FOR ELEMENT KL2 TO REFLECT
C  L4 TAKING L2'S PLACE
C
      DO 180 I = 1, 4
         IF (LXK (I, KL2) .EQ. L2) THEN
            LXK (I, KL2) = L4
            GOTO 190
         ENDIF
  180 CONTINUE
      CALL MESAGE ('** PROBLEMS IN TUCK FINDING L2 IN KL2 **')
      ERR = .TRUE.
      GOTO 290
  190 CONTINUE
C
C  RECONNECT ALL LINES CONNECTED TO N1 TO NC1 EXCEPT L5 AND L2
C
      CALL GETLXN (MXND, LXN, N1, L1LIST, NL, ERR)
      IF (ERR) THEN
         CALL MESAGE ('** PROBLEMS IN TUCK GETTING N1 LINES **')
         GOTO 290
      ENDIF
      IF (GRAPH) CALL LCOLOR ('BLACK')
      DO 200 I = 1, NL
         LL = L1LIST (I)
         IF ((GRAPH) .AND. (NXL (1, LL) .GT. 0) .AND.
     &      (NXL (2, LL) .GT. 0) )
     &      CALL D2NODE (MXND, XN, YN, NXL (1, LL), NXL (2, LL))
         IF (NXL (1, LL) .EQ. N1) THEN
            NXL (1, LL) = NC1
         ELSEIF (NXL (2, LL) .EQ. N1) THEN
            NXL (2, LL) = NC1
         ENDIF
  200 CONTINUE
      IF (GRAPH) THEN
         CALL LCOLOR ('WHITE')
         CALL SFLUSH
      ENDIF
C
C  FIX LXN ARRAY
C  UNHOOK L1, L2 AND L5 FROM N1
C
      CALL DELLXN (MXND, LXN, NUID, NAVAIL, IAVAIL, N1,
     &   L1, NNN, ERR, NOROOM)
      IF ((NOROOM) .OR. (ERR)) THEN
         CALL MESAGE ('** PROBLEMS IN TUCK UNHOOKING L1 FROM N1 **')
         GOTO 290
      ENDIF
      CALL DELLXN (MXND, LXN, NUID, NAVAIL, IAVAIL, N1,
     &   L2, NNN, ERR, NOROOM)
      IF ((NOROOM) .OR. (ERR)) THEN
         CALL MESAGE ('** PROBLEMS IN TUCK UNHOOKING L2 FROM N1 **')
         GOTO 290
      ENDIF
      CALL DELLXN (MXND, LXN, NUID, NAVAIL, IAVAIL, N1,
     &   L5, NNN, ERR, NOROOM)
      IF ((NOROOM) .OR. (ERR)) THEN
         CALL MESAGE ('** PROBLEMS IN TUCK UNHOOKING L5 FROM N1 **')
         GOTO 290
      ENDIF
C
C  ADD ALL LINES STILL HOOKED TO N1 TO THE LIST OF LINES FOR NC1
C
      DO 210 I = 1, NL
         LL = L1LIST (I)
         IF ((LL .NE. L2) .AND. (LL .NE. L5) .AND. (LL .NE. L1)) THEN
            CALL ADDLXN (MXND, LXN, NUID, NAVAIL, IAVAIL,
     &         NC1, LL, NNN, ERR, NOROOM)
            IF ((NOROOM) .OR. (ERR)) THEN
               CALL MESAGE ('** PROBLEMS IN TUCK HOOKING N1'' LINES'//
     &            ' TO NC1 **')
               GOTO 290
            ENDIF
         ENDIF
  210 CONTINUE
C
C  DELETE N1
C
      DO 220 I = 1, 3
         LXN (I, N1) = 0
  220 CONTINUE
      LXN (4, N1) = IAVAIL
      IAVAIL = N1
      NAVAIL = NAVAIL+1
      NUID (N1) = 0
C
C  RECONNECT ALL LINES CONNECTED TO N2 TO N0 (EXCEPT L2)
C
      CALL GETLXN (MXND, LXN, N2, L1LIST, NL, ERR)
      IF (ERR) THEN
         CALL MESAGE ('** PROBLEMS IN TUCK GETTING N2 LINES **')
         GOTO 290
      ENDIF
      IF (GRAPH) CALL LCOLOR ('BLACK')
      DO 230 I = 1, NL
         LL = L1LIST (I)
         IF ((GRAPH) .AND. (NXL (1, LL) .GT. 0) .AND.
     &      (NXL (2, LL) .GT. 0) )
     &      CALL D2NODE (MXND, XN, YN, NXL (1, LL), NXL (2, LL))
         IF (NXL (1, LL) .EQ. N2) THEN
            NXL (1, LL) = N0
         ELSEIF (NXL (2, LL) .EQ. N2) THEN
            NXL (2, LL) = N0
         ENDIF
  230 CONTINUE
      IF (GRAPH) THEN
         CALL LCOLOR ('WHITE')
         CALL SFLUSH
      ENDIF
C
C  FIX LXN ARRAY
C  UNHOOK L2 FROM N2, L1 FROM N0, AND L5 FROM NC2
C
      CALL DELLXN (MXND, LXN, NUID, NAVAIL, IAVAIL, N2,
     &   L2, NNN, ERR, NOROOM)
      IF ((NOROOM) .OR. (ERR)) THEN
         CALL MESAGE ('** PROBLEMS IN TUCK UNHOOKING L2 FROM N2 **')
         GOTO 290
      ENDIF
      CALL DELLXN (MXND, LXN, NUID, NAVAIL, IAVAIL, N0,
     &   L1, NNN, ERR, NOROOM)
      IF ((NOROOM) .OR. (ERR)) THEN
         CALL MESAGE ('** PROBLEMS IN TUCK UNHOOKING L1 FROM N0 **')
         GOTO 290
      ENDIF
      CALL DELLXN (MXND, LXN, NUID, NAVAIL, IAVAIL, NC2,
     &   L5, NNN, ERR, NOROOM)
      IF ((NOROOM) .OR. (ERR)) THEN
         CALL MESAGE ('** PROBLEMS IN TUCK UNHOOKING L1 FROM N0 **')
         GOTO 290
      ENDIF
C
C  ADD ALL LINES STILL HOOKED TO N2 TO THE LIST OF LINES FOR N0
C
      DO 240 I = 1, NL
         LL = L1LIST (I)
         IF (LL .NE. L2) THEN
            CALL ADDLXN (MXND, LXN, NUID, NAVAIL, IAVAIL,
     &         N0, LL, NNN, ERR, NOROOM)
            IF ((NOROOM) .OR. (ERR)) THEN
               CALL MESAGE ('** PROBLEMS IN TUCK HOOKING N2'' LINES'//
     &            ' TO N0 **')
               GOTO 290
            ENDIF
         ENDIF
  240 CONTINUE
C
C  DELETE N2
C
      DO 250 I = 1, 3
         LXN (I, N2) = 0
  250 CONTINUE
      LXN (4, N2) = IAVAIL
      IAVAIL = N2
      NAVAIL = NAVAIL+1
      NUID (N2) = 0
C
C  NOW DELETE THE OLD ELEMENT
C
      DO 260 I = 1, 4
         LXK (I, KOLD) = 0
  260 CONTINUE
C
C  NOW FIX THE LNODES ARRAY
C
      LNODES (3, N0) = LNODES (3, N2)
      LNODES (2, LNODES (3, N2) ) = N0
      LNODES (5, N0) = LNODES (5, N2)
C
      NLOOP = NLOOP - 2
      ERR = .FALSE.
C
C  NOW REDRAW THE ELEMENTS
C
      IF (GRAPH) THEN
         CALL LCOLOR ('BLACK')
         CALL D2NODE (MXND, XN, YN, N0, N1)
         CALL D2NODE (MXND, XN, YN, NC2, N1)
         CALL D2NODE (MXND, XN, YN, N2, N1)
         CALL LCOLOR ('WHITE')
         CALL GETLXN (MXND, LXN, N0, L1LIST, NL, ERR)
         IF (ERR) GOTO 290
         DO 270 II = 1, NL
            IDRAW = L1LIST (II)
            NODE1 = NXL (1, IDRAW)
            NODE2 = NXL (2, IDRAW)
            CALL D2NODE (MXND, XN, YN, NODE1, NODE2)
  270    CONTINUE
         CALL GETLXN (MXND, LXN, NC1, L1LIST, NL, ERR)
         IF (ERR) GOTO 290
         DO 280 II = 1, NL
            IDRAW = L1LIST (II)
            NODE1 = NXL (1, IDRAW)
            NODE2 = NXL (2, IDRAW)
            CALL D2NODE (MXND, XN, YN, NODE1, NODE2)
  280    CONTINUE
         CALL SFLUSH
      ENDIF
C
C  FLAG NODES FOR SMOOTHING
C
      CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES, NC1, ERR)
      IF (ERR) GOTO 290
      CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES, NC2, ERR)
      IF (ERR) GOTO 290
C
  290 CONTINUE
C
      RETURN
C
      END
