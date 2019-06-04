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

C $Id: undelm.f,v 1.2 1998/07/14 18:20:13 gdsjaar Exp $
C $Log: undelm.f,v $
C Revision 1.2  1998/07/14 18:20:13  gdsjaar
C Removed unused variables, cleaned up a little.
C
C Changed BLUE labels to GREEN to help visibility on black background
C (indirectly requested by a couple users)
C
C Revision 1.1.1.1  1990/11/30 11:17:29  gdsjaar
C FASTQ Version 2.0X
C
c Revision 1.1  90/11/30  11:17:27  gdsjaar
c Initial revision
c
C
CC* FILE: [.PAVING]UNDELM.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE UNDELM (MXND, MLN, LNODES, XN, YN, NUID, LXK, KXL, NXL,
     &   LXN, NNN, LLL, KKK, NAVAIL, IAVAIL, N0, N1, N2, N3, L1, L2, L3,
     &   K1, K2, NOROOM, ERR, GRAPH, VIDEO)
C***********************************************************************
C
C  SUBROUTINE UNDELM = UNDELETES AN ELEMENT BY EXPANDING N1 INTO A
C                      NEW ELEMENT
C
C***********************************************************************
C
      DIMENSION LXK(4, MXND), NXL(2, 3*MXND), KXL(2, 3*MXND)
      DIMENSION LXN(4, MXND), XN(MXND), YN(MXND), NUID(MXND)
      DIMENSION LNODES (MLN, MXND)
C
      LOGICAL ERR, NOROOM, GRAPH, VIDEO
C
      ERR = .FALSE.
C
C  MAKE SURE THAT N2 HAS AT LEAST FOUR LINES ATTACHED TO IT
C
      IF (LXN (4, N2) .EQ. 0) THEN
         ERR = .TRUE.
         CALL MESAGE ('** N2 IN UNDELM CANNOT BE USED'//
     &      ' TO EXPAND AN ELEMNT **')
         GOTO 140
      ENDIF
C
C  ERASE THE LINE L3
C
      IF ((GRAPH) .OR. (VIDEO)) THEN
         CALL LCOLOR ('BLACK')
         CALL D2NODE (MXND, XN, YN, N0, N2)
         IF (GRAPH) THEN
            CALL LCOLOR ('WHITE')
         ELSE
            CALL LCOLOR ('YELOW')
         ENDIF
         CALL SFLUSH
         IF (VIDEO) CALL SNAPIT (3)
      ENDIF
C
C  DEFINE THE NEW NODE AND THE TWO NEW LINES
C
      NNN = NNN + 1
      IF (NNN .GT. MXND) THEN
         NOROOM = .TRUE.
         GOTO 140
      ENDIF
      NNEW = NNN
      XN (NNEW) = (XN (N0) + XN (N2)) * .5
      YN (NNEW) = (YN (N0) + YN (N2)) * .5
C
      LLL = LLL + 2
      L4 = LLL -1
      L5 = LLL
      NXL (1, L4) = N1
      NXL (2, L4) = NNEW
      NXL (1, L5) = N3
      NXL (2, L5) = NNEW
C
C  NOW CHANGE LINE L3'S END POINT FROM N2 TO NNEW
C
      IF (NXL (1, L3) .EQ. N2) THEN
         NXL (1, L3) = NNEW
      ELSEIF (NXL (2, L3) .EQ. N2) THEN
         NXL (2, L3) = NNEW
      ELSE
         CALL MESAGE ('** PROBLEMS IN UNDLEM WITH L3''S END POINT **')
         ERR = .TRUE.
         GOTO 140
      ENDIF
C
C  NOW UPDATE THE LXN ARRAYS
C
      LXN (1, NNEW) = L3
      LXN (2, NNEW) = L4
      LXN (3, NNEW) = L5
      LXN (4, NNEW) = 0
      CALL FIXLXN (MXND, LXN, NXL, NUID, NAVAIL, IAVAIL, NNN, LLL,
     &   NNN, LLL, ERR, NOROOM)
      IF ((NOROOM) .OR. (ERR)) GOTO 140
C
C  REMOVE L3 FROM THE LIST OF LINES FOR N2
C
      CALL DELLXN (MXND, LXN, NUID, NAVAIL, IAVAIL, N2,
     &   L3, NNN, ERR, NOROOM)
      IF ((NOROOM) .OR. (ERR)) THEN
         CALL MESAGE ('** PROBLEMS IN UNDELM UNHOOKING L3 FROM N2 **')
         GOTO 140
      ENDIF
C
C  ADD LINE L4 TO N1
C
      CALL ADDLXN (MXND, LXN, NUID, NAVAIL, IAVAIL,
     &   N1, L4, NNN, ERR, NOROOM)
      IF ((NOROOM) .OR. (ERR)) THEN
         CALL MESAGE ('** PROBLEMS IN UNDELM HOOKING L4 TO N1 **')
         GOTO 140
      ENDIF
C
C  ADD LINE L5 TO N3
C
      CALL ADDLXN (MXND, LXN, NUID, NAVAIL, IAVAIL,
     &   N3, L5, NNN, ERR, NOROOM)
      IF ((NOROOM) .OR. (ERR)) THEN
         CALL MESAGE ('** PROBLEMS IN UNDELM HOOKING L5 TO N3 **')
         GOTO 140
      ENDIF
C
C  NOW ADD THE NEW ELEMENT
C
      KKK = KKK + 1
      LXK (1, KKK) = L1
      LXK (2, KKK) = L2
      LXK (3, KKK) = L5
      LXK (4, KKK) = L4
C
C  NOW FIX THE KXL ARRAY FOR LINE L1
C
      IF (KXL (1, L1) .EQ. K2) THEN
         KXL (1, L1) = KKK
      ELSEIF (KXL (2, L1) .EQ. K2) THEN
         KXL (2, L1) = KKK
      ELSE
         CALL MESAGE ('** PROBLEMS IN UNDELM REPLACING K2 FOR L1 **')
         ERR = .TRUE.
         GOTO 140
      ENDIF
C
C  NOW FIX THE KXL ARRAY FOR LINE L2
C
      IF (KXL (1, L2) .EQ. K1) THEN
         KXL (1, L2) = KKK
      ELSEIF (KXL (2, L2) .EQ. K1) THEN
         KXL (2, L2) = KKK
      ELSE
         CALL MESAGE ('** PROBLEMS IN UNDELM REPLACING K1 FOR L2 **')
         ERR = .TRUE.
         GOTO 140
      ENDIF
C
C  ADD THE KXL ENTRIES FOR THE NEW LINES
C
      KXL (1, L4) = K2
      KXL (2, L4) = KKK
      KXL (1, L5) = K1
      KXL (2, L5) = KKK
C
C  NOW FIX THE LXK ARRAY FOR THE ELEMENT K1
C
      DO 100 I = 1, 4
         IF (LXK (I, K1) .EQ. L2) THEN
            LXK (I, K1) = L5
            GOTO 110
         ENDIF
  100 CONTINUE
      CALL MESAGE ('** PROBLEMS IN UNDELM REPLACING L2 WITH L5 IN '//
     &   'K1 **')
      ERR = .TRUE.
      GOTO 140
  110 CONTINUE
C
C  NOW FIX THE LXK ARRAY FOR THE ELEMENT K2
C
      DO 120 I = 1, 4
         IF (LXK (I, K2) .EQ. L1) THEN
            LXK (I, K2) = L4
            GOTO 130
         ENDIF
  120 CONTINUE
      CALL MESAGE ('** PROBLEMS IN UNDELM REPLACING L1 WITH L4 IN '//
     &   'K2 **')
      ERR = .TRUE.
      GOTO 140
  130 CONTINUE
C
C  NOW REDRAW THE LINES
C
      IF ((GRAPH) .OR. (VIDEO)) THEN
         CALL D2NODE (MXND, XN, YN, N0, NNEW)
         CALL D2NODE (MXND, XN, YN, N1, NNEW)
         CALL D2NODE (MXND, XN, YN, N3, NNEW)
         CALL SFLUSH
         IF (VIDEO) CALL SNAPIT (3)
      ENDIF
C
      LNODES (4, NNEW) = 2
      CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES,
     &   NNEW, ERR)
      IF (ERR) GOTO 140
      CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES,
     &   N0, ERR)
      IF (ERR) GOTO 140
      CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES,
     &   N1, ERR)
      IF (ERR) GOTO 140
      CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES,
     &   N2, ERR)
      IF (ERR) GOTO 140
      CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES,
     &   N3, ERR)
      IF (ERR) GOTO 140
  140 CONTINUE
      RETURN
C
      END
