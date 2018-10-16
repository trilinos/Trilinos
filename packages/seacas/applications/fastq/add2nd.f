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

C $Id: add2nd.f,v 1.1 1990/11/30 11:02:46 gdsjaar Exp $
C $Log: add2nd.f,v $
C Revision 1.1  1990/11/30 11:02:46  gdsjaar
C Initial revision
C
C
CC* FILE: [.PAVING]ADD2ND.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/31/90
CC* MODIFICATION: ADDED ARGUMENTS TO CALL TO ADD2ND TO PASS MINIMUM
CC**              ELEMENT SIZE (SIZMIN) AND GETSIZ PARAMETERS OF
CC**              EMIN AND EMAX
C
      SUBROUTINE ADD2ND (MXND, MLN, XN, YN, LXK, KXL, NXL, LXN,
     &   BNSIZE, LNODES, X1, Y1, X2, Y2, DIST1, DIST2, NNN, LLL, KKK,
     &   N1, N2, NLOOP, SIZEIT, ERR, NOROOM, XNOLD, YNOLD, NXKOLD,
     &   LINKEG, LISTEG, BMESUR, MLINK, NPNOLD, NPEOLD, NNXK, REMESH,
     &   REXMIN, REXMAX, REYMIN, REYMAX, IDIVIS, SIZMIN, EMAX, EMIN)
C***********************************************************************
C
C  SUBROUTINE ADD2ND = ADDS A NEW ELEMENT JUTTING OUT FROM AN EXISTING
C                      LINE
C
C***********************************************************************
C
      DIMENSION XN (MXND), YN (MXND)
      DIMENSION LXK (4, MXND), KXL (2, 3*MXND)
      DIMENSION NXL (2, 3*MXND), LXN (4, MXND)
      DIMENSION LNODES (MLN, MXND), BNSIZE (2, MXND)
C
      DIMENSION XNOLD(NPNOLD), YNOLD(NPNOLD)
      DIMENSION NXKOLD(NNXK, NPEOLD)
      DIMENSION LINKEG(2, MLINK), LISTEG(4 * NPEOLD), BMESUR(NPNOLD)
C
      LOGICAL SIZEIT, ERR, NOROOM
C
      NNN = NNN+1
      IF (NNN .GT. MXND) THEN
         NOROOM = .TRUE.
         GOTO 110
      ENDIF
      XN (NNN) = X1
      YN (NNN) = Y1
      NODE1 = NNN
      NNN = NNN+1
      IF (NNN .GT. MXND) THEN
         NOROOM = .TRUE.
         GOTO 110
      ENDIF
      XN (NNN) = X2
      YN (NNN) = Y2
      NODE2 = NNN
C
C  PUT THE BEGINNING BOUNDARY DISTANCE IN PLACE
C
      IF (LXN (2, N1) .LT. 0) THEN
         BNSIZE (1, NODE1) = DIST1
         BNSIZE (2, NODE1) = 1.
      ELSE
         IF (SIZEIT) THEN
C
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/21/90
CC* MODIFICATION: CHANGED PROJECTION SIZE TO CHOOSE MINIMUM OF CURRENT
C**               LOCATION SIZE AND PROJECTING FROM LOCATION SIZE.
C
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/31/90
CC* MODIFICATION: ADDED ARGUMENTS TO CALL TO GETSIZ TO PASS MINIMUM
CC**              ELEMENT SIZE (SIZMIN) AND GETSIZ PARAMETERS OF
CC**              EMIN AND EMAX
C
            CALL GETSIZ (XNOLD, YNOLD, NXKOLD, LINKEG, LISTEG, BMESUR,
     &         MLINK, NPNOLD, NPEOLD, NNXK, REMESH, REXMIN, REXMAX,
     &         REYMIN, REYMAX, IDIVIS, SIZMIN, EMAX, EMIN, X1, Y1,
     &         SIZE1)
            CALL GETSIZ (XNOLD, YNOLD, NXKOLD, LINKEG, LISTEG, BMESUR,
     &         MLINK, NPNOLD, NPEOLD, NNXK, REMESH, REXMIN, REXMAX,
     &         REYMIN, REYMAX, IDIVIS, SIZMIN, EMAX, EMIN, XN(N1),
     &         YN(N1), SIZE2)
            SIZNEW = AMIN1 (SIZE1, SIZE2)
         ELSE
            SIZNEW = BNSIZE (1, N1)
         ENDIF
         BNSIZE (1, NODE1) = SIZNEW
         IF ((BNSIZE (1, N1) .EQ. 0.) .OR. (SIZEIT)) THEN
            BNSIZE (2, NODE1) = 1.
         ELSE
            BNSIZE (2, NODE1) = DIST1 / SIZNEW
         ENDIF
      ENDIF
      IF (LXN (2, N2) .LT. 0) THEN
         BNSIZE (1, NODE2) = DIST2
         BNSIZE (2, NODE2) = 1.
      ELSE
         IF (SIZEIT) THEN
            CALL GETSIZ (XNOLD, YNOLD, NXKOLD, LINKEG, LISTEG, BMESUR,
     &         MLINK, NPNOLD, NPEOLD, NNXK, REMESH, REXMIN, REXMAX,
     &         REYMIN, REYMAX, IDIVIS, SIZMIN, EMAX, EMIN, X2, Y2,
     &         SIZNEW)
         ELSE
            SIZNEW = BNSIZE (1, N2)
         ENDIF
         BNSIZE (1, NODE2) = SIZNEW
         IF ((BNSIZE (1, N2) .EQ. 0.) .OR. (SIZEIT)) THEN
            BNSIZE (2, NODE2) = 1.
         ELSE
            BNSIZE (2, NODE2) = DIST2 / SIZNEW
         ENDIF
      ENDIF
C
C  MAKE NXL ARRAY
C  ADD THE THREE NEW LINES
C
      LLL = LLL+1
      L1 = LLL
      NXL (1, L1) = N1
      NXL (2, L1) = NODE1
C
      LLL = LLL+1
      L2 = LLL
      NXL (1, L2) = NODE1
      NXL (2, L2) = NODE2
C
      LLL = LLL+1
      L3 = LLL
      NXL (1, L3) = NODE2
      NXL (2, L3) = N2
C
C  MAKE THE NEW ELEMENT
C
      KKK = KKK+1
      LXK (1, KKK) = LNODES (5, N1)
      LXK (2, KKK) = L3
      LXK (3, KKK) = L2
      LXK (4, KKK) = L1
      CALL ADDKXL (MXND, KXL, KKK, L1)
      CALL ADDKXL (MXND, KXL, KKK, L2)
      CALL ADDKXL (MXND, KXL, KKK, L3)
      CALL ADDKXL (MXND, KXL, KKK, LNODES (5, N1))
C
C  ZERO OUT THE LXN ARRAY
C
      DO 100 I = 1, 4
         LXN (I, NODE1) = 0
         LXN (I, NODE2) = 0
  100 CONTINUE
C
C  REDO THE LNODES ARRAY
C
      LNODES (1, NODE1) = 0
      LNODES (1, NODE2) = 0
      LNODES (1, N1) = 0
      LNODES (1, N2) = 0
C
      LNODES (2, NODE1) = N1
      LNODES (2, NODE2) = NODE1
      LNODES (2, N2) = NODE2
C
      LNODES (3, N1) = NODE1
      LNODES (3, NODE1) = NODE2
      LNODES (3, NODE2) = N2
C
      LNODES (4, NODE1) = - 1
      LNODES (4, NODE2) = - 1
C
      LNODES (5, N1) = L1
      LNODES (5, NODE1) = L2
      LNODES (5, NODE2) = L3
C
      LNODES (8, NODE1) = LNODES (8, N1) + 1
      LNODES (8, NODE2) = LNODES (8, N2) + 1
C
      NLOOP = NLOOP + 2
C
      CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES, N1, ERR)
      CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES, N2, ERR)
C
  110 CONTINUE
      RETURN
C
      END
