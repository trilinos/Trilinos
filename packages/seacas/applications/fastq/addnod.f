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

C $Id: addnod.f,v 1.1 1990/11/30 11:02:59 gdsjaar Exp $
C $Log: addnod.f,v $
C Revision 1.1  1990/11/30 11:02:59  gdsjaar
C Initial revision
C
C
CC* FILE: [.PAVING]ADDNOD.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/31/90
CC* MODIFICATION: ADDED ARGUMENTS TO CALL TO ADDNOD TO PASS MINIMUM
CC**              ELEMENT SIZE (SIZMIN) AND GETSIZ PARAMETERS OF
CC**              EMIN AND EMAX
C
      SUBROUTINE ADDNOD (MXND, MLN, XN, YN, LXK, KXL, NXL, LXN,
     &   ANGLE, BNSIZE, LNODES, XNEW, YNEW, DIST, NNN, KKK, LLL,
     &   N0, N1, N2, AMBIG, LAMBIG, SIZEIT, ERR, NOROOM, XNOLD,
     &   YNOLD, NXKOLD, LINKEG, LISTEG, BMESUR, MLINK, NPNOLD, NPEOLD,
     &   NNXK, REMESH, REXMIN, REXMAX, REYMIN, REYMAX, IDIVIS, SIZMIN,
     &   EMAX, EMIN)
C***********************************************************************
C
C  SUBROUTINE ADDNOD = ADDS A NEW ELEMENT TO A NEW NODE
C
C***********************************************************************
C
      DIMENSION XN (MXND), YN (MXND)
      DIMENSION LXK (4, MXND), KXL (2, 3*MXND)
      DIMENSION NXL (2, 3*MXND), LXN (4, MXND)
      DIMENSION ANGLE (MXND), LNODES (MLN, MXND), BNSIZE (2, MXND)
C
      DIMENSION XNOLD(NPNOLD), YNOLD(NPNOLD)
      DIMENSION NXKOLD(NNXK, NPEOLD)
      DIMENSION LINKEG(2, MLINK), LISTEG(4 * NPEOLD), BMESUR(NPNOLD)
C
      LOGICAL AMBIG, SIZEIT, ERR, NOROOM
C
      NNN = NNN+1
      IF (NNN .GT. MXND) THEN
         NOROOM = .TRUE.
         GOTO 110
      ENDIF
      XN (NNN) = XNEW
      YN (NNN) = YNEW
C
C  PUT THE BEGINNING BOUNDARY DISTANCE IN PLACE
C
      IF (LXN (2, N2) .LT. 0) THEN
         BNSIZE (1, NNN) = DIST
         BNSIZE (2, NNN) = 1.
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
     &         REYMIN, REYMAX, IDIVIS, SIZMIN, EMAX, EMIN, XNEW, YNEW,
     &         SIZE1)
            CALL GETSIZ (XNOLD, YNOLD, NXKOLD, LINKEG, LISTEG, BMESUR,
     &         MLINK, NPNOLD, NPEOLD, NNXK, REMESH, REXMIN, REXMAX,
     &         REYMIN, REYMAX, IDIVIS, SIZMIN, EMAX, EMIN, XN(N2),
     &         YN(N2), SIZE2)
            SIZNEW = AMIN1 (SIZE1, SIZE2)
         ELSE
            SIZNEW = BNSIZE (1, N2)
         ENDIF
         BNSIZE (1, NNN) = SIZNEW
         IF ((BNSIZE (1, N2) .EQ. 0.) .OR. (SIZEIT)) THEN
            BNSIZE (2, NNN) = 1.
         ELSE
            BNSIZE (2, NNN) = DIST / SIZNEW
         ENDIF
      ENDIF
C
C  MAKE LXN, NXL, KXL, AND LXK ARRAYS
C  FIRST, ADD THE NEW NODE'S LINES
C
      LLL = LLL+1
      NXL (1, LLL) = NNN
      NXL (2, LLL) = N0
      LLL = LLL+1
      NXL (1, LLL) = NNN
      NXL (2, LLL) = N2
      DO 100 I = 1, 4
         LXN (I, NNN) = 0
  100 CONTINUE
C
C  MAKE THE NEW ELEMENT
C
      KKK = KKK+1
      LXK (1, KKK) = LNODES (5, N0)
      IF (AMBIG) THEN
         LXK (2, KKK) = LAMBIG
      ELSE
         LXK (2, KKK) = LNODES (5, N1)
      ENDIF
      LXK (3, KKK) = LLL
      LXK (4, KKK) = LLL-1
      CALL ADDKXL (MXND, KXL, KKK, LLL)
      CALL ADDKXL (MXND, KXL, KKK, LLL-1)
      CALL ADDKXL (MXND, KXL, KKK, LNODES (5, N0))
      IF (AMBIG) THEN
         CALL ADDKXL (MXND, KXL, KKK, LAMBIG)
      ELSE
         CALL ADDKXL (MXND, KXL, KKK, LNODES (5, N1))
      ENDIF
C
C  REDO THE LNODES ARRAY
C
      LNODES (2, N2) = NNN
      LNODES (3, N0) = NNN
      LNODES (1, N0) = 0
      LNODES (1, NNN) = 0
      LNODES (1, N2) = 0
      LNODES (2, NNN) = N0
      LNODES (3, NNN) = N2
      LNODES (5, NNN) = LLL
      LNODES (5, N0) = LLL-1
      LNODES (4, NNN) = - 1
      CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES, N1, ERR)
      IF (ERR) GOTO 110
      IF (.NOT. AMBIG) LNODES (4, N1) = - 2
      LNODES (8, NNN) = LNODES (8, N2) + 1
C
  110 CONTINUE
      RETURN
C
      END
