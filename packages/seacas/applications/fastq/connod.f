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

C $Id: connod.f,v 1.1 1990/11/30 11:05:18 gdsjaar Exp $
C $Log: connod.f,v $
C Revision 1.1  1990/11/30 11:05:18  gdsjaar
C Initial revision
C
C
CC* FILE: [.PAVING]CONNOD.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE CONNOD (MXND, MLN, XN, YN, NUID, LXK, KXL, NXL, LXN,
     &   ANGLE, LNODES, NNN, KKK, LLL, NNNOLD, LLLOLD, N0, N1, N2,
     &   NODE, I, NLOOP, IAVAIL, NAVAIL, GRAPH, VIDEO, NOROOM, ERR)
C***********************************************************************
C
C  SUBROUTINE CONNOD = ADDS A NEW ELEMENT TO AN EXISTING NODE
C
C***********************************************************************
C
      DIMENSION XN (MXND), YN (MXND), NUID (MXND)
      DIMENSION LXK (4, MXND), KXL (2, 3*MXND)
      DIMENSION NXL (2, 3*MXND), LXN (4, MXND)
      DIMENSION ANGLE (MXND), LNODES (MLN, MXND)
C
      LOGICAL GRAPH, VIDEO, ERR, NOROOM
C
      ERR = .FALSE.
C
      NLOOP = NLOOP-2
C
C  TAKE CARE OF A COUNTERCLOCKWISE CONNECTION
C
      IF (LNODES (3, N2) .EQ. NODE) THEN
         LNODES (3, N0) = NODE
         LNODES (2, NODE) = N0
         LNODES (4, N1) = - 2
         LNODES (4, N2) = - 2
         ANGLE (N1) = 0.
         ANGLE (N2) = 0.
         LLL = LLL+1
         KKK = KKK+1
         NXL (1, LLL) = NODE
         NXL (2, LLL) = N0
C
C  MAKE THE NEW ELEMENT
C
         LXK (1, KKK) = LNODES (5, N0)
         LXK (2, KKK) = LNODES (5, N1)
         LXK (3, KKK) = LNODES (5, N2)
         LXK (4, KKK) = LLL
C
         CALL ADDKXL (MXND, KXL, KKK, LLL)
         CALL ADDKXL (MXND, KXL, KKK, LNODES (5, N0))
         CALL ADDKXL (MXND, KXL, KKK, LNODES (5, N1))
         CALL ADDKXL (MXND, KXL, KKK, LNODES (5, N2))
C
         LNODES (5, N0) = LLL
         IF ((GRAPH) .OR. (VIDEO)) THEN
            CALL D2NODE (MXND, XN, YN, N0, NODE)
            CALL SFLUSH
            IF (VIDEO) CALL SNAPIT (1)
         ENDIF
         I = NODE
C
C  TAKE CARE OF A CLOCKWISE CONNECTION
C
      ELSE IF (LNODES (2, N0) .EQ. NODE) THEN
         LNODES (3, NODE) = N2
         LNODES (2, N2) = NODE
         LNODES (4, N0) = - 2
         LNODES (4, N1) = - 2
         ANGLE (N0) = 0.
         ANGLE (N1) = 0.
         LLL = LLL+1
         KKK = KKK+1
         NXL (1, LLL) = NODE
         NXL (2, LLL) = N2
C
C  MAKE THE NEW ELEMENT
C
         LXK (1, KKK) = LNODES (5, NODE)
         LXK (2, KKK) = LNODES (5, N0)
         LXK (3, KKK) = LNODES (5, N1)
         LXK (4, KKK) = LLL
C
         CALL ADDKXL (MXND, KXL, KKK, LNODES (5, NODE))
         CALL ADDKXL (MXND, KXL, KKK, LNODES (5, N0))
         CALL ADDKXL (MXND, KXL, KKK, LNODES (5, N1))
         CALL ADDKXL (MXND, KXL, KKK, LLL)
C
         LNODES (5, NODE) = LLL
C
         IF ((GRAPH) .OR. (VIDEO)) THEN
            CALL D2NODE (MXND, XN, YN, NODE, N2)
            CALL SFLUSH
            IF (VIDEO) CALL SNAPIT (1)
         ENDIF
         I = N2
C
      ELSE
         CALL MESAGE ('PROBLEM HANDLING CONNECTION IN CONNOD')
         CALL MESAGE ('CONNECTION DOES NOT FORM A SIMPLE CORNER')
         ERR = .TRUE.
         GOTO 100
      ENDIF
      CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES, N1, ERR)
      IF (ERR) GOTO 100
      CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES, N2, ERR)
      IF (ERR) GOTO 100
      CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES, NODE, ERR)
      IF (ERR) GOTO 100
C
      CALL FIXLXN (MXND, LXN, NXL, NUID, NAVAIL, IAVAIL, NNN, LLL,
     &   NNNOLD, LLLOLD, ERR, NOROOM)
      IF (ERR) GOTO 100
      LLLOLD = LLL
      NNNOLD = NNN
C
  100 CONTINUE
C
      RETURN
C
      END
