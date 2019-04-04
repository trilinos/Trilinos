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

C $Id: mak2el.f,v 1.1 1990/11/30 11:11:48 gdsjaar Exp $
C $Log: mak2el.f,v $
C Revision 1.1  1990/11/30 11:11:48  gdsjaar
C Initial revision
C
C
CC* FILE: [.QMESH]MAK2EL.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE MAK2EL (MP, MXNPER, MXND, NNN0, NNN, KKK, X, Y, NID,
     &   XN, YN, NUID, LXK, COOR, IP3)
C***********************************************************************
C
C  SUBROUTINE MAK2EL = GENERATES  (ADDS) ELEMENT CONNECTIVITY FOR 2 NODES
C
C***********************************************************************
C
      DIMENSION COOR (2, MP), X (MXNPER), Y (MXNPER), NID (MXNPER)
      DIMENSION XN (MXND), YN (MXND), NUID (MXND), LXK (4, MXND)
C
C  PUT NODES AND NUID'S INTO THE PROPER LOCATIONS
C
      KOUNT = 0
      DO 100 I = NNN0, NNN-1
         KOUNT = KOUNT + 1
         KKK = KKK + 1
         XN (I) = X (KOUNT)
         YN (I) = Y (KOUNT)
         NUID (I) = NID (KOUNT)
         LXK (1, KKK) = I
         LXK (2, KKK) = I + 1
         LXK (3, KKK) = 0
         LXK (4, KKK) = 0
         IF (IP3.GT.0)THEN
            X1 = X (I + 1) - X (I)
            Y1 = Y (I + 1) - Y (I)
            X2 = X (I + 1) - COOR (1, IP3)
            Y2 = Y (I + 1) - COOR (2, IP3)
            CROSSP =  (X1 * Y2) -  (Y1 * X2)
            IF (CROSSP.GT.0)THEN
               LXK (1, KKK) = I + 1
               LXK (2, KKK) = I
            ENDIF
            IF (CROSSP * CROSSP .LT.
     &         (.01 * ((X1 * X1) +  (Y1 * Y1)) *
     &         ((X2 * X2) + (Y2 * Y2)))) WRITE (*, 10000) KKK
         ENDIF
  100 CONTINUE
C
      XN (NNN) = X (KOUNT + 1)
      YN (NNN) = Y (KOUNT + 1)
      NUID (NNN) = NID (KOUNT + 1)
C
      RETURN
C
10000 FORMAT (' ** WARNING **  -  COLINEAR REFERENCE NODE FOR ELEMENT:',
     &   I5)
      END
