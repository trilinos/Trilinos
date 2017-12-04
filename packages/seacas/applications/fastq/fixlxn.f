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

C $Id: fixlxn.f,v 1.1 1990/11/30 11:07:27 gdsjaar Exp $
C $Log: fixlxn.f,v $
C Revision 1.1  1990/11/30 11:07:27  gdsjaar
C Initial revision
C
C
CC* FILE: [.PAVING]FIXLXN.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE FIXLXN (MXND, LXN, NXL, NUID, NAVAIL, IAVAIL, NNN, LLL,
     &   NNNOLD, LLLOLD, ERR, NOROOM)
C***********************************************************************
C
C  SUBROUTINE FIXLXN = FIXES THE ADDITIONS TO LXN
C
C***********************************************************************
C
      DIMENSION NXL (2, 3*MXND), LXN (4, MXND), NUID (MXND)
C
      LOGICAL ERR, NOROOM
C
C     RE-SETUP AVAILABLE LXN-SPACE LINKS
C
      IOLD = 0
      NAVAIL = 0
      DO 100 I = 1, NNNOLD
         IF (LXN (1, I).EQ.0)THEN
            IF (IOLD.LE.0)THEN
               IAVAIL = I
            ELSE
               LXN (4, IOLD) = I
            ENDIF
            IOLD = I
            NAVAIL = NAVAIL + 1
         ENDIF
  100 CONTINUE
      IF (IOLD.LE.0)THEN
         IAVAIL = NNN + 1
      ELSE
         LXN (4, IOLD) = NNN + 1
      ENDIF
      NAVAIL = NAVAIL + (MXND - NNN)
      IF (NNN.LT.MXND)THEN
         NNN1 = NNN + 1
         DO 110 I = NNN1, MXND
            LXN (1, I) = 0
            LXN (2, I) = 0
            LXN (3, I) = 0
            LXN (4, I) = I + 1
  110    CONTINUE
      ENDIF
C
C     COMPLETE LXN ARRAYS FOR ANY NEW LINES
C
      DO 130 L = LLLOLD + 1, LLL
         DO 120 I = 1, 2
            CALL ADDLXN (MXND, LXN, NUID, NAVAIL, IAVAIL, NXL (I, L),
     &         L, NNN, ERR, NOROOM)
            IF (ERR)THEN
               CALL MESAGE ('ERROR IN FIXLXN - NXL TABLE GENERATION')
               GOTO 140
            ELSEIF (NOROOM) THEN
               GOTO 140
            ENDIF
  120    CONTINUE
  130 CONTINUE
C
  140 CONTINUE
      RETURN
C
      END
