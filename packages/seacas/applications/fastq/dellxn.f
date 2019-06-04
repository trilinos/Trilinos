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

C $Id: dellxn.f,v 1.1 1990/11/30 11:05:58 gdsjaar Exp $
C $Log: dellxn.f,v $
C Revision 1.1  1990/11/30 11:05:58  gdsjaar
C Initial revision
C

CC* FILE: [.QMESH]DELLXN.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE DELLXN (MXND, LXN, NUID, NAVAIL, IAVAIL, NODE, LINE,
     &   NNN, ERR, NOROOM)
C***********************************************************************
C
C  SUBROUTINE DELLXN = DELETE LINE FROM THE LIST OF LINES FOR THIS NODE
C
C***********************************************************************
C
      DIMENSION LINES (20), LXN (4, MXND), NUID (MXND)
C
      LOGICAL ERR, NOROOM
C
      CALL GETLXN (MXND, LXN, NODE, LINES, NL, ERR)
      IF (NL.LT.1) THEN
         WRITE (*, 10000)NODE
         GOTO 110
      ENDIF
      IF (ERR) GOTO 110
C
      K = 0
      DO 100 I = 1, NL
         IF (LINES (I) .NE. LINE) THEN
            K = K + 1
            LINES (K) = LINES (I)
         ENDIF
  100 CONTINUE
C
      IF (K .NE. NL - 1) THEN
         WRITE (*, 10010) NODE, (LINES (I), I = 1, NL)
         ERR = .TRUE.
         GOTO 110
      ENDIF
      NL = NL-1
      CALL PUTLXN (MXND, NL, LXN, NUID, NODE, LINES, NAVAIL, IAVAIL,
     &   NNN, ERR, NOROOM)
C
  110 CONTINUE
      RETURN
C
10000 FORMAT (' ERROR IN DELLXN - NODE', I5, ' HAS NO LINES')
10010 FORMAT (' ERROR IN DELLXN - NODE:', I5, /, ' LINES:', 20I5)
C
      END
