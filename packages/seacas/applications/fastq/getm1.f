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

C $Id: getm1.f,v 1.1 1990/11/30 11:08:22 gdsjaar Exp $
C $Log: getm1.f,v $
C Revision 1.1  1990/11/30 11:08:22  gdsjaar
C Initial revision
C
C
CC* FILE: [.QMESH]GETM1.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE GETM1(ML, MS, MNNPS, NS, ISLIST, NINT, IFLINE, NLPS,
     &   ILLIST, LINKL, LINKS, X, Y, NID, NNPS, ANGLE, NPER, SCHSTR,
     &   M1, CCW, NORM, REAL, ERR)
C***********************************************************************
C
C  SUBROUTINE GETM1 = GETS THE APPROPRIATE M1 VALUE
C
C***********************************************************************
C
C  SUBROUTINE CALLED BY:
C     QMESH = GENERATES QUAD ELEMENTS
C
C***********************************************************************
C
C  VARIABLES USED:
C     NNPS  = ARRAY OF NUMBER OF NODES PER SIDE
C     CCW   = .TRUE. IF THE SIDE IS ORIENTED CCW
C     NORM  = .TRUE. IF THE FIRST SIDE IS TO BE TRIED AS THE BASE
C
C***********************************************************************
C
      DIMENSION NNPS(MNNPS), ISLIST(NS), LINKL(2, ML), LINKS(MS*2)
      DIMENSION NINT(ML), NLPS(MS), IFLINE(MS), ILLIST(MS*3)
      DIMENSION X(NPER), Y(NPER), NID(NPER), ANGLE(NPER)
C
      LOGICAL CCW, ERR, NORM, REAL
C
      CHARACTER*72 SCHSTR
C
C  SEE IF LETTER M OCCURS IN THE SCHEME (NOT THE "NORMAL" CASE)
C
      NORM = .TRUE.
      DO 100 J = 1, 72
         IF ( (SCHSTR (J:J) .EQ. 'M') .OR.
     &      (SCHSTR (J:J) .EQ. 'm') ) THEN
            NORM = .FALSE.
            GOTO 110
         ENDIF
  100 CONTINUE
  110 CONTINUE
C
C  CALCULATE THE NUMBER OF NODES PER SIDE
C
      CALL NPS (ML, MS, MNNPS, NS, ISLIST, NINT, IFLINE, NLPS, ILLIST,
     &   LINKL, LINKS, NNPS, ERR)
      IF (ERR) RETURN
C
C  NORMAL CASE - TRY THE FIRST SIDE IN THE SIDE LIST AS THE M1 BASE
C  IN CCW ORDER IF IT IS NOT RIDICULOUS
C
      IF (NORM) THEN
         M1 = NNPS (1) - 1
         IF (.NOT. CCW) M1 = NPER / 2 - M1
         IF ( (M1 .GT. 0).AND. (M1 .LE. NPER/2) ) RETURN
      ENDIF
      IF (.NOT. CCW) CALL IREVER (NNPS, NS)
C
C  IF THE BOUNDARY IS A LOGICAL RECTANGLE,  USE IT
C
      IF ( (NS .EQ. 4) .AND. (NNPS (1) .EQ. NNPS (3)) .AND.
     &   (NNPS (2) .EQ. NNPS (4) ) ) THEN
         M1 = NNPS (1) - 1
C
C  OTHERWISE,  FIND AN INITIAL M1 FOR A NON-LOGICAL RECTANGLE
C
      ELSE
         CALL PICKM1 (NPER, X, Y, ANGLE, M1, IFIRST, REAL)
         IF (IFIRST .NE. 1) CALL FQ_ROTATE (NPER, X, Y, NID, IFIRST)
      ENDIF
C
      RETURN
C
      END
