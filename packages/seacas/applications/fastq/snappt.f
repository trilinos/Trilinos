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

C $Id: snappt.f,v 1.1 1990/11/30 11:15:56 gdsjaar Exp $
C $Log: snappt.f,v $
C Revision 1.1  1990/11/30 11:15:56  gdsjaar
C Initial revision
C
C
CC* FILE: [.MAIN]SNAPPT.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE SNAPPT (MSNAP, SNAPDX, NSNAP, X, Y)
C***********************************************************************
C
C  SUBROUTINE SNAPPT = GETS THE X, Y TO THE CLOSEST GRID POINT
C
C***********************************************************************
C
C  VARIABLES USED:
C     SNAPDX = SNAP GRID LINE ARRAY
C     NSNAP  = ARRAY OF GRID LINE COUNTERS
C     MSNAP  = THE DIMENSION OF THE GRID LINE ARRAY
C     XMID   = .TRUE. IF THE X VALUE FALLS BETWEEN TWO X GRID LINES
C     YMID   = .TRUE. IF THE Y VALUE FALLS BETWEEN TWO Y GRID LINES
C
C***********************************************************************
C
      DIMENSION SNAPDX (2, MSNAP), NSNAP (2)
C
      LOGICAL XMID,  YMID
C
      XMID = .FALSE.
      YMID = .FALSE.
C
C  GET THE BOUNDING X SNAP LINES
C
      IF (X .LE. SNAPDX (1, 1)) THEN
         XP = SNAPDX (1, 1)
      ELSEIF (X .GE. SNAPDX (1, NSNAP (1))) THEN
         XP = SNAPDX (1, NSNAP (1))
      ELSE
         XMID = .TRUE.
         DO 100 I = 2, NSNAP (1)
            IF (X .LE. SNAPDX (1, I)) THEN
               XP1 = SNAPDX (1, I - 1)
               XP2 = SNAPDX (1, I)
               GOTO 110
            ENDIF
  100    CONTINUE
  110    CONTINUE
      ENDIF
C
C  GET THE BOUNDING Y SNAP LINES
C
      IF (Y .LE. SNAPDX (2, 1)) THEN
         YP = SNAPDX (2, 1)
      ELSEIF (Y .GE. SNAPDX (2, NSNAP (2))) THEN
         YP = SNAPDX (2, NSNAP (2))
      ELSE
         YMID = .TRUE.
         DO 120 I = 2, NSNAP (2)
            IF (Y .LE. SNAPDX (2, I)) THEN
               YP1 = SNAPDX (2, I - 1)
               YP2 = SNAPDX (2, I)
               GOTO 130
            ENDIF
  120    CONTINUE
  130    CONTINUE
      ENDIF
C
C  NOW GET THE APPROPRIATE COMBINATION OF XLOW, XHIGH, XMID, YLOW,  ETC.
C
C
C  FIRST THE MOST COMMON CASE OF FITTING BETWEEN X AND Y GRIDS
C
      IF ( (YMID) .AND. (XMID)) THEN
C
C  GET THE SHORTEST DISTANCE TO THIS COMBINATION
C
         DIST1 = SQRT ( ( (XP1 - X) ** 2) +  ( (YP1 - Y) ** 2))
         DIST2 = SQRT ( ( (XP1 - X) ** 2) +  ( (YP2 - Y) ** 2))
         DIST3 = SQRT ( ( (XP2 - X) ** 2) +  ( (YP1 - Y) ** 2))
         DIST4 = SQRT ( ( (XP2 - X) ** 2) +  ( (YP2 - Y) ** 2))
C
         IF (DIST1 .LE. AMIN1 (DIST2, DIST3, DIST4)) THEN
            X = XP1
            Y = YP1
         ELSEIF (DIST2 .LE. AMIN1 (DIST1, DIST3, DIST4)) THEN
            X = XP1
            Y = YP2
         ELSEIF (DIST3 .LE. AMIN1 (DIST1, DIST2, DIST4)) THEN
            X = XP2
            Y = YP1
         ELSE
            X = XP2
            Y = YP2
         ENDIF
C
C  NOW THE CORNER CASES OF XLOW,  XHIGH,  YLOW,  AND YHIGH COMBINATIONS
C
      ELSEIF ( (.NOT.XMID) .AND. (.NOT.YMID)) THEN
         X = XP
         Y = YP
C
C  NOW THE EDGE CASES OF XLOW OR XHIGH AND YMID
C
      ELSEIF (.NOT.XMID) THEN
         X = XP
         DIST1 = SQRT ( ( (XP - X) ** 2) +  ( (YP1 - Y) ** 2))
         DIST2 = SQRT ( ( (XP - X) ** 2) +  ( (YP2 - Y) ** 2))
         IF (DIST1 .LT. DIST2) THEN
            Y = YP1
         ELSE
            Y = YP2
         ENDIF
C
C  NOW THE EDGE CASES OF XMID AND YHIGH OR YLOW
C
      ELSEIF (.NOT.YMID) THEN
         Y = YP
         DIST1 = SQRT ( ( (XP1 - X) ** 2) +  ( (YP - Y) ** 2))
         DIST2 = SQRT ( ( (XP2 - X) ** 2) +  ( (YP - Y) ** 2))
         IF (DIST1 .LT. DIST2) THEN
            X = XP1
         ELSE
            X = XP2
         ENDIF
C
C  NOW A CHECK TO MAKE SURE THAT SOMETHING DIDN'T FALL THROUGH
C
      ELSE
         CALL MESAGE (' **  ERROR  -  IMPOSSIBLE CASE IN SNAPPT  ** ')
         CALL PLTBEL
      ENDIF
C
      RETURN
C
      END
