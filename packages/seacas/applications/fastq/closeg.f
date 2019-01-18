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

C $Id: closeg.f,v 1.1 1990/11/30 11:04:57 gdsjaar Exp $
C $Log: closeg.f,v $
C Revision 1.1  1990/11/30 11:04:57  gdsjaar
C Initial revision
C
C
CC* FILE: [.MAIN]CLOSEG.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE CLOSEG (MSNAP, SNAPDX, NSNAP, X, Y, II, INDEX, XBOT,
     &   XTOP, YBOT, YTOP)
C***********************************************************************
C
C  SUBROUTINE CLOSEG = SUBROUTINE TO RETURN CLOSEST GRID LINE
C
C***********************************************************************
C
C  SUBROUTINE CALLED BY:
C     DIGIT = A SUBROUTINE TO INPUT GEOMETRY
C
C***********************************************************************
C
C  VARIABLES USED:
C     X      = THE X LOCATION IN USER COORDINATES
C     Y      = THE Y LOCATION IN USER COORDINATES
C
C***********************************************************************
C
      DIMENSION SNAPDX(2, MSNAP), NSNAP(2)
C
C  FIND CLOSEST GRID CROSSING IN X OR Y
C
      XHOLD = X
      YHOLD = Y
      CALL SNAPPT (MSNAP, SNAPDX, NSNAP, XHOLD, YHOLD)
      IF (ABS(XHOLD - X) .LT. ABS(YHOLD - Y)) THEN
         INDEX = 1
      ELSE
         INDEX = 2
         XHOLD = YHOLD
      END IF
C
C  FIND INDEX TO GRID LINE
C
      DO 100 I = 1, NSNAP(INDEX)
         IF (SNAPDX(INDEX, I) .GE. XHOLD) THEN
            II = I
            GO TO 110
         END IF
  100 CONTINUE
      II = NSNAP(INDEX)
  110 CONTINUE
C
C  SET GRID LINE LIMITS
C
      IF (INDEX .EQ. 1) THEN
         XBOT = SNAPDX(1, II)
         XTOP = XBOT
         YBOT = SNAPDX(2, 1)
         YTOP = SNAPDX(2, NSNAP(2))
      ELSE
         XBOT = SNAPDX(1, 1)
         XTOP = SNAPDX(1, NSNAP(1))
         YBOT = SNAPDX(2, II)
         YTOP = YBOT
      END IF
C
      RETURN
C
      END
