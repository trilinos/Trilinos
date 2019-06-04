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

C $Id: closep.f,v 1.1 1990/11/30 11:05:03 gdsjaar Exp $
C $Log: closep.f,v $
C Revision 1.1  1990/11/30 11:05:03  gdsjaar
C Initial revision
C
C
CC* FILE: [.MAIN]CLOSEP.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE CLOSEP (MP, N15, X, Y, IPOINT, COOR, LINKP, JJ)
C***********************************************************************
C
C  SUBROUTINE CLOSE = FINDS THE CLOSEST EXISTING POINT TO THE MOUSE
C
C***********************************************************************
C
C  SUBROUTINE CALLED BY:
C     INPUT  = INPUTS MESH DEFINITIONS FROM THE LIGHT TABLE
C
C***********************************************************************
C
C  VARIABLES USED:
C     X      = THE X LOCATION IN USER COORDINATES
C     Y      = THE Y LOCATION IN USER COORDINATES
C     POINT  = ARRAY OF VALUES DEFINING A POINT
C               (I, 1) = THE NUMBER OF THE POINT
C               (I, 2) = THE X COORDINATE OF THE POINT
C               (I, 3) = THE Y COORDINATE OF THE POINT
C               (I, 4) = THE BOUNDARY FLAG OF THE POINT
C     I      = THE NUMBER OF THE CLOSEST POINT FOUND
C     K      = THE NUMBER OF POINTS IN THE DATABASE
C
C***********************************************************************
C
      DIMENSION IPOINT (MP), COOR (2, MP), LINKP (2, MP)
C
      LOGICAL ADDLNK
C
      ADDLNK = .FALSE.
      DMIN = 100000.
C
      DO 100 I = 1, N15
         CALL LTSORT (MP, LINKP, I, IPNTR, ADDLNK)
         IF (IPNTR .GT. 0) THEN
            DIST = SQRT (((COOR (1, IPNTR) - X) **2) +
     &         ((COOR (2, IPNTR) - Y) **2))
            IF (DIST .LT. DMIN) THEN
               DMIN = DIST
               JJ = IPOINT (IPNTR)
            ENDIF
         ENDIF
  100 CONTINUE
      RETURN
      END
