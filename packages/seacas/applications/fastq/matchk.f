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

C $Id: matchk.f,v 1.1 1990/11/30 11:11:57 gdsjaar Exp $
C $Log: matchk.f,v $
C Revision 1.1  1990/11/30 11:11:57  gdsjaar
C Initial revision
C
C
CC* FILE: [.PAVING]MATCHK.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      LOGICAL FUNCTION MATCHK (MXND, I1, I2, J1, J2, LXN)
C***********************************************************************
C
C  FUNCTION MATCHK = CHECKS THE CURRENT COLAPSED LINES TO SEE IF THEY
C                    CAN BE JOINED WITHOUT AFFECTING THE BOUNDARY.
C                    I1 & I2 MAY END UP SWITCHED WITH J1 & J2.
C
C***********************************************************************
C
      DIMENSION LXN (4, MXND)
C
      IF ( (LXN (2, I1) .LT. 0) .OR. (LXN (2, I2) .LT. 0) .OR.
     &   (LXN (2, J1) .LT. 0) .OR. (LXN (2, J2) .LT. 0) ) THEN
C
C  FIRST CHECK FOR COMPLETELY HOOKED BOUNDARY LINES.
C
         IF ((LXN (2, J1) .LT. 0) .AND. (LXN (2, J2) .LT. 0)) THEN
            MATCHK = .FALSE.
         ELSEIF ( ((LXN (2, I1) .LT. 0) .AND. (LXN (2, J2) .LT. 0)) .OR.
     &      ((LXN (2, I2) .LT. 0) .AND. (LXN (2, J1) .LT. 0)))
     &      THEN
            MATCHK = .FALSE.
         ELSE
            MATCHK = .TRUE.
         ENDIF
      ELSE
         MATCHK = .TRUE.
      ENDIF
C
      RETURN
C
      END
