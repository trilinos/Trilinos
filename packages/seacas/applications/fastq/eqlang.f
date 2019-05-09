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

C $Id: eqlang.f,v 1.2 1991/03/21 15:44:42 gdsjaar Exp $
C $Log: eqlang.f,v $
C Revision 1.2  1991/03/21 15:44:42  gdsjaar
C Changed all 3.14159... to atan2(0.0, -1.0)
C
c Revision 1.1.1.1  1990/11/30  11:06:41  gdsjaar
c FASTQ Version 2.0X
c
c Revision 1.1  90/11/30  11:06:39  gdsjaar
c Initial revision
c
C
CC* FILE: [.PAVING]EQLANG.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE EQLANG (MXND, XN, YN, LXN, NODE, N0, N2, NFROM, DIST,
     &   VRO, XDEL, YDEL)
C***********************************************************************
C
C  SUBROUTINE EQLANG = CALCULATES A VECTOR SUM THAT ATTEMPTS TO
C                      MAINTAIN EQUAL ANGLES FOR A NODE
C
C***********************************************************************
C
      DIMENSION XN(MXND), YN(MXND), LXN(4, MXND)
C
      LOGICAL EXPAND
C
      PI = ATAN2(0.0, -1.0)
      TWOPI = 2.0 * PI

      IF (NFROM .GT. 0) THEN
C
C  TEST FOR THE EXPANSION CASE
C
         IF ( ( ((LXN (4, NFROM) .NE. 0) .AND.
     &      (LXN (2, NFROM) .LT. 0)) .OR.
     &      ((LXN (4, NFROM) .LT. 0) .AND.
     &      (LXN (2, NFROM) .GT. 0)) )
     &      .AND.
     &      ((LXN (3, N0) .EQ. 0) .OR. (LXN (3, N2) .EQ. 0)) ) THEN
            EXPAND = .TRUE.
         ELSE
            EXPAND = .FALSE.
         ENDIF
C
         ANG1 = ATAN2 ( YN (N2) - YN (NFROM), XN (N2) - XN (NFROM))
         IF (ANG1 .LT. 0.) ANG1 = ANG1 + TWOPI
         ANG2 = ATAN2 ( YN (N0) - YN (NFROM), XN (N0) - XN (NFROM))
         IF (ANG2 .LT. 0.) ANG2 = ANG2 + TWOPI
         ANG3 = ATAN2 ( YN (NODE) - YN (NFROM), XN (NODE) - XN (NFROM))
         IF (ANG3 .LT. 0.) ANG3 = ANG3 + TWOPI
C
C  GET THE APPROPRIATE ANGLE BETWEN ANGLE 1 AND 2
C
         ANG12D = ANG2 - ANG1
         IF (ANG12D .LT. 0.) ANG12D = ANG12D + TWOPI
C
C  IF THIS IS AN EXPANSION, THEN ADJUST THE ANGLE ACCORDINGLY
C
         IF (EXPAND) THEN
            IF (LXN (3, N2) .EQ. 0) THEN
               ANG12 = ANG1 + (ANG12D * .6)
            ELSEIF (LXN (3, N0) .EQ. 0) THEN
               ANG12 = ANG1 + (ANG12D * .4)
            ELSE
               ANG12 = ANG1 + (ANG12D * .5)
            ENDIF
         ELSE
            ANG12 = ANG1 + (ANG12D * .5)
         ENDIF
         IF (ANG12 .GT. TWOPI) ANG12 = ANG12 - TWOPI
C
C  GET THE AVERAGE ANGLE BETWEEN ANGLE 12 AND 3
C
         IF (ANG12 .GT. ANG3) THEN
            ANG3D = ANG12 - ANG3
            IF (ANG3D .GT. PI) THEN
               ANG = ANG12 + ((TWOPI - ANG3D) * .5)
            ELSE
               ANG = ANG12 - (ANG3D * .5)
            ENDIF
         ELSE
            ANG3D = ANG3 - ANG12
            IF (ANG3D .GT. PI) THEN
               ANG = ANG3 + ((TWOPI - ANG3D) * .5)
            ELSE
               ANG = ANG3 - (ANG3D * .5)
            ENDIF
         ENDIF
C
C  GET THE DISTANCE TO MAKE THE OUTSIDE FLAT AT THIS ANGLE
C
         D1 = SQRT ( ((XN (NFROM) - XN (N0)) ** 2) +
     &      ((YN (NFROM) - YN (N0)) ** 2) )
         D2 = SQRT ( ((XN (N2) - XN (N0)) ** 2) +
     &      ((YN (N2) - YN (N0)) ** 2) )
         D3 = SQRT ( ((XN (NFROM) - XN (N2)) ** 2) +
     &      ((YN (NFROM) - YN (N2)) ** 2) )
         ARG = (SIN (ANG12D) * D1) / D2
         IF (ARG .GT. 1.0) ARG = 1.0
         IF (ARG .LT. -1.0) ARG = -1.0
         BETA = ASIN (ARG)
         D0 = (D3 * SIN (BETA)) / SIN (PI - BETA - (ANG12D * .5))
C
         IF (D0 .GT. DIST) THEN
            IF (EXPAND) THEN
               DIST0 = D0
            ELSE
               DIST0 = (DIST + D0) * .5
            ENDIF
         ELSE
            DIST0 = DIST
         ENDIF
C
C  CALCULATE THE NEW COORDINATES
C
         X0 = XN (NFROM) + (COS (ANG) * DIST0)
         Y0 = YN (NFROM) + (SIN (ANG) * DIST0)
         XDEL = (X0 - XN (NODE)) * VRO
         YDEL = (Y0 - YN (NODE)) * VRO
C
      ELSE
         XDEL = 0.
         YDEL = 0.
      ENDIF
C
      RETURN
C
      END
