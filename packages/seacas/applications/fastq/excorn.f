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

C $Id: excorn.f,v 1.2 1991/03/21 15:44:44 gdsjaar Exp $
C $Log: excorn.f,v $
C Revision 1.2  1991/03/21 15:44:44  gdsjaar
C Changed all 3.14159... to atan2(0.0, -1.0)
C
c Revision 1.1.1.1  1990/11/30  11:06:58  gdsjaar
c FASTQ Version 2.0X
c
c Revision 1.1  90/11/30  11:06:57  gdsjaar
c Initial revision
c
C
CC* FILE: [.PAVING]EXCORN.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE EXCORN (MXND, XN, YN, LNODES, ANGLE, N0, N1, N2, XNEW,
     &   YNEW)
C***********************************************************************
C
C  SUBROUTINE EXCORN = CALCULATES A POSITION AN AVERAGE LENGTH AWAY
C                      FROM A CORNER NODE
C
C***********************************************************************
C
      DIMENSION XN (MXND), YN (MXND), LNODES (7, MXND), ANGLE (MXND)
C
      LOGICAL SIDEP
C
C      XNEW = XN (N0) + XN (N2) - XN (N1)
C      YNEW = YN (N0) + YN (N2) - YN (N1)

      PID2 = 0.5 * ATAN2(0.0, -1.0)
C
C      ANG2 = ATAN2 (YN (N1)-YN (N0), XN (N1)-XN (N0))+PID2
      BANG1 = ATAN2 (YN (N1) - YN (N2), XN (N1) - XN (N2))
      BANG2 = ATAN2 (YN (N1) - YN (N0), XN (N1) - XN (N0))
      IF (SIDEP (ANGLE (N2)))THEN
         ANG1 = BANG1 - (ANGLE (N2) * .5)
      ELSE
         ANG1 = BANG1 - PID2
      ENDIF
      IF (SIDEP (ANGLE (N0)))THEN
         ANG2 = BANG2 + (ANGLE (N0) * .5)
      ELSE
         ANG2 = BANG2 + PID2
      ENDIF
      DIST1 = SQRT ((YN (N2)-YN (N1)) **2 +  (XN (N2)-XN (N1)) **2)
      DIST2 = SQRT ((YN (N0)-YN (N1)) **2 +  (XN (N0)-XN (N1)) **2)
      DIST = (DIST1 + DIST2) * .5
      XNEW = ( (DIST * COS (ANG1) + XN (N2)) +
     &   (DIST * COS (ANG2) + XN (N0)) ) * .5
      YNEW = ( (DIST * SIN (ANG1) + YN (N2)) +
     &   (DIST * SIN (ANG2) + YN (N0)) ) * .5
C      XNEW = DIST * COS (ANG1) + XN (N2)
C      YNEW = DIST * SIN (ANG1) + YN (N2)
C
      RETURN
C
      END
