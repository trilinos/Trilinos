C    Copyright (c) 2014, Sandia Corporation.
C    Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
C    the U.S. Governement retains certain rights in this software.
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

C $Id: extnd1.f,v 1.2 1998/07/14 18:18:50 gdsjaar Exp $
C $Log: extnd1.f,v $
C Revision 1.2  1998/07/14 18:18:50  gdsjaar
C Removed unused variables, cleaned up a little.
C
C Changed BLUE labels to GREEN to help visibility on black background
C (indirectly requested by a couple users)
C
C Revision 1.1.1.1  1990/11/30 11:07:09  gdsjaar
C FASTQ Version 2.0X
C
c Revision 1.1  90/11/30  11:07:08  gdsjaar
c Initial revision
c 
C
CC* FILE: [.PAVING]EXTND1.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE EXTND1 (MXND, XN, YN, ANGLE, N1, N2, N3, X, Y, DIST)
C***********************************************************************
C
C  SUBROUTINE EXCORN = CALCULATES A POSITION AN AVERAGE LENGTH AWAY
C                      FROM A CORNER NODE
C
C***********************************************************************
C
      DIMENSION XN (MXND), YN (MXND), ANGLE (MXND)
      DIMENSION X(1), Y(1)
C
      CANG = (ANGLE (N2) * .5)
      ANG = ATAN2 (YN (N1) - YN (N2), XN (N1) - XN (N2)) - CANG
      DIST1 = SQRT ((YN (N2) - YN (N1)) **2 +  (XN (N2) - XN (N1)) **2)
      DIST2 = SQRT ((YN (N3) - YN (N2)) **2 +  (XN (N3) - XN (N2)) **2)
      DIST = (DIST1 + DIST2) * .5
      IF (CANG .EQ. 0.) THEN
         ADIST = DIST
      ELSE
         ADIST = DIST / SIN (CANG)
      ENDIF
C
      X(1) = ADIST * COS (ANG) + XN (N2)
      Y(1) = ADIST * SIN (ANG) + YN (N2)
C
      RETURN
C
      END
