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

C $Id: vinter.f,v 1.1 1990/11/30 11:17:36 gdsjaar Exp $
C $Log: vinter.f,v $
C Revision 1.1  1990/11/30 11:17:36  gdsjaar
C Initial revision
C
C
CC* FILE: [.PAVING]VINTER.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE VINTER (MXND, XN, YN, N1, N2, N3, XOLD, YOLD,
     &   XNEW, YNEW, VCROSS)
C***********************************************************************
C
C  SUBROUTINE VINTER = FINDS WHERE A VECTOR FROM N1 TO N2
C                      INTERSECTS THE VECTOR FROM N3 TO (XOLD, YOLD)
C
C***********************************************************************
C
C  NOTE:  THIS INTERSECTION ROUTINE IS BASED ON AN ALGORITHM GIVEN
C         IN THE BOOK "GEOMETRIC MODELING" BY MICHAEL E. MORTENSON ON
C         PAGES 319 - 320.
C
C***********************************************************************
C
      DIMENSION XN (MXND), YN (MXND)
C
      LOGICAL VCROSS
C
      VCROSS = .FALSE.
C
C  SET UP THE FIRST LINE'S VECTORS (A AND B)
C
      XA = XN (N1)
      YA = YN (N1)
      XB = XN (N2) - XN (N1)
      YB = YN (N2) - YN (N1)
C
C  SET UP THE SECOND LINE'S VECTORS (C AND D)
C
      XC = XN (N3)
      YC = YN (N3)
      XD = XOLD - XN (N3)
      YD = YOLD - YN (N3)
C
C  NOW USE THE VECTORS AND SOLVE FOR W.
C  W IS THE PROPORTION OF THE DISTANCE ALONG THE VECTOR D
C  WHERE THE INTERSECTION OCCURS.  LIKEWISE U IS THE PROPORTIONAL
C  DISTANCE ALONG THE VECTOR B FOR THE INTERSECTION.
C
      DENOM = (YB * XD) - (XB * YD)
C
C  CHECK FOR SPECIAL PARALLEL CASE - THE DENOMINATOR IS EQUAL TO ZERO.
C
      IF (DENOM .NE. 0.) THEN
C
C  GET INTERSECTION LOCATION
C
         W = ( (YC * XB) - (XB * YA) - (XC * YB) + (YB * XA) ) / DENOM
C
C  GET THE U VALUE TO CONFIRM.
C
         IF (XB .NE. 0.) THEN
            U = ( XC + (W * XD) - XA ) / XB
         ELSE
            U = ( YC + (W * YD) - YA ) / YB
         ENDIF
C
C  CALCULATE THE INTERSECTION POINT BASED ON SIMILAR TRIANGLES
C
         XNEW = ( (XA + (XB * U)) + (XC + (XD * W)) ) * .5
         YNEW = ( (YA + (YB * U)) + (YC + (YD * W)) ) * .5
         VCROSS = .TRUE.
      ENDIF
C
      RETURN
C
      END
