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

C $Id: gettrn.f,v 1.1 1990/11/30 11:08:48 gdsjaar Exp $
C $Log: gettrn.f,v $
C Revision 1.1  1990/11/30 11:08:48  gdsjaar
C Initial revision
C
C
CC* FILE: [.QMESH]GETTRN.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE GETTRN (ML, MS, MNNPS, NS, ISLIST, NINT, IFLINE, NLPS,
     &   ILLIST, LINKL, LINKS, X, Y, NID, NNPS, ANGLE, NPER, I1, I2,
     &   I3, I4, I5, I6, I7, I8, XCEN1, YCEN1, XCEN2, YCEN2, XMID1,
     &   YMID1, XMID2, YMID2, CCW, HALFC, ERR)
C***********************************************************************
C
C  SUBROUTINE GETTRN = GETS THE APPROPRIATE SIDES,  CORNERS,  AND MIDPOINT
C                      VALUES FOR A TRANSITION REGION
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
      DIMENSION NNPS(MNNPS), ISLIST(NS), LINKL(2, ML), LINKS(MS * 2)
      DIMENSION NINT (ML), NLPS (MS), IFLINE (MS), ILLIST (MS * 3)
      DIMENSION X (NPER), Y (NPER), NID (NPER), ANGLE (NPER)
C
      LOGICAL CCW, ERR, HALFC
C
C  CALCULATE THE NUMBER OF NODES PER SIDE
C
      CALL NPS (ML, MS, MNNPS, NS, ISLIST, NINT, IFLINE, NLPS, ILLIST,
     &   LINKL, LINKS, NNPS, ERR)
      IF (ERR) RETURN
      IF (.NOT.CCW)CALL IREVER (NNPS, NS)
C
C  FIND THE BEST CORNER NODES IN THE LIST
C
      CALL PICKTR (NPER, X, Y, NID, ANGLE, HALFC, I1, I2, I3, I4, I5,
     &   I6, I7, I8)
C
C  DEFINE THE MIDDLE POINT OF BOTH TRIANGLES AS THE AVERAGE
C  OF PROPORIONAL DIVISIONS OF SIDE DIVISION POINT TO OPPOSITE
C  TRIANGLE CORNER LINES
C  FOR THE FIRST TRIANGLE,
C  FIND DISTANCES FROM CORNER TO CORNER,  AND CORNERS TO SIDE DIVISIONS
C
      INT = I6 - I4
      PROP = FLOAT (I6 - I5) / FLOAT (INT)
      XMID1 = X (I3) +  (PROP *  (X (I7) - X (I3)))
      YMID1 = Y (I3) +  (PROP *  (Y (I7) - Y (I3)))
      D1 = SQRT ( (X (I5) - X (I3)) **2 +  (Y (I5) - Y (I3)) **2)
      D2 = SQRT ( (X (I7) - X (I5)) **2 +  (Y (I7) - Y (I5)) **2)
      D3 = SQRT ( (X (I3) - X (I7)) **2 +  (Y (I3) - Y (I7)) **2)
      D1A = SQRT ( (X (I4) - X (I3)) **2 +  (Y (I4) - Y (I3)) **2)
      D1B = SQRT ( (X (I5) - X (I4)) **2 +  (Y (I5) - Y (I4)) **2)
      D2A = SQRT ( (X (I6) - X (I5)) **2 +  (Y (I6) - Y (I5)) **2)
      D2B = SQRT ( (X (I7) - X (I6)) **2 +  (Y (I7) - Y (I6)) **2)
      D3A = SQRT ( (XMID1 - X (I7)) **2 +  (YMID1 - Y (I7)) **2)
      D3B = SQRT ( (X (I3) - XMID1) **2 +  (Y (I3) - YMID1) **2)
C
C  GET MIDPOINT TRIALS 1,  2,  AND 3 AS PROPORTIONS
C
      PRO1 = .5 *  ( (D3A / D3) +  (D1B / D1))
      X1 = X (I6) -  (PRO1 *  (X (I6) - X (I3)))
      Y1 = Y (I6) -  (PRO1 *  (Y (I6) - Y (I3)))
      PRO2 = .5 *  ( (D2B / D2) +  (D1A / D1))
      X2 = XMID1 -  (PRO2 *  (XMID1 - X (I5)))
      Y2 = YMID1 -  (PRO2 *  (YMID1 - Y (I5)))
      PRO3 = .5 *  ( (D2A / D2) +  (D3B / D3))
      X3 = X (I4) -  (PRO3 *  (X (I4) - X (I7)))
      Y3 = Y (I4) -  (PRO3 *  (Y (I4) - Y (I7)))
C
C  AVERAGE POINTS TO GET THE FIRST CENTER
C
      XCEN1 =  (X1 + X2 + X3) / 3.
      YCEN1 =  (Y1 + Y2 + Y3) / 3.
C
C  FOR THE SECOND TRIANGLE,
C  FIND DISTANCES FROM CORNER TO CORNER,  AND CORNERS TO SIDE DIVISIONS
C
      INT = I6 - I4
      PROP = FLOAT (NPER + 1 - I8) / FLOAT (INT)
      XMID2 = X (I3) +  (PROP *  (X (I7) - X (I3)))
      YMID2 = Y (I3) +  (PROP *  (Y (I7) - Y (I3)))
      D1 = SQRT ( (X (I3) - X (I1)) **2 +  (Y (I3) - Y (I1)) **2)
      D2 = SQRT ( (X (I7) - X (I3)) **2 +  (Y (I7) - Y (I3)) **2)
      D3 = SQRT ( (X (I1) - X (I7)) **2 +  (Y (I1) - Y (I7)) **2)
      D1A = SQRT ( (X (I2) - X (I1)) **2 +  (Y (I2) - Y (I1)) **2)
      D1B = SQRT ( (X (I3) - X (I2)) **2 +  (Y (I3) - Y (I2)) **2)
      D2A = SQRT ( (XMID2 - X (I3)) **2 +  (YMID2 - Y (I3)) **2)
      D2B = SQRT ( (X (I7) - XMID2) **2 +  (Y (I7) - YMID2) **2)
      D3A = SQRT ( (X (I8) - X (I7)) **2 +  (Y (I8) - Y (I7)) **2)
      D3B = SQRT ( (X (I1) - X (I8)) **2 +  (Y (I1) - Y (I8)) **2)
C
C  GET MIDPOINT TRIALS 1,  2,  AND 3 AS PROPORTIONS
C
      PRO1 = .5 *  ((D3A / D3) +  (D1B / D1))
      X1 = XMID2 -  (PRO1 *  (XMID2 - X (I1)))
      Y1 = YMID2 -  (PRO1 *  (YMID2 - Y (I1)))
      PRO2 = .5 *  ((D2B / D2) +  (D1A / D1))
      X2 = X (I8) -  (PRO2 *  (X (I8) - X (I3)))
      Y2 = Y (I8) -  (PRO2 *  (Y (I8) - Y (I3)))
      PRO3 = .5 *  ((D2A / D2) +  (D3B / D3))
      X3 = X (I2) -  (PRO3 *  (X (I2) - X (I7)))
      Y3 = Y (I2) -  (PRO3 *  (Y (I2) - Y (I7)))
C
C  AVERAGE POINTS TO GET THE CENTER
C
      XCEN2 =  (X1 + X2 + X3) / 3.
      YCEN2 =  (Y1 + Y2 + Y3) / 3.
C
      ERR = .FALSE.
      RETURN
C
      END
