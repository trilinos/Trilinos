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

C $Id: dvol2d.f,v 1.1 1991/02/21 15:42:59 gdsjaar Exp $
C $Log: dvol2d.f,v $
C Revision 1.1  1991/02/21 15:42:59  gdsjaar
C Initial revision
C
      SUBROUTINE DVOL2D( COORD, DISP, LSTSN, NSEG, DELVOL,
     *   NDIM, AXI, NUMNP)
C
C***********************************************************************
C
C     DESCRIPTION:
C       This routine computes the change in volume of a cavity formed
C       by the boundary of an element side set flag
C
C     FORMAL PARAMETERS:
C       COORD   REAL      Nodal Coordinates
C       DISP    REAL      Nodal Displacements
C       LSTSN   INTEGER   List of nodes on this boundary
C       LSTLEN  INTEGER   Length of node list
C       NSEG    INTEGER   Number of segments in this boundary
C       DELVOL  REAL      Change in volume of this cavity
C       NUMNP   INTEGER   Number of Nodes
C       AXI     LOGICAL   TRUE if axisymmetric mesh
C
C     CALLED BY: 
C
C***********************************************************************
C
      DIMENSION COORD(NUMNP, *), DISP(NUMNP, *), LSTSN(*)
      LOGICAL AXI
      PI = ATAN2(0.0, -1.0)
C
      DELVOL = 0.0
C
      DO 100 KSEG = 1 , NSEG
         J = LSTSN(2*KSEG)
         I = LSTSN(2*KSEG - 1)
C
         X1  = COORD(I,1) 
         X2  = COORD(J,1) 
         DX1 = DISP(I,1)
         DX2 = DISP(J,1)
C
         Y1  = COORD(I,2) 
         Y2  = COORD(J,2) 
         DY1 = DISP(I,2)
         DY2 = DISP(J,2)
C             
         X12 = X1 - X2
         Y12 = Y1 - Y2
C
         VOL = (X12 * (DY2 + DY1) - Y12 * (DX2 + DX1) + DX1 * DY2 -
     *      DX2 * DY1 ) / 2.0
C
         IF (AXI) THEN
            XC = (DX2 + DX1) / 4.0 + (X2 + X1) / 2.0
            VOL = 2.0 * PI * XC * VOL
         END IF
C
         DELVOL = DELVOL + VOL
  100 CONTINUE
C
      RETURN
      END
