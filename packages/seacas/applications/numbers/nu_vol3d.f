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

C $Id: vol3d.f,v 1.2 1992/12/11 22:34:16 gdsjaar Exp $
C $Log: vol3d.f,v $
C Revision 1.2  1992/12/11 22:34:16  gdsjaar
C Fixed problem with incorrect determination of cavity center in 2d
C
c Revision 1.1.1.1  1991/02/21  15:46:20  gdsjaar
c NUMBERS: Greg Sjaardema, initial Unix release
c
c Revision 1.1  1991/02/21  15:46:19  gdsjaar
c Initial revision
c
      SUBROUTINE VOL3D( COORD, LSTSN, NSEG, VOLUME, NDIM, NUMESS,
     *   CENT, NUMNP, CENTER)
C
C***********************************************************************
C
C     DESCRIPTION:
C       This routine computes the volume of a cavity formed
C       by the boundary of an element side set flag
C
C     FORMAL PARAMETERS:
C       COORD   REAL      Nodal Coordinates
C       LSTSN   INTEGER   List of nodes on this boundary
C       NSEG    INTEGER   Number of segments in this boundary
C       VOLUME  REAL      Volume of this cavity
C       NDIM    INTEGER   Number of Nodes
C       CENT    REAL      Apex of Cavity Volume pentahedra
C
C     CALLED BY: 
C
C***********************************************************************
C
      DIMENSION COORD(NUMNP, *), LSTSN(*), CENT(3)
      LOGICAL CENTER
C
      VOLUME = 0.0
C
      X5 = CENT(1)
      Y5 = CENT(2)
      Z5 = CENT(3)
C
      DO 100 KSEG = 1 , NSEG
         L = LSTSN(4*KSEG)
         K = LSTSN(4*KSEG - 1)
         J = LSTSN(4*KSEG - 2)
         I = LSTSN(4*KSEG - 3)
C
         X1 = COORD(I,1) 
         X2 = COORD(J,1) 
         X3 = COORD(K,1) 
         X4 = COORD(L,1) 
C             
         Y1 = COORD(I,2) 
         Y2 = COORD(J,2) 
         Y3 = COORD(K,2) 
         Y4 = COORD(L,2) 
C             
         Z1 = COORD(I,3) 
         Z2 = COORD(J,3) 
         Z3 = COORD(K,3) 
         Z4 = COORD(L,3) 
C             
         Z13 = Z1 - Z3
         Z24 = Z2 - Z4
         Z31 = Z3 - Z1
         Z42 = Z4 - Z2
         Z51 = Z5 - Z1
         Z52 = Z5 - Z2
         Z53 = Z5 - Z3
         Z54 = Z5 - Z4
C
C         VP =(( Y3 * Z24 - Y2 * (Z3 + Z4) + Y4 * (Z2 + Z3) ) * X1 +
C     *        ( Y4 * Z31 - Y3 * (Z1 + Z4) + Y1 * (Z3 + Z4) ) * X2 + 
C     *        ( Y1 * Z42 - Y4 * (Z1 + Z2) + Y2 * (Z1 + Z4) ) * X3 +
C     *        ( Y2 * Z13 - Y1 * (Z2 + Z3) + Y3 * (Z1 + Z2) ) * X4)/12.0
         VP = ((2.*Y5 - Y3) * Z42 + Y2 * (Z53 + Z54) -
     *      Y4 * (Z53 + Z52) ) * X1 +
     *       ( (Y4 - 2.*Y5) * Z31 + Y3 * (Z54 + Z51) -
     *      Y1 * (Z54 + Z53) ) * X2 +
     *       ( (Y1 - 2.*Y5) * Z42 + Y4 * (Z51 + Z52) -
     *      Y2 * (Z54 + Z51) ) * X3 +
     *       ( (2.*Y5 - Y2) * Z31 + Y1 * (Z52 + Z53) -
     *      Y3 * (Z52 + Z51) ) * X4 +
     *      ((Y2 - Y4) * (Z3 - Z1) + (Y3 - Y1) *(Z4 - Z2)) * 2. * X5
         VOLUME = VOLUME + VP / 12.0
  100 CONTINUE
C
      RETURN
      END
