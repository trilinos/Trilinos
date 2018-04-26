C Copyright(C) 2011-2017 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C Redistribution and use in source and binary forms, with or without
C modification, are permitted provided that the following conditions are
C met:
C 
C * Redistributions of source code must retain the above copyright
C    notice, this list of conditions and the following disclaimer.
C           
C * Redistributions in binary form must reproduce the above
C   copyright notice, this list of conditions and the following
C   disclaimer in the documentation and/or other materials provided
C   with the distribution.
C                         
C * Neither the name of NTESS nor the names of its
C   contributors may be used to endorse or promote products derived
C   from this software without specific prior written permission.
C                                                 
C THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
C "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
C LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
C A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
C OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
C SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
C LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
C DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
C THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
C (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
C OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

C=======================================================================
      SUBROUTINE WRPXYZ (XN, YN, XN3, YN3, ZN3, ATRIB)
C=======================================================================
 
C   $Id: wrpxyz.f,v 1.5 1993/05/27 22:17:06 gdsjaar Exp $
C   $Log: wrpxyz.f,v $
C   Revision 1.5  1993/05/27 22:17:06  gdsjaar
C   Added new ellipse transformation code
C
c Revision 1.4  1991/07/31  17:30:53  gdsjaar
c Added WARP AXIS VERTICAL command to map to surface
c without changing input X and Y coordinates.  Updated Version to X0.01.00
c
c Revision 1.3  1991/01/11  08:39:58  gdsjaar
c Removed DEBUG comment lines
c
c Revision 1.2  91/01/09  12:59:47  gdsjaar
c Initial conversion from GEN3D to GENSHELL, no BC yet
c
c Revision 1.1.1.1  90/08/20  12:23:33  gdsjaar
c Gen3D Mesh Generation Program
c
c Revision 1.1  90/08/20  12:23:31  gdsjaar
c Initial revision
c
 
C   --*** WRPXYZ *** (GEN3D) Calculate 3D coordinates
C   --   Written by Amy Gilkey - revised 05/09/88
C   --   Modified by Greg Sjaardema - 02/06/89
C   --       Added Warp Function
C   --       Added Gradient to Rotations (not for center blocks)
C   --
C   --WRPXYZ calculates the coordinate array for the 3D warp translations.
C   --
C   --Parameters:
C   --   XN, YN - IN - the 2D coordinates, destroyed
C   --   XN3, YN3, ZN3 - OUT - the 3D coordinates
C   --
C   --Common Variables:
C   --   Uses NDIM, NUMNP of /DBNUMS/
C   --   Uses NDIM3, NUMNP3 of /DBNUM3/
C   --   Uses DOTRAN, NNREPL, DIM3, NRTRAN, D3TRAN, ZGRAD,
C   --      CENTER, NUMCOL, NUMROW of /PARAMS/
 
      INCLUDE 'gs_dbnums.blk'
      INCLUDE 'gs_dbnum3.blk'
      INCLUDE 'gs_params.blk'
 
      REAL XN(NUMNP), YN(NUMNP),
     &   XN3(NUMNP3), YN3(NUMNP3), ZN3(NUMNP3)
      REAL ATRIB(NUMEL)
 
      IF (IWARP .EQ. 1) THEN
C
C ... Warp type 1: Point Centered
C
         DO 100 INP = 1, NUMNP
            XN3(INP) = XN(INP)
            YN3(INP) = YN(INP)
            ZN3(INP) = DWARP - SQRT(DWARP**2 - XN(INP)**2 - YN(INP)**2)
  100    CONTINUE
 
         CONTINUE
      ELSE IF (IWARP .EQ. -1) THEN
C
C ... Warp type -1: X Axis Centered
C
         DO 110 INP = 1, NUMNP
            THET = YN(INP) / DWARP
            XN3(INP) = XN(INP)
            YN3(INP) = SIN(THET) * DWARP
            ZN3(INP) = DWARP - COS(THET) * DWARP
  110    CONTINUE
 
      ELSE IF (IWARP .EQ. -2) THEN
C
C ... Warp type -2: Y Axis Centered
C
         DO 120 INP = 1, NUMNP
            THET = XN(INP) / DWARP
            XN3(INP) = SIN(THET) * DWARP
            YN3(INP) = YN(INP)
            ZN3(INP) = DWARP - COS(THET) * DWARP
  120    CONTINUE
 
      ELSE IF (IWARP .EQ. -3) THEN
C
C ... Warp type -3: X Axis Centered, Project straight up
C
         DO 130 INP = 1, NUMNP
            XN3(INP) = XN(INP)
            YN3(INP) = YN(INP)
            ZN3(INP) = DWARP - sqrt( dwarp**2 - yn(inp)**2 )
  130    CONTINUE
 
      ELSE IF (IWARP .EQ. -4) THEN
C
C ... Warp type -4: Y Axis Centered, Project straight up
C
         DO 140 INP = 1, NUMNP
            XN3(INP) = XN(INP)
            YN3(INP) = YN(INP)
            ZN3(INP) = DWARP - sqrt( dwarp**2 - xn(inp)**2 )
  140    CONTINUE
 
      ELSE IF (IWARP .EQ. 2) THEN
C
C ... Warp type 2: Point-Centered Ellipse
C
         DO 150 INP = 1, NUMNP
            DX = XN(INP)
            DY = YN(INP)
            ZT = DWARP * (1.0-SQRT(1.0-DX**2/XRAD**2-DY**2/YRAD**2))
            ZN3(INP) = ZT
            YN3(INP) = YN(INP)
            XN3(INP) = XN(INP)
  150    CONTINUE
      END IF
 
C ... Now do the attributes for all of the elements
 
      DO 160 IEL = 1, NUMEL
         ATRIB(IEL) = DIM3
  160 CONTINUE
 
      RETURN
      END
