C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C    
C    See packages/seacas/LICENSE for details

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
