C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE DVOL3D( COORD, DISP, LSTSN, NSEG, DELVOL, NDIM, NUMNP)

C***********************************************************************

C     DESCRIPTION:
C       This routine computes the change in volume of a cavity formed
C       by the boundary of an element side set flag

C     FORMAL PARAMETERS:
C       COORD   REAL      Nodal Coordinates
C       DISP    REAL      Nodal Displacements
C       LSTSN   INTEGER   List of nodes on this boundary
C       LSTLEN  INTEGER   Length of node list
C       NSEG    INTEGER   Number of segments in this boundary
C       DELVOL  REAL      Change in volume of this cavity
C       NDIM   INTEGER   Number of Nodes

C     CALLED BY:

C***********************************************************************

      DIMENSION COORD(NUMNP, *), DISP(NUMNP, *), LSTSN(*)

      DELVOL = 0.0

      DO 100 KSEG = 1 , NSEG
         L = LSTSN(4*KSEG)
         K = LSTSN(4*KSEG - 1)
         J = LSTSN(4*KSEG - 2)
         I = LSTSN(4*KSEG - 3)

         Y1 = COORD(I,2)
         Y2 = COORD(J,2)
         Y3 = COORD(K,2)
         Y4 = COORD(L,2)

         Y5 = COORD(I,2) + DISP(I,2)
         Y6 = COORD(J,2) + DISP(J,2)
         Y7 = COORD(K,2) + DISP(K,2)
         Y8 = COORD(L,2) + DISP(L,2)

         Z1 = COORD(I,3)
         Z2 = COORD(J,3)
         Z3 = COORD(K,3)
         Z4 = COORD(L,3)

         Z5 = COORD(I,3) + DISP(I,3)
         Z6 = COORD(J,3) + DISP(J,3)
         Z7 = COORD(K,3) + DISP(K,3)
         Z8 = COORD(L,3) + DISP(L,3)

         Z24 = Z2 - Z4
         Z52 = Z5 - Z2
         Z45 = Z4 - Z5
         B11 = ( Y2*(Z6-Z3-Z45) + Y3*Z24 + Y4*(Z3-Z8-Z52)
     *         + Y5*(Z8-Z6-Z24) + Y6*Z52 + Y8*Z45 ) / 12.
         Z31 = Z3 - Z1
         Z63 = Z6 - Z3
         Z16 = Z1 - Z6
         B21 = ( Y3*(Z7-Z4-Z16) + Y4*Z31 + Y1*(Z4-Z5-Z63)
     *         + Y6*(Z5-Z7-Z31) + Y7*Z63 + Y5*Z16 ) / 12.
         Z42 = Z4 - Z2
         Z74 = Z7 - Z4
         Z27 = Z2 - Z7
         B31 = ( Y4*(Z8-Z1-Z27) + Y1*Z42 + Y2*(Z1-Z6-Z74)
     *         + Y7*(Z6-Z8-Z42) + Y8*Z74 + Y6*Z27 ) / 12.
         Z13 = Z1 - Z3
         Z81 = Z8 - Z1
         Z38 = Z3 - Z8
         B41 = ( Y1*(Z5-Z2-Z38) + Y2*Z13 + Y3*(Z2-Z7-Z81)
     *         + Y8*(Z7-Z5-Z13) + Y5*Z81 + Y7*Z38 ) / 12.
         Z86 = Z8 - Z6
         Z18 = Z1 - Z8
         Z61 = Z6 - Z1
         B51 = ( Y8*(Z4-Z7-Z61) + Y7*Z86 + Y6*(Z7-Z2-Z18)
     *         + Y1*(Z2-Z4-Z86) + Y4*Z18 + Y2*Z61 ) / 12.
         Z57 = Z5 - Z7
         Z25 = Z2 - Z5
         Z72 = Z7 - Z2
         B61 = ( Y5*(Z1-Z8-Z72) + Y8*Z57 + Y7*(Z8-Z3-Z25)
     *         + Y2*(Z3-Z1-Z57) + Y1*Z25 + Y3*Z72 ) / 12.
         Z68 = Z6 - Z8
         Z36 = Z3 - Z6
         Z83 = Z8 - Z3
         B71 = ( Y6*(Z2-Z5-Z83) + Y5*Z68 + Y8*(Z5-Z4-Z36)
     *         + Y3*(Z4-Z2-Z68) + Y2*Z36 + Y4*Z83 ) / 12.
         Z75 = Z7 - Z5
         Z47 = Z4 - Z7
         Z54 = Z5 - Z4
         B81 = ( Y7*(Z3-Z6-Z54) + Y6*Z75 + Y5*(Z6-Z1-Z47)
     *         + Y4*(Z1-Z3-Z75) + Y3*Z47 + Y1*Z54 ) / 12.

C Calculate volume of displaced element face

         VOL=  COORD(I,1)              * B11
     *      +  COORD(J,1)              * B21
     *      +  COORD(K,1)              * B31
     *      +  COORD(L,1)              * B41
     *      + (COORD(I,1) + DISP(I,1)) * B51
     *      + (COORD(J,1) + DISP(J,1)) * B61
     *      + (COORD(K,1) + DISP(K,1)) * B71
     *      + (COORD(L,1) + DISP(L,1)) * B81
         DELVOL = DELVOL + VOL
  100 CONTINUE

      RETURN
      END
