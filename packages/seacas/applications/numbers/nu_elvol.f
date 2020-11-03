C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE ELVOL (NDIM, CRD, DISP, IX, NUMNP, NUMEL, NNODE,
     &     VOLUME)

C     ... ESTIMATE TIMESTEP FOR MESH --- BRICKS ONLY

      DIMENSION CRD(NUMNP, *), DISP(numnp, *), IX(NNODE,*)
      DIMENSION GRADOP(8,8)
      REAL volume(*)

      IF (NDIM .EQ. 3) THEN
         DO 20 IEL = 1, numel
            x1 = crd(ix(1,iel),1) + disp(ix(1,iel),1)
            x2 = crd(ix(2,iel),1) + disp(ix(2,iel),1)
            x3 = crd(ix(3,iel),1) + disp(ix(3,iel),1)
            x4 = crd(ix(4,iel),1) + disp(ix(4,iel),1)
            x5 = crd(ix(5,iel),1) + disp(ix(5,iel),1)
            x6 = crd(ix(6,iel),1) + disp(ix(6,iel),1)
            x7 = crd(ix(7,iel),1) + disp(ix(7,iel),1)
            x8 = crd(ix(8,iel),1) + disp(ix(8,iel),1)

            y1 = crd(ix(1,iel),2) + disp(ix(1,iel),2)
            y2 = crd(ix(2,iel),2) + disp(ix(2,iel),2)
            y3 = crd(ix(3,iel),2) + disp(ix(3,iel),2)
            y4 = crd(ix(4,iel),2) + disp(ix(4,iel),2)
            y5 = crd(ix(5,iel),2) + disp(ix(5,iel),2)
            y6 = crd(ix(6,iel),2) + disp(ix(6,iel),2)
            y7 = crd(ix(7,iel),2) + disp(ix(7,iel),2)
            y8 = crd(ix(8,iel),2) + disp(ix(8,iel),2)

            Z1 = crd(ix(1,iel),3) + disp(ix(1,iel),3)
            Z2 = crd(ix(2,iel),3) + disp(ix(2,iel),3)
            Z3 = crd(ix(3,iel),3) + disp(ix(3,iel),3)
            Z4 = crd(ix(4,iel),3) + disp(ix(4,iel),3)
            Z5 = crd(ix(5,iel),3) + disp(ix(5,iel),3)
            Z6 = crd(ix(6,iel),3) + disp(ix(6,iel),3)
            Z7 = crd(ix(7,iel),3) + disp(ix(7,iel),3)
            Z8 = crd(ix(8,iel),3) + disp(ix(8,iel),3)

            Z24 = Z2 - Z4
            Z52 = Z5 - Z2
            Z45 = Z4 - Z5
            GRADOP(1,1) = ( Y2*(Z6-Z3-Z45) + Y3*Z24 + Y4*(Z3-Z8-Z52)
     *           + Y5*(Z8-Z6-Z24) + Y6*Z52 + Y8*Z45 ) / 12.
            Z31 = Z3 - Z1
            Z63 = Z6 - Z3
            Z16 = Z1 - Z6
            GRADOP(2,1) = ( Y3*(Z7-Z4-Z16) + Y4*Z31 + Y1*(Z4-Z5-Z63)
     *           + Y6*(Z5-Z7-Z31) + Y7*Z63 + Y5*Z16 ) / 12.
            Z42 = Z4 - Z2
            Z74 = Z7 - Z4
            Z27 = Z2 - Z7
            GRADOP(3,1) = ( Y4*(Z8-Z1-Z27) + Y1*Z42 + Y2*(Z1-Z6-Z74)
     *           + Y7*(Z6-Z8-Z42) + Y8*Z74 + Y6*Z27 ) / 12.
            Z13 = Z1 - Z3
            Z81 = Z8 - Z1
            Z38 = Z3 - Z8
            GRADOP(4,1) = ( Y1*(Z5-Z2-Z38) + Y2*Z13 + Y3*(Z2-Z7-Z81)
     *           + Y8*(Z7-Z5-Z13) + Y5*Z81 + Y7*Z38 ) / 12.
            Z86 = Z8 - Z6
            Z18 = Z1 - Z8
            Z61 = Z6 - Z1
            GRADOP(5,1) = ( Y8*(Z4-Z7-Z61) + Y7*Z86 + Y6*(Z7-Z2-Z18)
     *           + Y1*(Z2-Z4-Z86) + Y4*Z18 + Y2*Z61 ) / 12.
            Z57 = Z5 - Z7
            Z25 = Z2 - Z5
            Z72 = Z7 - Z2
            GRADOP(6,1) = ( Y5*(Z1-Z8-Z72) + Y8*Z57 + Y7*(Z8-Z3-Z25)
     *           + Y2*(Z3-Z1-Z57) + Y1*Z25 + Y3*Z72 ) / 12.
            Z68 = Z6 - Z8
            Z36 = Z3 - Z6
            Z83 = Z8 - Z3
            GRADOP(7,1) = ( Y6*(Z2-Z5-Z83) + Y5*Z68 + Y8*(Z5-Z4-Z36)
     *           + Y3*(Z4-Z2-Z68) + Y2*Z36 + Y4*Z83 ) / 12.
            Z75 = Z7 - Z5
            Z47 = Z4 - Z7
            Z54 = Z5 - Z4
            GRADOP(8,1) = ( Y7*(Z3-Z6-Z54) + Y6*Z75 + Y5*(Z6-Z1-Z47)
     *           + Y4*(Z1-Z3-Z75) + Y3*Z47 + Y1*Z54 ) / 12.
            X24 = X2 - X4
            X52 = X5 - X2
            X45 = X4 - X5
            GRADOP(1,2) = ( Z2*(X6-X3-X45) + Z3*X24 + Z4*(X3-X8-X52)
     *           + Z5*(X8-X6-X24) + Z6*X52 + Z8*X45 ) / 12.
            X31 = X3 - X1
            X63 = X6 - X3
            X16 = X1 - X6
            GRADOP(2,2) = ( Z3*(X7-X4-X16) + Z4*X31 + Z1*(X4-X5-X63)
     *           + Z6*(X5-X7-X31) + Z7*X63 + Z5*X16 ) / 12.
            X42 = X4 - X2
            X74 = X7 - X4
            X27 = X2 - X7
            GRADOP(3,2) = ( Z4*(X8-X1-X27) + Z1*X42 + Z2*(X1-X6-X74)
     *           + Z7*(X6-X8-X42) + Z8*X74 + Z6*X27 ) / 12.
            X13 = X1 - X3
            X81 = X8 - X1
            X38 = X3 - X8
            GRADOP(4,2) = ( Z1*(X5-X2-X38) + Z2*X13 + Z3*(X2-X7-X81)
     *           + Z8*(X7-X5-X13) + Z5*X81 + Z7*X38 ) / 12.
            X86 = X8 - X6
            X18 = X1 - X8
            X61 = X6 - X1
            GRADOP(5,2) = ( Z8*(X4-X7-X61) + Z7*X86 + Z6*(X7-X2-X18)
     *           + Z1*(X2-X4-X86) + Z4*X18 + Z2*X61 ) / 12.
            X57 = X5 - X7
            X25 = X2 - X5
            X72 = X7 - X2
            GRADOP(6,2) = ( Z5*(X1-X8-X72) + Z8*X57 + Z7*(X8-X3-X25)
     *           + Z2*(X3-X1-X57) + Z1*X25 + Z3*X72 ) / 12.
            X68 = X6 - X8
            X36 = X3 - X6
            X83 = X8 - X3
            GRADOP(7,2) = ( Z6*(X2-X5-X83) + Z5*X68 + Z8*(X5-X4-X36)
     *           + Z3*(X4-X2-X68) + Z2*X36 + Z4*X83 ) / 12.
            X75 = X7 - X5
            X47 = X4 - X7
            X54 = X5 - X4
            GRADOP(8,2) = ( Z7*(X3-X6-X54) + Z6*X75 + Z5*(X6-X1-X47)
     *           + Z4*(X1-X3-X75) + Z3*X47 + Z1*X54 ) / 12.
            Y24 = Y2 - Y4
            Y52 = Y5 - Y2
            Y45 = Y4 - Y5
            GRADOP(1,3) = ( X2*(Y6-Y3-Y45) + X3*Y24 + X4*(Y3-Y8-Y52)
     *           + X5*(Y8-Y6-Y24) + X6*Y52 + X8*Y45 ) / 12.
            Y31 = Y3 - Y1
            Y63 = Y6 - Y3
            Y16 = Y1 - Y6
            GRADOP(2,3) = ( X3*(Y7-Y4-Y16) + X4*Y31 + X1*(Y4-Y5-Y63)
     *           + X6*(Y5-Y7-Y31) + X7*Y63 + X5*Y16 ) / 12.
            Y42 = Y4 - Y2
            Y74 = Y7 - Y4
            Y27 = Y2 - Y7
            GRADOP(3,3) = ( X4*(Y8-Y1-Y27) + X1*Y42 + X2*(Y1-Y6-Y74)
     *           + X7*(Y6-Y8-Y42) + X8*Y74 + X6*Y27 ) / 12.
            Y13 = Y1 - Y3
            Y81 = Y8 - Y1
            Y38 = Y3 - Y8
            GRADOP(4,3) = ( X1*(Y5-Y2-Y38) + X2*Y13 + X3*(Y2-Y7-Y81)
     *           + X8*(Y7-Y5-Y13) + X5*Y81 + X7*Y38 ) / 12.
            Y86 = Y8 - Y6
            Y18 = Y1 - Y8
            Y61 = Y6 - Y1
            GRADOP(5,3) = ( X8*(Y4-Y7-Y61) + X7*Y86 + X6*(Y7-Y2-Y18)
     *           + X1*(Y2-Y4-Y86) + X4*Y18 + X2*Y61 ) / 12.
            Y57 = Y5 - Y7
            Y25 = Y2 - Y5
            Y72 = Y7 - Y2
            GRADOP(6,3) = ( X5*(Y1-Y8-Y72) + X8*Y57 + X7*(Y8-Y3-Y25)
     *           + X2*(Y3-Y1-Y57) + X1*Y25 + X3*Y72 ) / 12.
            Y68 = Y6 - Y8
            Y36 = Y3 - Y6
            Y83 = Y8 - Y3
            GRADOP(7,3) = ( X6*(Y2-Y5-Y83) + X5*Y68 + X8*(Y5-Y4-Y36)
     *           + X3*(Y4-Y2-Y68) + X2*Y36 + X4*Y83 ) / 12.
            Y75 = Y7 - Y5
            Y47 = Y4 - Y7
            Y54 = Y5 - Y4
            GRADOP(8,3) = ( X7*(Y3-Y6-Y54) + X6*Y75 + X5*(Y6-Y1-Y47)
     *           + X4*(Y1-Y3-Y75) + X3*Y47 + X1*Y54 ) / 12.

C     Calculate element volume and characteristic element aspect ratio
C     (used in time step and hourglass control) -

            VOLUME(iel) = crd(ix(1,iel),1) * GRADOP(1,1)
     *           + crd(ix(2,iel),1) * GRADOP(2,1)
     *           + crd(ix(3,iel),1) * GRADOP(3,1)
     *           + crd(ix(4,iel),1) * GRADOP(4,1)
     *           + crd(ix(5,iel),1) * GRADOP(5,1)
     *           + crd(ix(6,iel),1) * GRADOP(6,1)
     *           + crd(ix(7,iel),1) * GRADOP(7,1)
     *           + crd(ix(8,iel),1) * GRADOP(8,1)
 20      CONTINUE
      ELSE
         STOP 'Not Implemented'
      END IF
      RETURN
      END
