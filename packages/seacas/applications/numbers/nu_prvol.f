C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE PRVOL (NDIM, CRD, IX, NUMNP, NUMEL, NNODE,
     &  VOLUME, IUNIT)

C     ... ESTIMATE TIMESTEP FOR MESH --- BRICKS ONLY

      DIMENSION CRD(NUMNP, *), IX(NNODE,*)
      DIMENSION GRADOP(8)
      REAL volume(*)

      IF (NDIM .EQ. 3 .AND. NNODE .EQ. 8) THEN
        DO 20 IEL = 1, numel
          y1 = crd(ix(1,iel),2)
          y2 = crd(ix(2,iel),2)
          y3 = crd(ix(3,iel),2)
          y4 = crd(ix(4,iel),2)
          y5 = crd(ix(5,iel),2)
          y6 = crd(ix(6,iel),2)
          y7 = crd(ix(7,iel),2)
          y8 = crd(ix(8,iel),2)

          Z1 = crd(ix(1,iel),3)
          Z2 = crd(ix(2,iel),3)
          Z3 = crd(ix(3,iel),3)
          Z4 = crd(ix(4,iel),3)
          Z5 = crd(ix(5,iel),3)
          Z6 = crd(ix(6,iel),3)
          Z7 = crd(ix(7,iel),3)
          Z8 = crd(ix(8,iel),3)

          Z24 = Z2 - Z4
          Z52 = Z5 - Z2
          Z45 = Z4 - Z5
          GRADOP(1) = ( Y2*(Z6-Z3-Z45) + Y3*Z24 + Y4*(Z3-Z8-Z52)
     *      + Y5*(Z8-Z6-Z24) + Y6*Z52 + Y8*Z45 ) / 12.
          Z31 = Z3 - Z1
          Z63 = Z6 - Z3
          Z16 = Z1 - Z6
          GRADOP(2) = ( Y3*(Z7-Z4-Z16) + Y4*Z31 + Y1*(Z4-Z5-Z63)
     *      + Y6*(Z5-Z7-Z31) + Y7*Z63 + Y5*Z16 ) / 12.
          Z42 = Z4 - Z2
          Z74 = Z7 - Z4
          Z27 = Z2 - Z7
          GRADOP(3) = ( Y4*(Z8-Z1-Z27) + Y1*Z42 + Y2*(Z1-Z6-Z74)
     *      + Y7*(Z6-Z8-Z42) + Y8*Z74 + Y6*Z27 ) / 12.
          Z13 = Z1 - Z3
          Z81 = Z8 - Z1
          Z38 = Z3 - Z8
          GRADOP(4) = ( Y1*(Z5-Z2-Z38) + Y2*Z13 + Y3*(Z2-Z7-Z81)
     *      + Y8*(Z7-Z5-Z13) + Y5*Z81 + Y7*Z38 ) / 12.
          Z86 = Z8 - Z6
          Z18 = Z1 - Z8
          Z61 = Z6 - Z1
          GRADOP(5) = ( Y8*(Z4-Z7-Z61) + Y7*Z86 + Y6*(Z7-Z2-Z18)
     *      + Y1*(Z2-Z4-Z86) + Y4*Z18 + Y2*Z61 ) / 12.
          Z57 = Z5 - Z7
          Z25 = Z2 - Z5
          Z72 = Z7 - Z2
          GRADOP(6) = ( Y5*(Z1-Z8-Z72) + Y8*Z57 + Y7*(Z8-Z3-Z25)
     *      + Y2*(Z3-Z1-Z57) + Y1*Z25 + Y3*Z72 ) / 12.
          Z68 = Z6 - Z8
          Z36 = Z3 - Z6
          Z83 = Z8 - Z3
          GRADOP(7) = ( Y6*(Z2-Z5-Z83) + Y5*Z68 + Y8*(Z5-Z4-Z36)
     *      + Y3*(Z4-Z2-Z68) + Y2*Z36 + Y4*Z83 ) / 12.
          Z75 = Z7 - Z5
          Z47 = Z4 - Z7
          Z54 = Z5 - Z4
          GRADOP(8) = ( Y7*(Z3-Z6-Z54) + Y6*Z75 + Y5*(Z6-Z1-Z47)
     *      + Y4*(Z1-Z3-Z75) + Y3*Z47 + Y1*Z54 ) / 12.

C     Calculate element volume

          VOLUME(iel) = crd(ix(1,iel),1) * GRADOP(1)
     *      + crd(ix(2,iel),1) * GRADOP(2)
     *      + crd(ix(3,iel),1) * GRADOP(3)
     *      + crd(ix(4,iel),1) * GRADOP(4)
     *      + crd(ix(5,iel),1) * GRADOP(5)
     *      + crd(ix(6,iel),1) * GRADOP(6)
     *      + crd(ix(7,iel),1) * GRADOP(7)
     *      + crd(ix(8,iel),1) * GRADOP(8)
          if (volume(iel) .lt. 0.0) then
            write (*,*) 'Zero or negative volume at element',
     &        iel
          endif
 20     CONTINUE

C ... Special for Frank Dempsey -- Print volumes and connectivity to file
        DO 30 IEL=1, NUMEL
          write (IUNIT,999) iel, volume(iel), (ix(i,iel),i=1,8)
 999      format(i8,1PE17.8,8i8)
 30     CONTINUE
      ELSE
        STOP 'Not Implemented'
      END IF
      RETURN
      END

