C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C $Log: ptsnrm.f,v $
C Revision 1.2  2009/03/25 12:36:46  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.1  1994/04/07 20:08:16  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:55:30  gdsjaar
c Added RCS Id and Log to all files
c
C=======================================================================
      SUBROUTINE PTSNRM( PLAPTS, INPT,  CUTPT, CUTNRM , IERR)
C=======================================================================

C   --*** PTRNRM *** (MESH) Turns three point plane to point-normal form
C   --   Written by Ray J. Meyers 29 May, 1990
C   --
C   -- Given three non-colinear points, returns a description of the plane
C   -- consisting of a point on the plane and a normal vector to the plane.
C   --
C   --Parameters:
C   --   PLAPTS - IN - the 3 points defining the cutting plane;
C   --   INPT - IN - a point in the mesh, ie, on the opposite side
C                    from the cutting plane normal
C   --   CUTPT - OUT - a point on the cutting plane;
C   --   CUTNRM - OUT - the normal to the cutting plane;
C   --   IERR - OUT - non-zero if the input points are colinear (no plane
C                     defined)

      REAL PLAPTS(3,3)
      REAL INPT(3)
      REAL CUTPT(3)
      REAL CUTNRM(3)
      INTEGER IERR
      REAL VEC1(3), VEC2(3)
      REAL TOL, DOT
C
C DEFINE TWO VECTORS IN THE PLANE : VEC1 = P1-P2, VEC2 = P3-P2
C
      VEC1(1) = PLAPTS(1,1) - PLAPTS(1,2)
      VEC1(2) = PLAPTS(2,1) - PLAPTS(2,2)
      VEC1(3) = PLAPTS(3,1) - PLAPTS(3,2)
      VEC2(1) = PLAPTS(1,3) - PLAPTS(1,2)
      VEC2(2) = PLAPTS(2,3) - PLAPTS(2,2)
      VEC2(3) = PLAPTS(3,3) - PLAPTS(3,2)
C
C DEFINE A TOLERANCE BASED ON THE LENGTH OF  VEC1
C
      TOL = 1E-06 * AMAX1 ( VEC1(1), VEC1(2), VEC1(3) )
C
C TAKE THE CROSS PRODUCT OF VEC1 AND VEC2 AS THE NORMAL OF THE PLANE AND
C NORMALIZE
C
      CUTNRM(1) = VEC1(2)*VEC2(3) - VEC1(3)*VEC2(2)
      CUTNRM(2) = VEC1(3)*VEC2(1) - VEC1(1)*VEC2(3)
      CUTNRM(3) = VEC1(1)*VEC2(2) - VEC1(2)*VEC2(1)

      DIST = SQRT( CUTNRM(1)*CUTNRM(1) + CUTNRM(2)*CUTNRM(2) +
     $             CUTNRM(3)*CUTNRM(3) )
      CUTNRM(1) = CUTNRM(1)/DIST
      CUTNRM(2) = CUTNRM(2)/DIST
      CUTNRM(3) = CUTNRM(3)/DIST

C
C IF THE NORMAL IS (0,0,0), THE ORIGINAL POINTS WERE COLINEAR
C
      IF( ABS(CUTNRM(1)) .LT. TOL .AND. ABS(CUTNRM(2)) .LT. TOL
     $    .AND. ABS(CUTNRM(3)) .LT. TOL) THEN
          CALL PRTERR ('CMDERR', 'Points do not define a plane')
          IERR = 1
      ELSE
          IERR = 0
C
C USE THE SECOND INPUT POINT AS THE CHOSEN CUTPT
C
          CUTPT(1) = PLAPTS(1,2)
          CUTPT(2) = PLAPTS(2,2)
          CUTPT(3) = PLAPTS(3,2)
C
C CHECK THE INPT TO SEE IF THE NORMAL IS THE CORRECT SENSE
C
          VEC1(1) = INPT(1) - CUTPT(1)
          VEC1(2) = INPT(2) - CUTPT(2)
          VEC1(3) = INPT(3) - CUTPT(3)
          DOT = VEC1(1)*CUTNRM(1) + VEC1(2)*CUTNRM(2) +
     $          VEC1(3)*CUTNRM(3)
          IF(DOT .GT. 0) THEN
             CUTNRM(1) = -CUTNRM(1)
             CUTNRM(2) = -CUTNRM(2)
             CUTNRM(3) = -CUTNRM(3)
          END IF
      END IF

      RETURN
      END
