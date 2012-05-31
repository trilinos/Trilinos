C Copyright(C) 2009 Sandia Corporation. Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
C certain rights in this software.
C         
C Redistribution and use in source and binary forms, with or without
C modification, are permitted provided that the following conditions are
C met:
C 
C     * Redistributions of source code must retain the above copyright
C       notice, this list of conditions and the following disclaimer.
C 
C     * Redistributions in binary form must reproduce the above
C       copyright notice, this list of conditions and the following
C       disclaimer in the documentation and/or other materials provided
C       with the distribution.
C     * Neither the name of Sandia Corporation nor the names of its
C       contributors may be used to endorse or promote products derived
C       from this software without specific prior written permission.
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
