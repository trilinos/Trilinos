C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE GETTRN (ML, MS, MNNPS, NS, ISLIST, NINT, IFLINE, NLPS,
     &   ILLIST, LINKL, LINKS, X, Y, NID, NNPS, ANGLE, NPER, I1, I2,
     &   I3, I4, I5, I6, I7, I8, XCEN1, YCEN1, XCEN2, YCEN2, XMID1,
     &   YMID1, XMID2, YMID2, CCW, HALFC, ERR)
C***********************************************************************

C  SUBROUTINE GETTRN = GETS THE APPROPRIATE SIDES,  CORNERS,  AND MIDPOINT
C                      VALUES FOR A TRANSITION REGION

C***********************************************************************

C  SUBROUTINE CALLED BY:
C     QMESH = GENERATES QUAD ELEMENTS

C***********************************************************************

C  VARIABLES USED:
C     NNPS  = ARRAY OF NUMBER OF NODES PER SIDE
C     CCW   = .TRUE. IF THE SIDE IS ORIENTED CCW
C     NORM  = .TRUE. IF THE FIRST SIDE IS TO BE TRIED AS THE BASE

C***********************************************************************

      DIMENSION NNPS(MNNPS), ISLIST(NS), LINKL(2, ML), LINKS(MS * 2)
      DIMENSION NINT (ML), NLPS (MS), IFLINE (MS), ILLIST (MS * 3)
      DIMENSION X (NPER), Y (NPER), NID (NPER), ANGLE (NPER)

      LOGICAL CCW, ERR, HALFC

C  CALCULATE THE NUMBER OF NODES PER SIDE

      CALL NPS (ML, MS, MNNPS, NS, ISLIST, NINT, IFLINE, NLPS, ILLIST,
     &   LINKL, LINKS, NNPS, ERR)
      IF (ERR) RETURN
      IF (.NOT.CCW)CALL IREVER (NNPS, NS)

C  FIND THE BEST CORNER NODES IN THE LIST

      CALL PICKTR (NPER, X, Y, NID, ANGLE, HALFC, I1, I2, I3, I4, I5,
     &   I6, I7, I8)

C  DEFINE THE MIDDLE POINT OF BOTH TRIANGLES AS THE AVERAGE
C  OF PROPORIONAL DIVISIONS OF SIDE DIVISION POINT TO OPPOSITE
C  TRIANGLE CORNER LINES
C  FOR THE FIRST TRIANGLE,
C  FIND DISTANCES FROM CORNER TO CORNER,  AND CORNERS TO SIDE DIVISIONS

      INT = I6 - I4
      PROP = DBLE(I6 - I5) / DBLE(INT)
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

C  GET MIDPOINT TRIALS 1,  2,  AND 3 AS PROPORTIONS

      PRO1 = .5 *  ( (D3A / D3) +  (D1B / D1))
      X1 = X (I6) -  (PRO1 *  (X (I6) - X (I3)))
      Y1 = Y (I6) -  (PRO1 *  (Y (I6) - Y (I3)))
      PRO2 = .5 *  ( (D2B / D2) +  (D1A / D1))
      X2 = XMID1 -  (PRO2 *  (XMID1 - X (I5)))
      Y2 = YMID1 -  (PRO2 *  (YMID1 - Y (I5)))
      PRO3 = .5 *  ( (D2A / D2) +  (D3B / D3))
      X3 = X (I4) -  (PRO3 *  (X (I4) - X (I7)))
      Y3 = Y (I4) -  (PRO3 *  (Y (I4) - Y (I7)))

C  AVERAGE POINTS TO GET THE FIRST CENTER

      XCEN1 =  (X1 + X2 + X3) / 3.
      YCEN1 =  (Y1 + Y2 + Y3) / 3.

C  FOR THE SECOND TRIANGLE,
C  FIND DISTANCES FROM CORNER TO CORNER,  AND CORNERS TO SIDE DIVISIONS

      INT = I6 - I4
      PROP = DBLE(NPER + 1 - I8) / DBLE(INT)
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

C  GET MIDPOINT TRIALS 1,  2,  AND 3 AS PROPORTIONS

      PRO1 = .5 *  ((D3A / D3) +  (D1B / D1))
      X1 = XMID2 -  (PRO1 *  (XMID2 - X (I1)))
      Y1 = YMID2 -  (PRO1 *  (YMID2 - Y (I1)))
      PRO2 = .5 *  ((D2B / D2) +  (D1A / D1))
      X2 = X (I8) -  (PRO2 *  (X (I8) - X (I3)))
      Y2 = Y (I8) -  (PRO2 *  (Y (I8) - Y (I3)))
      PRO3 = .5 *  ((D2A / D2) +  (D3B / D3))
      X3 = X (I2) -  (PRO3 *  (X (I2) - X (I7)))
      Y3 = Y (I2) -  (PRO3 *  (Y (I2) - Y (I7)))

C  AVERAGE POINTS TO GET THE CENTER

      XCEN2 =  (X1 + X2 + X3) / 3.
      YCEN2 =  (Y1 + Y2 + Y3) / 3.

      ERR = .FALSE.
      RETURN

      END
