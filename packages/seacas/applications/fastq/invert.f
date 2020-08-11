C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE INVERT_FQ (MXND, MLN, XN, YN, ZN, LXK, KXL, NXL, LXN,
     *  LLL, LNODES, XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX, DEV1, KREG,
     *  NODE, XDEL, YDEL)
C***********************************************************************

C  SUBROUTINE INVERT = CHECKS FOR AN INVERSION OR CROSSING OF A BOUNDARY
C                      UPON ITSELF AND CORRECTS IT WHERE NECESSARY

C***********************************************************************

      DIMENSION XN(MXND), YN(MXND), ZN(MXND)
      DIMENSION LXN(4, MXND), NXL(2, 3*MXND)
      DIMENSION LXK(4, MXND), KXL(2, 3*MXND)
      DIMENSION LNODES (MLN, MXND)

      LOGICAL ERR, VCROSS

      CHARACTER*3 DEV1

      ERR = .FALSE.

      XOLD = XN (NODE)
      YOLD = YN (NODE)

      N2 = NODE
      N3 = LNODES (3, N2)
      N4 = LNODES (3, N3)
      N1 = LNODES (2, N2)
      N0 = LNODES (2, N1)

C  GET THE ANGLES BEFORE MOVEMENT

      IF (LXN (4, N1) .EQ. 0)
     &   CALL GETANG (MXND, MLN, XN, YN, LNODES, LXK, KXL, NXL,
     &   LXN, N0, N1, N2, ANG1A, ERR)
      IF (LXN (4, N2) .EQ. 0)
     &   CALL GETANG (MXND, MLN, XN, YN, LNODES, LXK, KXL, NXL,
     &   LXN, N1, N2, N3, ANG2A, ERR)
      IF (LXN (4, N3) .EQ. 0)
     &   CALL GETANG (MXND, MLN, XN, YN, LNODES, LXK, KXL, NXL,
     &   LXN, N2, N3, N4, ANG3A, ERR)

C  NOW PLACE THE NODE TEMPORARILY AT THE NEW PROPOSED LOCATION

      XN (NODE) = XN (NODE) + XDEL
      YN (NODE) = YN (NODE) + YDEL

C  GET THE ANGLE BEING ADJUSTED AT THE NODE ITSELF

      IF ((LXN (4, N2) .EQ. 0) .AND. (ANG2A .GT. 0.)) THEN
         CALL GETANG (MXND, MLN, XN, YN, LNODES, LXK, KXL, NXL,
     &      LXN, N1, N2, N3, ANG2B, ERR)

C  ADJUST THE NODE LOCATION IF NECESSARY

         IF (ANG2B .LT. 0.) THEN
            CALL VINTER (MXND, XN, YN, N1, N3, N2, XOLD, YOLD,
     &         XNEW, YNEW, VCROSS)
            IF (VCROSS) THEN
               XN (NODE) = XNEW
               YN (NODE) = YNEW
            ENDIF
         ENDIF
      ENDIF

C  GET THE ANGLE BEING ADJUSTED ON THE CCW SIDE OF THIS NODE

      IF ((LXN (4, N1) .EQ. 0) .AND. (ANG1A .GT. 0.)) THEN
         CALL GETANG (MXND, MLN, XN, YN, LNODES, LXK, KXL, NXL,
     &      LXN, N0, N1, N2, ANG1B, ERR)

C  ADJUST THE NODE LOCATION IF NECESSARY

         IF (ANG1B .LT. 0.) THEN
            CALL VINTER (MXND, XN, YN, N1, N0, N2, XOLD, YOLD,
     &         XNEW, YNEW, VCROSS)
            IF (VCROSS) THEN
               XN (NODE) = XNEW
               YN (NODE) = YNEW
            ENDIF
         ENDIF
      ENDIF

C  GET THE ANGLE BEING ADJUSTED ON THE CW SIDE OF THIS NODE

      IF ((LXN (4, N3) .EQ. 0) .AND. (ANG3A .GT. 0.)) THEN
         CALL GETANG (MXND, MLN, XN, YN, LNODES, LXK, KXL, NXL,
     &      LXN, N2, N3, N4, ANG3B, ERR)

C  ADJUST THE NODE LOCATION IF NECESSARY

         IF (ANG3B .LT. 0.) THEN
            CALL VINTER (MXND, XN, YN, N3, N4, N2, XOLD, YOLD,
     &         XNEW, YNEW, VCROSS)
            IF (VCROSS) THEN
               XN (NODE) = XNEW
               YN (NODE) = YNEW
            ENDIF
         ENDIF
      ENDIF

C  RESTORE THE OLD LOCATION AND THE XDEL AND YDEL TO THE CORRECTED
C  VALUES

      XDEL = XN (NODE) - XOLD
      YDEL = YN (NODE) - YOLD
      XN (NODE) = XOLD
      YN (NODE) = YOLD

      RETURN

      END
