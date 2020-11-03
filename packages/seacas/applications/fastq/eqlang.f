C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE EQLANG (MXND, XN, YN, LXN, NODE, N0, N2, NFROM, DIST,
     &   VRO, XDEL, YDEL)
C***********************************************************************

C  SUBROUTINE EQLANG = CALCULATES A VECTOR SUM THAT ATTEMPTS TO
C                      MAINTAIN EQUAL ANGLES FOR A NODE

C***********************************************************************

      DIMENSION XN(MXND), YN(MXND), LXN(4, MXND)

      LOGICAL EXPAND

      PI = ATAN2(0.0, -1.0)
      TWOPI = 2.0 * PI

      IF (NFROM .GT. 0) THEN

C  TEST FOR THE EXPANSION CASE

         IF ( ( ((LXN (4, NFROM) .NE. 0) .AND.
     &      (LXN (2, NFROM) .LT. 0)) .OR.
     &      ((LXN (4, NFROM) .LT. 0) .AND.
     &      (LXN (2, NFROM) .GT. 0)) )
     &      .AND.
     &      ((LXN (3, N0) .EQ. 0) .OR. (LXN (3, N2) .EQ. 0)) ) THEN
            EXPAND = .TRUE.
         ELSE
            EXPAND = .FALSE.
         ENDIF

         ANG1 = ATAN2 ( YN (N2) - YN (NFROM), XN (N2) - XN (NFROM))
         IF (ANG1 .LT. 0.) ANG1 = ANG1 + TWOPI
         ANG2 = ATAN2 ( YN (N0) - YN (NFROM), XN (N0) - XN (NFROM))
         IF (ANG2 .LT. 0.) ANG2 = ANG2 + TWOPI
         ANG3 = ATAN2 ( YN (NODE) - YN (NFROM), XN (NODE) - XN (NFROM))
         IF (ANG3 .LT. 0.) ANG3 = ANG3 + TWOPI

C  GET THE APPROPRIATE ANGLE BETWEEN ANGLE 1 AND 2

         ANG12D = ANG2 - ANG1
         IF (ANG12D .LT. 0.) ANG12D = ANG12D + TWOPI

C  IF THIS IS AN EXPANSION, THEN ADJUST THE ANGLE ACCORDINGLY

         IF (EXPAND) THEN
            IF (LXN (3, N2) .EQ. 0) THEN
               ANG12 = ANG1 + (ANG12D * .6)
            ELSEIF (LXN (3, N0) .EQ. 0) THEN
               ANG12 = ANG1 + (ANG12D * .4)
            ELSE
               ANG12 = ANG1 + (ANG12D * .5)
            ENDIF
         ELSE
            ANG12 = ANG1 + (ANG12D * .5)
         ENDIF
         IF (ANG12 .GT. TWOPI) ANG12 = ANG12 - TWOPI

C  GET THE AVERAGE ANGLE BETWEEN ANGLE 12 AND 3

         IF (ANG12 .GT. ANG3) THEN
            ANG3D = ANG12 - ANG3
            IF (ANG3D .GT. PI) THEN
               ANG = ANG12 + ((TWOPI - ANG3D) * .5)
            ELSE
               ANG = ANG12 - (ANG3D * .5)
            ENDIF
         ELSE
            ANG3D = ANG3 - ANG12
            IF (ANG3D .GT. PI) THEN
               ANG = ANG3 + ((TWOPI - ANG3D) * .5)
            ELSE
               ANG = ANG3 - (ANG3D * .5)
            ENDIF
         ENDIF

C  GET THE DISTANCE TO MAKE THE OUTSIDE FLAT AT THIS ANGLE

         D1 = SQRT ( ((XN (NFROM) - XN (N0)) ** 2) +
     &      ((YN (NFROM) - YN (N0)) ** 2) )
         D2 = SQRT ( ((XN (N2) - XN (N0)) ** 2) +
     &      ((YN (N2) - YN (N0)) ** 2) )
         D3 = SQRT ( ((XN (NFROM) - XN (N2)) ** 2) +
     &      ((YN (NFROM) - YN (N2)) ** 2) )
         ARG = (SIN (ANG12D) * D1) / D2
         IF (ARG .GT. 1.0) ARG = 1.0
         IF (ARG .LT. -1.0) ARG = -1.0
         BETA = ASIN (ARG)
         D0 = (D3 * SIN (BETA)) / SIN (PI - BETA - (ANG12D * .5))

         IF (D0 .GT. DIST) THEN
            IF (EXPAND) THEN
               DIST0 = D0
            ELSE
               DIST0 = (DIST + D0) * .5
            ENDIF
         ELSE
            DIST0 = DIST
         ENDIF

C  CALCULATE THE NEW COORDINATES

         X0 = XN (NFROM) + (COS (ANG) * DIST0)
         Y0 = YN (NFROM) + (SIN (ANG) * DIST0)
         XDEL = (X0 - XN (NODE)) * VRO
         YDEL = (Y0 - YN (NODE)) * VRO

      ELSE
         XDEL = 0.
         YDEL = 0.
      ENDIF

      RETURN

      END
