C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE HIDMID (NLNKF, LINKF1, XN, YN, ZN, LINSET,
     &   IPSET, TVHMAX, ICROSS, XP, YP, ZP, NPART)
C=======================================================================

C   --*** HIDMID *** (MESH) Hide nodes beneath a face
C   --   Written by Amy Gilkey - revised 02/24/88
C   --
C   --HIDMID marks a node as hidden if it is beneath a face.
C   --
C   --Parameters:
C   --   NLNKF - IN - the number of nodes per face
C   --   LINKF1 - IN - the connectivity for the face
C   --   XN, YN, ZN - IN - the nodal coordinates
C   --   LINSET - IN - the sorted line set
C   --   IPSET - IN/OUT - the indices of the partial line set
C   --   NPART - IN/OUT - the number of lines in the partial line set
C   --
C   --Common Variables:

      include 'debug.blk'
      include 'dbnums.blk'
      include 'd3nums.blk'

      INTEGER LINKF1(NLNKF)
      REAL XN(*), YN(*), ZN(*)
      INTEGER LINSET(LLNSET,*)
      INTEGER IPSET(*)
      REAL TVHMAX(*)
      INTEGER ICROSS(*)
      REAL XP(*), YP(*), ZP(*)

      REAL X(4), Y(4), Z(4)
      REAL XNORM(4), YNORM(4), ZNORM(4)
      LOGICAL FIRSIN, FIRSTZ, FIRST(4)

C   --Calculate enclosing X-Y-Z box for face

      DO 100 ILINK = 1, NLNKF
         X(ILINK) = XN(LINKF1(ILINK))
         Y(ILINK) = YN(LINKF1(ILINK))
         Z(ILINK) = ZN(LINKF1(ILINK))
  100 CONTINUE

      XMIN = MIN (X(1), X(2), X(3), X(4))
      XMAX = MAX (X(1), X(2), X(3), X(4))
      YMIN = MIN (Y(1), Y(2), Y(3), Y(4))
      YMAX = MAX (Y(1), Y(2), Y(3), Y(4))
      ZMIN = MIN (Z(1), Z(2), Z(3), Z(4))
      ZMAX = MAX (Z(1), Z(2), Z(3), Z(4))

      FIRSIN = .TRUE.
      FIRSTZ = .TRUE.

C   --Check for nodes behind the face, mark as hidden

      IP = 1
  110 CONTINUE
      IF (IP .LE. NPART) THEN

C      --Check if outside of box enclosing face

         IF (XMAX .LT. XP(IP)) GOTO 120
         IF (XMIN .GT. XP(IP)) GOTO 120
         IF (YMAX .LT. YP(IP)) GOTO 120
         IF (YMIN .GT. YP(IP)) GOTO 120
         IF (ZMAX .LT. ZP(IP)) GOTO 120

         IV = LINSET(1,IPSET(IP))
         IH = LINSET(2,IPSET(IP))
         IF (((IV .EQ. LINKF1(1)) .OR. (IV .EQ. LINKF1(2)) .OR.
     &        (IV .EQ. LINKF1(3)) .OR. (IV .EQ. LINKF1(4))) .AND.
     &       ((IH .EQ. LINKF1(1)) .OR. (IH .EQ. LINKF1(2)) .OR.
     &        (IH .EQ. LINKF1(3)) .OR. (IH .EQ. LINKF1(4)))) GOTO 120

         IF (FIRSIN) THEN
            FIRSIN = .FALSE.

C         --Pre-compute factors for "within" test

            X2X1 = X(2) - X(1)
            Y2Y1 = Y(2) - Y(1)
            XY21 = X(1)*Y(2) - X(2)*Y(1)
            X3X2 = X(3) - X(2)
            Y3Y2 = Y(3) - Y(2)
            XY32 = X(2)*Y(3) - X(3)*Y(2)
            X4X3 = X(4) - X(3)
            Y4Y3 = Y(4) - Y(3)
            XY43 = X(3)*Y(4) - X(4)*Y(3)
            X1X4 = X(1) - X(4)
            Y1Y4 = Y(1) - Y(4)
            XY14 = X(4)*Y(1) - X(1)*Y(4)
         END IF

C      --Check if all 4 triangles formed by 2 nodes of the face
C      --and the node have a positive area

         X0 = XP(IP)
         Y0 = YP(IP)

         A12 = X2X1*Y0 - X0*Y2Y1 + XY21
         IF (A12 .LT. 0.0) GOTO 120
         A23 = X3X2*Y0 - X0*Y3Y2 + XY32
         IF (A23 .LT. 0.0) GOTO 120
         A34 = X4X3*Y0 - X0*Y4Y3 + XY43
         IF (A34 .LT. 0.0) GOTO 120
         A41 = X1X4*Y0 - X0*Y1Y4 + XY14
         IF (A41 .LT. 0.0) GOTO 120

C      --Check if the node's Z is behind the face's Z at that point

         IF (FIRSTZ) THEN
            FIRSTZ = .FALSE.
            FIRST(1) = .TRUE.
            FIRST(2) = .TRUE.
            FIRST(3) = .TRUE.
            FIRST(4) = .TRUE.

C         --Calculate epsilon based on length of "standard" side

            A = A12 + A23 + A34 + A41
            EPS = A * SQRT (A) * .001

C         --Calculate center of face

            XCEN = 0.25 * (X(1) + X(2) + X(3) + X(4))
            YCEN = 0.25 * (Y(1) + Y(2) + Y(3) + Y(4))
            ZCEN = 0.25 * (Z(1) + Z(2) + Z(3) + Z(4))
         END IF

C      --Check if the node's Z is behind the face's minimum Z

         IF ((ZMIN - ZP(IP)) .LT. EPS) THEN

C         --Calculate normal for center of face with two nodes
C         --that have the smallest area with the given node

            AMIN = MIN (A12, A23, A34, A41)

            IF (AMIN .EQ. A12) THEN
               ISIDE = 1
            ELSE IF (AMIN .EQ. A23) THEN
               ISIDE = 2
            ELSE IF (AMIN .EQ. A34) THEN
               ISIDE = 3
            ELSE IF (AMIN .EQ. A41) THEN
               ISIDE = 4
            END IF
            IF (FIRST(ISIDE)) THEN
               FIRST(ISIDE) = .FALSE.
               IF (ISIDE .LT. 4) THEN
                  N2 = ISIDE + 1
               ELSE
                  N2 = 1
               END IF
               AX = X(ISIDE) - XCEN
               AY = Y(ISIDE) - YCEN
               AZ = Z(ISIDE) - ZCEN
               BX = X(N2) - XCEN
               BY = Y(N2) - YCEN
               BZ = Z(N2) - ZCEN
               XNORM(ISIDE) = 0.5 * (AY*BZ - BY*AZ)
               YNORM(ISIDE) = 0.5 * (AZ*BX - BZ*AX)
               ZNORM(ISIDE) = 0.5 * (AX*BY - BX*AY)
            END IF

C         --Check if the node's Z is behind the face's Z at that point
            IF ((XNORM(ISIDE)*(X0-XCEN) + YNORM(ISIDE)*(Y0-YCEN)
     &         + ZNORM(ISIDE)*(ZP(IP)-ZCEN)) .GE. -EPS) GOTO 120
         END IF

C      --Hide the entire line (as is), and move the line to hidden section
         I = IPSET(IP)
         IPSET(IP) = IPSET(NPART)
         IPSET(NPART) = I
         T = TVHMAX(IP)
         TVHMAX(IP) = TVHMAX(NPART)
         TVHMAX(NPART) = T
         I = ICROSS(IP)
         ICROSS(IP) = ICROSS(NPART)
         ICROSS(NPART) = I
         XP(IP) = XP(NPART)
         YP(IP) = YP(NPART)
         ZP(IP) = ZP(NPART)
         NPART = NPART - 1
         IP = IP - 1

  120    CONTINUE
         IP = IP + 1
         GOTO 110
      END IF

      RETURN
      END
