C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE HIDNOD (NLNKF, LINKF1, XN, YN, ZN,
     &   COLMIN, ROWMIN, CRDELT, CREPS, NUMCOL, NUMROW,
     &   IXNCRS, IXNCRE, NPCR,
     &   HIDENP, nhid)
C=======================================================================

C   --*** HIDNOD *** (MESH) Hide nodes beneath a face
C   --   Written by Amy Gilkey - revised 10/22/87
C   --
C   --HIDNOD marks a node as hidden if it is beneath a face.
C   --
C   --Parameters:
C   --   NLNKF - IN - the number of nodes per face
C   --   LINKF1 - IN - the connectivity for the face
C   --   XN, YN, ZN - IN - the nodal coordinates
C   --   COLMIN, ROWMIN - IN - the minimum row and column value and 1 / interval
C   --   CRDELT - the row and column interval reciprical (1 / interval)
C   --   NUMCOL, NUMROW - IN - the number of rows and columns
C   --   IXNCRS, IXNCRE - OUT - the starting and ending NPCR index
C   --      for each column and row
C   --   NPCR - OUT - the nodes indexed by IXNCRS, IXNCRE
C   --   HIDENP - IN/OUT - node status (as in HIDDEN)
C   --
C   --Common Variables:

      PARAMETER (KNVIS=0, KNFOVR=10, KNHID=100)

      common /debugc/ cdebug
      common /debugn/ idebug
      character*8 cdebug

      include 'dbnums.blk'
      include 'd3nums.blk'

      INTEGER LINKF1(NLNKF)
      REAL XN(*), YN(*), ZN(*)
      INTEGER IXNCRS(0:NUMCOL,0:NUMROW), IXNCRE(0:NUMCOL,0:NUMROW)
      INTEGER NPCR(*)
      INTEGER HIDENP(*)

      REAL X(4), Y(4), Z(4)
      REAL XNORM(4), YNORM(4), ZNORM(4)
      LOGICAL FIRSTZ, FIRST(4)

C   --Skip face if all nodes hidden

      IF ((HIDENP(LINKF1(1)) .GT. KNVIS) .AND.
     &   (HIDENP(LINKF1(2)) .GT. KNVIS) .AND.
     &   (HIDENP(LINKF1(3)) .GT. KNVIS) .AND.
     &   (HIDENP(LINKF1(4)) .GT. KNVIS)) RETURN

C   --Calculate enclosing X-Y-Z box for face

      DO 100 ILINK = 1, 4
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

C   --Find the minimum and maximum row and column for this face

      MINCOL = INT ((XMIN - CREPS - COLMIN) * CRDELT)
      MAXCOL = INT ((XMAX + CREPS - COLMIN) * CRDELT)
      MINROW = INT ((YMIN - CREPS - ROWMIN) * CRDELT)
      MAXROW = INT ((YMAX + CREPS - ROWMIN) * CRDELT)

C   --Pre-compute factors for "within" test

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

      FIRSTZ = .TRUE.

C   --Check for nodes behind the face, mark as hidden

      DO 140 ICOL = MINCOL, MAXCOL
         DO 130 IROW = MINROW, MAXROW
            IXN = IXNCRS(ICOL,IROW)
  110       CONTINUE
            IF (IXN .LE. IXNCRE(ICOL,IROW)) THEN
               INP = NPCR(IXN)

C            --Check if outside of box enclosing face

               IF (XMAX .LT. XN(INP)) GOTO 120
               IF (XMIN .GT. XN(INP)) GOTO 120
               IF (YMAX .LT. YN(INP)) GOTO 120
               IF (YMIN .GT. YN(INP)) GOTO 120
               IF (ZMAX .LT. ZN(INP)) GOTO 120

C            --Check if node is part of the face

               IF ((INP .EQ. LINKF1(1)) .OR. (INP .EQ. LINKF1(2)) .OR.
     &            (INP .EQ. LINKF1(3)) .OR. (INP .EQ. LINKF1(4)))
     &            GOTO 120

C            --Check if all 4 triangles formed by 2 nodes of the face
C            --and the node have a positive area

               X0 = XN(INP)
               Y0 = YN(INP)
               A12 = X2X1*Y0 - X0*Y2Y1 + XY21
               IF (A12 .LT. 0.0) GOTO 120
               A23 = X3X2*Y0 - X0*Y3Y2 + XY32
               IF (A23 .LT. 0.0) GOTO 120
               A34 = X4X3*Y0 - X0*Y4Y3 + XY43
               IF (A34 .LT. 0.0) GOTO 120
               A41 = X1X4*Y0 - X0*Y1Y4 + XY14
               IF (A41 .LT. 0.0) GOTO 120

               IF (FIRSTZ) THEN
                  FIRSTZ = .FALSE.
                  FIRST(1) = .TRUE.
                  FIRST(2) = .TRUE.
                  FIRST(3) = .TRUE.
                  FIRST(4) = .TRUE.

C               --Calculate epsilon based on length of "standard" side

                  A = A12 + A23 + A34 + A41
                  EPS = A * SQRT (A) * .001

C               --Calculate center of face

                  XCEN = 0.25 * (X(1) + X(2) + X(3) + X(4))
                  YCEN = 0.25 * (Y(1) + Y(2) + Y(3) + Y(4))
                  ZCEN = 0.25 * (Z(1) + Z(2) + Z(3) + Z(4))
               END IF

C            --Check if the node's Z is behind the face's minimum Z

               IF ((ZMIN - ZN(INP)) .LT. EPS) THEN

C               --Calculate normal for center of face with two nodes
C               --that have the smallest area with the given node

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

C               --Check if the node's Z is behind the face's Z at that point

                  Z0C = XNORM(ISIDE)*(X0-XCEN) + YNORM(ISIDE)*(Y0-YCEN)
     &               + ZNORM(ISIDE)*(ZN(INP)-ZCEN)

                  IF (Z0C .GT. EPS) GOTO 120
                  IF (Z0C .GE. -EPS) THEN
C                  --Mark as slide line
                     HIDENP(INP) = -1
                     GOTO 120
                  END IF
               END IF

C            --Mark the node as hidden and delete from visible set
               HIDENP(INP) = KNFOVR
               I = IXNCRE(ICOL,IROW)
               NPCR(IXN) = NPCR(I)
               IXNCRE(ICOL,IROW) = I-1
               IXN = IXN - 1
               nhid = nhid + 1

  120          CONTINUE
               IXN = IXN + 1
               GOTO 110
            END IF
  130    CONTINUE
  140 CONTINUE

      RETURN
      END
