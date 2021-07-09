C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE CON3D (CRD, NDIM, NUMNP, IX, NNODES, NUMEL, MAT,
     *   NELBLK, SELECT, ASPECT1, ASPECT2, ASPECT3, AREA,
     *   SUMRY, ISUMRY, SKEWX, SKEWY, SKEWZ, TAPERX, TAPERY, TAPERZ,
     *   JAC, DEBUG)
C=======================================================================

C *** CON3D *** Calculate state of mesh -- Aspect ratio, Skewness,
C               and Taper  (Three-dimensional mesh)

C    (Greg Sjaardema, 16 April, 1989)

C Based on article by John Robinson, "CRE Method of element testing
C    and the Jacobian shape parameters," Eng. Comput., 1987, Vol. 4,
C    June, pp 113 - 118

C     NOTE:  This is an experimental routine.

C -- ARRAYS:
C     CRD(NUMNP, NDIM)  - IN -
C     IX(NNODES, NUMEL) - IN -
C     MAT(5, NELBLK)    - IN -
C     SELECT(NUMEL)     - IN -
C     ASPECT(NUMEL)     - OUT- Aspect ratio (1.0 <= AR <= infinity)
C     SKEW(NUMEL)       - OUT- Skewness of mesh, degrees (0 <= skew <= ?)
C     TAPER(NUMEL)      - OUT- Taper of mesh, combination of X and Y taper
C     AREA(NUMEL)       - OUT- Area of element

C -- SCALARS:
C     NDIM   - Number of spatial dimensions
C     NUMNP  - Number of nodal points
C     NNODES - Number of nodes per element
C     NUMEL  - Number of elements
C     NELBLK - Number of material/element blocks

C         E2                                               E4
C     +----------+          +-----------+       +---------+ |
C     |          | F3      /           /       /           \
C     |          |        / A         /       /             \
C     +----------+       +-----------+       +---------------+

C      AR = E2/F3        SKEW = SIN(A)            TAPER Y

C=======================================================================

C      REAL   CRD(NUMNP, NDIM), ASPECT(*), SKEW(*), TAPER(*), AREA(*)
      REAL   CRD(NUMNP, NDIM)
      REAL   ASPECT1(*), ASPECT2(*), ASPECT3(*), AREA(*)
      INTEGER  IX(NNODES, NUMEL), MAT(6, NELBLK), ISUMRY(2,4,NELBLK)
      REAL   SUMRY(4,4,NELBLK)
      LOGICAL  SELECT(*), DEBUG, ISABRT
      REAL   SKEWX(*), SKEWY(*), SKEWZ(*)
      REAL   TAPERX(*), TAPERY(*), TAPERZ(*)
      REAL   JAC(*)

C      REAL   CRD(8,3)
C      REAL   ASPECT1(1), ASPECT2(1), ASPECT3(1), AREA(1)
C      INTEGER  IX(8,1), MAT(5, 1), ISUMRY(2,4,1)
C      REAL   SUMRY(4,4,1)
C      LOGICAL  SELECT(1)
C      REAL   SKEWX(1), SKEWY(1)
C      REAL   TAPERX(1), TAPERY(1), TAPERZ(1)

C      DATA IX/1,2,3,4,5,6,7,8/
C      DATA MAT /1,1,1,1,1/
C      DATA SELECT /.TRUE./
C      DATA NELBLK /1/

      include 'nu_io.blk'

C      IOMIN = 6
C      IOMAX = 6

      DO 50 IBLK = 1, NELBLK
         IF (MAT(5, IBLK) .EQ. 1) THEN
            IELBEG = MAT(3, IBLK)
            IELEND = MAT(4, IBLK)
            DO 40 IEL = IELBEG, IELEND
               IF (ISABRT()) RETURN
c .. if doesn't vectorize, remove next line and do select on summary only
               IF (SELECT(IEL)) THEN
                  X1 = CRD(IX(1,IEL), 1)
                  X2 = CRD(IX(2,IEL), 1)
                  X3 = CRD(IX(3,IEL), 1)
                  X4 = CRD(IX(4,IEL), 1)
                  X5 = CRD(IX(5,IEL), 1)
                  X6 = CRD(IX(6,IEL), 1)
                  X7 = CRD(IX(7,IEL), 1)
                  X8 = CRD(IX(8,IEL), 1)

                  Y1 = CRD(IX(1,IEL), 2)
                  Y2 = CRD(IX(2,IEL), 2)
                  Y3 = CRD(IX(3,IEL), 2)
                  Y4 = CRD(IX(4,IEL), 2)
                  Y5 = CRD(IX(5,IEL), 2)
                  Y6 = CRD(IX(6,IEL), 2)
                  Y7 = CRD(IX(7,IEL), 2)
                  Y8 = CRD(IX(8,IEL), 2)

                  Z1 = CRD(IX(1,IEL), 3)
                  Z2 = CRD(IX(2,IEL), 3)
                  Z3 = CRD(IX(3,IEL), 3)
                  Z4 = CRD(IX(4,IEL), 3)
                  Z5 = CRD(IX(5,IEL), 3)
                  Z6 = CRD(IX(6,IEL), 3)
                  Z7 = CRD(IX(7,IEL), 3)
                  Z8 = CRD(IX(8,IEL), 3)

                  E2 = -X1 + X2 + X3 - X4 - X5 + X6 + X7 - X8
                  F3 = -Y1 - Y2 + Y3 + Y4 - Y5 - Y6 + Y7 + Y8
                  G4 = -Z1 - Z2 - Z3 - Z4 + Z5 + Z6 + Z7 + Z8

C ... Make centroid of element the center of coordinate system

                  XS = (X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8) / 8.0
                  YS = (Y1 + Y2 + Y3 + Y4 + Y5 + Y6 + Y7 + Y8) / 8.0
                  ZS = (Z1 + Z2 + Z3 + Z4 + Z5 + Z6 + Z7 + Z8) / 8.0

                  X1 = X1 - XS
                  X2 = X2 - XS
                  X3 = X3 - XS
                  X4 = X4 - XS
                  X5 = X5 - XS
                  X6 = X6 - XS
                  X7 = X7 - XS
                  X8 = X8 - XS

                  Y1 = Y1 - YS
                  Y2 = Y2 - YS
                  Y3 = Y3 - YS
                  Y4 = Y4 - YS
                  Y5 = Y5 - YS
                  Y6 = Y6 - YS
                  Y7 = Y7 - YS
                  Y8 = Y8 - YS

                  Z1 = Z1 - ZS
                  Z2 = Z2 - ZS
                  Z3 = Z3 - ZS
                  Z4 = Z4 - ZS
                  Z5 = Z5 - ZS
                  Z6 = Z6 - ZS
                  Z7 = Z7 - ZS
                  Z8 = Z8 - ZS

C ... Rotate element such that center of side 2-3 and 4-1 define X axis

                  X2376 = X2 + X3 + X6 + X7
                  X1485 = X1 + X4 + X8 + X5

                  Y2376 = Y2 + Y3 + Y6 + Y7
                  Y1485 = Y1 + Y4 + Y8 + Y5

                  Z2376 = Z2 + Z3 + Z7 + Z6
                  Z1485 = Z1 + Z4 + Z8 + Z5

                  DX = X2376 - X1485
                  DY = Y2376 - Y1485
                  DZ = Z2376 - Z1485

                  AMAGX = SQRT(DX**2 + DZ**2)
                  AMAGY = SQRT(DX**2 + DY**2 + DZ**2)
                  FMAGX = SIGN(0.5, AMAGX) + SIGN(0.5, -AMAGX)
                  FMAGY = SIGN(0.5, AMAGY) + SIGN(0.5, -AMAGY)

                  CZ = DX / (AMAGX + FMAGX) + FMAGX
                  SZ = DZ / (AMAGX + FMAGX)
                  CY = SQRT(DX**2 + DZ**2) / (AMAGY + FMAGY) + FMAGY
                  SY = DY / (AMAGY + FMAGY)

                  XT =  CY * CZ * X1 + CY * SZ * Z1 + SY * Y1
                  Y1 = -SY * CZ * X1 - SY * SZ * Z1 + CY * Y1
                  Z1 = -SZ * X1      +      CZ * Z1
                  X1 = XT

                  XT =  CY * CZ * X2 + CY * SZ * Z2 + SY * Y2
                  Y2 = -SY * CZ * X2 - SY * SZ * Z2 + CY * Y2
                  Z2 = -SZ * X2      +      CZ * Z2
                  X2 = XT

                  XT =  CY * CZ * X3 + CY * SZ * Z3 + SY * Y3
                  Y3 = -SY * CZ * X3 - SY * SZ * Z3 + CY * Y3
                  Z3 = -SZ * X3      +      CZ * Z3
                  X3 = XT

                  XT =  CY * CZ * X4 + CY * SZ * Z4 + SY * Y4
                  Y4 = -SY * CZ * X4 - SY * SZ * Z4 + CY * Y4
                  Z4 = -SZ * X4      +      CZ * Z4
                  X4 = XT

                  XT =  CY * CZ * X5 + CY * SZ * Z5 + SY * Y5
                  Y5 = -SY * CZ * X5 - SY * SZ * Z5 + CY * Y5
                  Z5 = -SZ * X5      +      CZ * Z5
                  X5 = XT

                  XT =  CY * CZ * X6 + CY * SZ * Z6 + SY * Y6
                  Y6 = -SY * CZ * X6 - SY * SZ * Z6 + CY * Y6
                  Z6 = -SZ * X6      +      CZ * Z6
                  X6 = XT

                  XT =  CY * CZ * X7 + CY * SZ * Z7 + SY * Y7
                  Y7 = -SY * CZ * X7 - SY * SZ * Z7 + CY * Y7
                  Z7 = -SZ * X7      +      CZ * Z7
                  X7 = XT

                  XT =  CY * CZ * X8 + CY * SZ * Z8 + SY * Y8
                  Y8 = -SY * CZ * X8 - SY * SZ * Z8 + CY * Y8
                  Z8 = -SZ * X8      +      CZ * Z8
                  X8 = XT

C ... Now, rotate about Y

                  X3487 = X3 + X4 + X8 + X7
                  Y3487 = Y3 + Y4 + Y8 + Y7
                  Z3487 = Z3 + Z4 + Z8 + Z7

                  X1562 = X1 + X5 + X6 + X2
                  Y1562 = Y1 + Y5 + Y6 + Y2
                  Z1562 = Z1 + Z5 + Z6 + Z2

                  DX = X3487 - X1562
                  DY = Y3487 - Y1562
                  DZ = Z3487 - Z1562

                  AMAGY = SQRT(DY**2 + DZ**2)
                  FMAGY = SIGN(0.5, AMAGY) + SIGN(0.5, -AMAGY)

                  CX = DY / (AMAGY + FMAGY) + FMAGY
                  SX = DZ / (AMAGY + FMAGY)

                  YT =  CX * Y1 + SX * Z1
                  Z1 = -SX * Y1 + CX * Z1
                  Y1 = YT

                  YT =  CX * Y2 + SX * Z2
                  Z2 = -SX * Y2 + CX * Z2
                  Y2 = YT

                  YT =  CX * Y3 + SX * Z3
                  Z3 = -SX * Y3 + CX * Z3
                  Y3 = YT

                  YT =  CX * Y4 + SX * Z4
                  Z4 = -SX * Y4 + CX * Z4
                  Y4 = YT

                  YT =  CX * Y5 + SX * Z5
                  Z5 = -SX * Y5 + CX * Z5
                  Y5 = YT

                  YT =  CX * Y6 + SX * Z6
                  Z6 = -SX * Y6 + CX * Z6
                  Y6 = YT

                  YT =  CX * Y7 + SX * Z7
                  Z7 = -SX * Y7 + CX * Z7
                  Y7 = YT

                  YT =  CX * Y8 + SX * Z8
                  Z8 = -SX * Y8 + CX * Z8
                  Y8 = YT

C ... Calculate ``Shape function'' parameters - E1, F1, F2 = 0.0

                  E1 =  X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8
                  E2 = -X1 + X2 + X3 - X4 - X5 + X6 + X7 - X8
                  E3 = -X1 - X2 + X3 + X4 - X5 - X6 + X7 + X8
                  E4 = -X1 - X2 - X3 - X4 + X5 + X6 + X7 + X8
                  E5 =  X1 + X2 - X3 - X4 - X5 - X6 + X7 + X8
                  E6 =  X1 - X2 - X3 + X4 - X5 + X6 + X7 - X8
                  E7 =  X1 - X2 + X3 - X4 + X5 - X6 + X7 - X8
                  E8 = -X1 + X2 - X3 + X4 + X5 - X6 + X7 - X8

                  F1 =  Y1 + Y2 + Y3 + Y4 + Y5 + Y6 + Y7 + Y8
                  F2 = -Y1 + Y2 + Y3 - Y4 - Y5 + Y6 + Y7 - Y8
                  F3 = -Y1 - Y2 + Y3 + Y4 - Y5 - Y6 + Y7 + Y8
                  F4 = -Y1 - Y2 - Y3 - Y4 + Y5 + Y6 + Y7 + Y8
                  F5 =  Y1 + Y2 - Y3 - Y4 - Y5 - Y6 + Y7 + Y8
                  F6 =  Y1 - Y2 - Y3 + Y4 - Y5 + Y6 + Y7 - Y8
                  F7 =  Y1 - Y2 + Y3 - Y4 + Y5 - Y6 + Y7 - Y8
                  F8 = -Y1 + Y2 - Y3 + Y4 + Y5 - Y6 + Y7 - Y8

                  G1 =  Z1 + Z2 + Z3 + Z4 + Z5 + Z6 + Z7 + Z8
                  G2 = -Z1 + Z2 + Z3 - Z4 - Z5 + Z6 + Z7 - Z8
                  G3 = -Z1 - Z2 + Z3 + Z4 - Z5 - Z6 + Z7 + Z8
                  G4 = -Z1 - Z2 - Z3 - Z4 + Z5 + Z6 + Z7 + Z8
                  G5 =  Z1 + Z2 - Z3 - Z4 - Z5 - Z6 + Z7 + Z8
                  G6 =  Z1 - Z2 - Z3 + Z4 - Z5 + Z6 + Z7 - Z8
                  G7 =  Z1 - Z2 + Z3 - Z4 + Z5 - Z6 + Z7 - Z8
                  G8 = -Z1 + Z2 - Z3 + Z4 + Z5 - Z6 + Z7 - Z8

                  IF (DEBUG) THEN
                  WRITE (99, 20) 'Element Number ',IEL
                  WRITE (99, 10) 'COLM', 1., 2., 3., 4., 5., 6., 7., 8.
                  WRITE (99, 10) 'CSCS', CX, SX, CY, SY, CZ, SZ
                  WRITE (99, 10) 'X', X1,X2,X3,X4,X5,X6,X7,X8
                  WRITE (99, 10) 'Y', Y1,Y2,Y3,Y4,Y5,Y6,Y7,Y8
                  WRITE (99, 10) 'Z', Z1,Z2,Z3,Z4,Z5,Z6,Z7,Z8
                  WRITE (99, 10) 'E', E1,E2,E3,E4,E5,E6,E7,E8
                  WRITE (99, 10) 'F', F1,F2,F3,F4,F5,F6,F7,F8
                  WRITE (99, 10) 'G', G1,G2,G3,G4,G5,G6,G7,G8
                  END IF

   10             FORMAT (1X,A4,8(2X,F7.3))
   20             FORMAT (/,1X,A,I6)

                  IF (DEBUG) THEN
                  IF (E2 .EQ. 0.0) THEN
                     WRITE (*, 30) 'E2', IEL
                     E2 = 1.0
                  END IF
                  IF (F3 .EQ. 0.0) THEN
                     WRITE (*, 30) 'F3', IEL
                     F3 = 1.0
                  END IF
                  IF (G4 .EQ. 0.0) THEN
                     WRITE (*, 30) 'G4', IEL
                     G4 = 1.0
                  END IF
                  END IF
   30             FORMAT (1X,A2,' zero at element ',I6)

                  ASPECT1(IEL) = MAX (E2, F3) / MIN(E2, F3)
                  ASPECT2(IEL) = MAX (E2, G4) / MIN(E2, G4)
                  ASPECT3(IEL) = MAX (F3, G4) / MIN(F3, G4)
                  AREA(IEL)   = E2 * F3 * G4 / 64.0

                  SKEWX(IEL) = ABS(F4/G4) / SQRT((F4/G4)**2 + 1.0)
                  SKEWY(IEL) = ABS(E4/G4) / SQRT((E4/G4)**2 + 1.0)
                  SKEWZ(IEL) = ABS(E3/F3) / SQRT((E3/F3)**2 + 1.0)

                  TAPERX(IEL) = E5 / E2
                  TAPERY(IEL) = F5 / F3
                  TAPERZ(IEL) = G5 / G4

                  IF (DEBUG) THEN
                  WRITE (99, 10) 'ASPS', ASPECT1(IEL), ASPECT2(IEL),
     &               ASPECT3(IEL), AREA(IEL)
                  WRITE (99, 10) 'SKEW', SKEWX(IEL),   SKEWY(IEL),
     &               SKEWZ(IEL)
                  WRITE (99, 10) 'TAPR', TAPERX(IEL),  TAPERY(IEL),
     &               TAPERZ(IEL)
                  END IF

                  CALL JACOB(x1,x2,x3,x4,x5,x6,x7,x8,
     *                       y1,y2,y3,y4,y5,y6,y7,y8,
     *                       z1,z2,z3,z4,z5,z6,z7,z8, JAC(IEL))

               END IF
   40       CONTINUE
         END IF
   50 CONTINUE

C      IF (.FALSE.) THEN
      IF (.TRUE.) THEN

         DO 60 IO=IOMIN, IOMAX
            WRITE (IO, 70)
            WRITE (IO, 80)
   60    CONTINUE

   70    FORMAT (/'  **** EXPERIMENTAL ROUTINE - CHECK RESULTS ****')
   80    FORMAT (/'  Shape Parameters for Selected Elements',/
     *      /1X,' Mat   Minimum    Elem     Maximum    Elem',
     *     '     Average    Std. Dev.'/)

         DO 140 ITMP = 1, NELBLK
            IBLK = MAT(6, ITMP)
            IF (MAT(5, IBLK) .EQ. 1) THEN
               IELBEG = MAT(3, IBLK)
               IELEND = MAT(4, IBLK)
               NUMLST = IELEND - IELBEG + 1
C ... Determine mins, maxs, averages, and std. dev. for selected block/elem

               CALL SUMMRY (' ',NUMLST, SELECT(IELBEG), ASPECT1(IELBEG),
     *            SUMRY(1,1,IBLK), ISUMRY(1,1,IBLK),IELBEG-1)
               CALL SUMMRY (' ',NUMLST, SELECT(IELBEG), ASPECT2(IELBEG),
     *            SUMRY(1,2,IBLK), ISUMRY(1,2,IBLK),IELBEG-1)
               CALL SUMMRY (' ',NUMLST, SELECT(IELBEG), ASPECT3(IELBEG),
     *            SUMRY(1,3,IBLK), ISUMRY(1,3,IBLK),IELBEG-1)
               CALL SUMMRY (' ', NUMLST, SELECT(IELBEG), AREA(IELBEG),
     *            SUMRY(1,4,IBLK), ISUMRY(1,4,IBLK),IELBEG-1)

               DO 90 IO= IOMIN, IOMAX
                  WRITE (IO, 120) MAT(1,IBLK),
     *               SUMRY(1,1,IBLK), ISUMRY(1,1,IBLK),
     *               SUMRY(2,1,IBLK), ISUMRY(2,1,IBLK),
     *               SUMRY(3,1,IBLK),  SUMRY(4,1,IBLK), 'AspectXY',
     *               SUMRY(1,2,IBLK), ISUMRY(1,2,IBLK),
     *               SUMRY(2,2,IBLK), ISUMRY(2,2,IBLK),
     *               SUMRY(3,2,IBLK),  SUMRY(4,2,IBLK), 'AspectXZ',
     *               SUMRY(1,3,IBLK), ISUMRY(1,3,IBLK),
     *               SUMRY(2,3,IBLK), ISUMRY(2,3,IBLK),
     *               SUMRY(3,3,IBLK),  SUMRY(4,3,IBLK), 'AspectYZ',
     *               SUMRY(1,4,IBLK), ISUMRY(1,4,IBLK),
     *               SUMRY(2,4,IBLK), ISUMRY(2,4,IBLK),
     *               SUMRY(3,4,IBLK),  SUMRY(4,4,IBLK), 'Volume'
   90          CONTINUE

               CALL SUMMRY (' ', NUMLST, SELECT(IELBEG), SKEWX(IELBEG),
     *            SUMRY(1,1,IBLK), ISUMRY(1,1,IBLK),IELBEG-1)
               CALL SUMMRY (' ', NUMLST, SELECT(IELBEG), SKEWY(IELBEG),
     *            SUMRY(1,2,IBLK), ISUMRY(1,2,IBLK),IELBEG-1)
               CALL SUMMRY (' ', NUMLST, SELECT(IELBEG), SKEWZ(IELBEG),
     *            SUMRY(1,3,IBLK), ISUMRY(1,3,IBLK),IELBEG-1)
               DO 100 IO= IOMIN, IOMAX
                  WRITE (IO, 120) MAT(1,IBLK),
     *               SUMRY(1,1,IBLK), ISUMRY(1,1,IBLK),
     *               SUMRY(2,1,IBLK), ISUMRY(2,1,IBLK),
     *               SUMRY(3,1,IBLK),  SUMRY(4,1,IBLK), 'Skew X',
     *               SUMRY(1,2,IBLK), ISUMRY(1,2,IBLK),
     *               SUMRY(2,2,IBLK), ISUMRY(2,2,IBLK),
     *               SUMRY(3,2,IBLK),  SUMRY(4,2,IBLK), 'Skew Y',
     *               SUMRY(1,3,IBLK), ISUMRY(1,3,IBLK),
     *               SUMRY(2,3,IBLK), ISUMRY(2,3,IBLK),
     *               SUMRY(3,3,IBLK),  SUMRY(4,3,IBLK), 'Skew Z'
  100          CONTINUE

               CALL SUMMRY (' ', NUMLST, SELECT(IELBEG), TAPERX(IELBEG),
     *            SUMRY(1,1,IBLK), ISUMRY(1,1,IBLK),IELBEG-1)
               CALL SUMMRY (' ', NUMLST, SELECT(IELBEG), TAPERY(IELBEG),
     *            SUMRY(1,2,IBLK), ISUMRY(1,2,IBLK),IELBEG-1)
               CALL SUMMRY (' ', NUMLST, SELECT(IELBEG), TAPERZ(IELBEG),
     *            SUMRY(1,3,IBLK), ISUMRY(1,3,IBLK),IELBEG-1)
               CALL SUMMRY (' ', NUMLST, SELECT(IELBEG), JAC(IELBEG),
     *            SUMRY(1,4,IBLK), ISUMRY(1,4,IBLK),IELBEG-1)
               DO 110 IO= IOMIN, IOMAX
                  WRITE (IO, 120) MAT(1,IBLK),
     *               SUMRY(1,1,IBLK), ISUMRY(1,1,IBLK),
     *               SUMRY(2,1,IBLK), ISUMRY(2,1,IBLK),
     *               SUMRY(3,1,IBLK),  SUMRY(4,1,IBLK), 'TaperX',
     *               SUMRY(1,2,IBLK), ISUMRY(1,2,IBLK),
     *               SUMRY(2,2,IBLK), ISUMRY(2,2,IBLK),
     *               SUMRY(3,2,IBLK),  SUMRY(4,2,IBLK), 'TaperY',
     *               SUMRY(1,3,IBLK), ISUMRY(1,3,IBLK),
     *               SUMRY(2,3,IBLK), ISUMRY(2,3,IBLK),
     *               SUMRY(3,3,IBLK),  SUMRY(4,3,IBLK), 'TaperZ',
     *               SUMRY(1,4,IBLK), ISUMRY(1,4,IBLK),
     *               SUMRY(2,4,IBLK), ISUMRY(2,4,IBLK),
     *               SUMRY(3,4,IBLK),  SUMRY(4,4,IBLK), 'MinJacob'
                  WRITE (IO, 130)
  110          CONTINUE

  120          FORMAT (I4,(T7,2(1PE15.8,I9,1X),2(1PE15.8,2X),A8))
  130          FORMAT (' ------------------')
            END IF
  140    CONTINUE
      END IF
      END

