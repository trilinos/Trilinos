C    Copyright(C) 1988-2017 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    Redistribution and use in source and binary forms, with or without
C    modification, are permitted provided that the following conditions are
C    met:
C
C    * Redistributions of source code must retain the above copyright
C       notice, this list of conditions and the following disclaimer.
C
C    * Redistributions in binary form must reproduce the above
C      copyright notice, this list of conditions and the following
C      disclaimer in the documentation and/or other materials provided
C      with the distribution.
C
C    * Neither the name of NTESS nor the names of its
C      contributors may be used to endorse or promote products derived
C      from this software without specific prior written permission.
C
C    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
C    "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
C    LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
C    A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
C    OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
C    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
C    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
C    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
C    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
C    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
C    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

C $Id: con2d.f,v 1.5 2004/06/29 18:05:32 gdsjaar Exp $
C $Log: con2d.f,v $
C Revision 1.5  2004/06/29 18:05:32  gdsjaar
C General cleanup. Remove unused labels and variables.
C
C Revision 1.4  2000/07/06 18:07:42  gdsjaar
C Fix assumption that variables are saved between subroutine calls
C
C Revision 1.3  2000/07/06 16:49:57  gdsjaar
C Changed real*4 to real
C
C Revision 1.2  1991/04/10 19:28:18  gdsjaar
C Removed vax debug lines
C
c Revision 1.1.1.1  1991/02/21  15:42:40  gdsjaar
c NUMBERS: Greg Sjaardema, initial Unix release
c
c Revision 1.1  1991/02/21  15:42:39  gdsjaar
c Initial revision
c
C=======================================================================
      SUBROUTINE CON2D (CRD, NDIM, NUMNP, IX, NNODES, NUMEL, MAT,
     *   NELBLK, SELECT, ASPECT, SKEW, TAPER, AREA,
     *   SUMRY, ISUMRY, DEBUG)
C=======================================================================
C
C *** CON2D *** Calculate state of mesh -- Aspect ratio, Skewness,
C               and Taper
C
C    (Greg Sjaardema, 16 April, 1989)
C
C Based on article by John Robinson, "CRE Method of element testing
C    and the Jacobian shape parameters," Eng. Comput., 1987, Vol. 4,
C    June, pp 113 - 118
C
C -- ARRAYS:
C     CRD(NUMNP, NDIM)  - IN -
C     IX(NNODES, NUMEL) - IN -
C     MAT(5, NELBLK)    - IN -
C     SELECT(NUMEL)     - IN -
C     ASPECT(NUMEL)     - OUT- Aspect ratio (1.0 <= AR <= infinity)
C     SKEW(NUMEL)       - OUT- Skewness of mesh, degrees (0 <= skew <= ?)
C     TAPER(NUMEL)      - OUT- Taper of mesh, combination of X and Y taper
C     AREA(NUMEL)       - OUT- Area of element
C
C -- SCALARS:
C     NDIM   - Number of spatial dimensions
C     NUMNP  - Number of nodal points
C     NNODES - Number of nodes per element
C     NUMEL  - Number of elements
C     NELBLK - Number of material/element blocks
C
C         E2                                               E4
C     +----------+          +-----------+       +---------+ |
C     |          | F3      /           /       /           \
C     |          |        / A         /       /             \
C     +----------+       +-----------+       +---------------+
C
C      AR = E2/F3        SKEW = SIN(A)            TAPER Y
C
C=======================================================================
C
      REAL   CRD(NUMNP, NDIM), ASPECT(*), SKEW(*), TAPER(*), AREA(*)
      INTEGER  IX(NNODES, NUMEL), MAT(6, NELBLK), ISUMRY(2,4,NELBLK)
      REAL   SUMRY(4,4,NELBLK)
      LOGICAL  SELECT(*), DEBUG, ISABRT
      include 'nu_io.blk'

      CALL INIREA (NUMEL, 0.0, ASPECT)
      CALL INIREA (NUMEL, 0.0, SKEW)
      CALL INIREA (NUMEL, 0.0, TAPER)
      CALL INIREA (NUMEL, 0.0, AREA)

      DO 40 IBLK = 1, NELBLK
         IF (MAT(5, IBLK) .EQ. 1) THEN
            IELBEG = MAT(3, IBLK)
            IELEND = MAT(4, IBLK)
            IF (ISABRT()) RETURN
            DO 30 IEL = IELBEG, IELEND
c .. if doesn't vectorize, remove next line and do select on summary only
               IF (SELECT(IEL)) THEN
                  X1 = CRD(IX(1,IEL), 1)
                  X2 = CRD(IX(2,IEL), 1)
                  X3 = CRD(IX(3,IEL), 1)
                  X4 = CRD(IX(4,IEL), 1)

                  Y1 = CRD(IX(1,IEL), 2)
                  Y2 = CRD(IX(2,IEL), 2)
                  Y3 = CRD(IX(3,IEL), 2)
                  Y4 = CRD(IX(4,IEL), 2)

C ... Make centroid of element the center of coordinate system

                  XS = (X1 + X2 + X3 + X4) / 4.0
                  YS = (Y1 + Y2 + Y3 + Y4) / 4.0

                  X1 = X1 - XS
                  X2 = X2 - XS
                  X3 = X3 - XS
                  X4 = X4 - XS

                  Y1 = Y1 - YS
                  Y2 = Y2 - YS
                  Y3 = Y3 - YS
                  Y4 = Y4 - YS

C ... Rotate element such that center of side 2-3 and 4-1 define X axis

                  AMAG = SQRT((X2+X3-X4-X1)**2 + (Y2+Y3-Y4-Y1)**2)
                  C = (X2 + X3 - X4 - X1) / AMAG
                  S = (Y2 + Y3 - Y4 - Y1) / AMAG

                  XT =  C * X1 + S * Y1
                  Y1 = -S * X1 + C * Y1
                  X1 = XT

                  XT =  C * X2 + S * Y2
                  Y2 = -S * X2 + C * Y2
                  X2 = XT

                  XT =  C * X3 + S * Y3
                  Y3 = -S * X3 + C * Y3
                  X3 = XT

                  XT =  C * X4 + S * Y4
                  Y4 = -S * X4 + C * Y4
                  X4 = XT

C ... Calculate ``Shape function'' parameters - E1, F1, F2 = 0.0

                  E2 = -X1 + X2 + X3 - X4
                  E3 = -X1 - X2 + X3 + X4
                  E4 =  X1 - X2 + X3 - X4

                  F3 = -Y1 - Y2 + Y3 + Y4
                  F4 =  Y1 - Y2 + Y3 - Y4

                  ASPECT(IEL) = MAX (E2 / F3, F3 / E2)
                  SKEW(IEL)   = ABS(E3/F3) / SQRT((E3/F3)**2 + 1.0)
                  TAPER(IEL)  = SQRT((F4 / F3)**2 + (E4 / E2)**2)
                  AREA(IEL)   = E2 * F3 / 4.0

               END IF
   30       CONTINUE
         END IF
   40 CONTINUE
      DO 50 IO=IOMIN, IOMAX
         WRITE (IO, 60)
   50 CONTINUE
   60 FORMAT (/'  Shape Parameters for Selected Elements',/
     *   /2X,'Mat    Minimum   Elem    Maximum   Elem    Average  ',
     *   '  Std. Dev.'/)

      DO 100 ITMP = 1, NELBLK
         IBLK = MAT(6, ITMP)
         IF (MAT(5, IBLK) .EQ. 1) THEN
            IELBEG = MAT(3, IBLK)
            IELEND = MAT(4, IBLK)
            NUMLST = IELEND - IELBEG + 1
C ... Determine mins, maxs, averages, and std. dev. for selected block/elem

            CALL SUMMRY (' ', NUMLST, SELECT(IELBEG), ASPECT(IELBEG),
     *         SUMRY(1,1,IBLK), ISUMRY(1,1,IBLK),IELBEG-1)
            CALL SUMMRY (' ', NUMLST, SELECT(IELBEG), SKEW(IELBEG),
     *         SUMRY(1,2,IBLK), ISUMRY(1,2,IBLK),IELBEG-1)
            CALL SUMMRY (' ', NUMLST, SELECT(IELBEG), TAPER(IELBEG),
     *         SUMRY(1,3,IBLK), ISUMRY(1,3,IBLK),IELBEG-1)
            CALL SUMMRY (' ', NUMLST, SELECT(IELBEG), AREA(IELBEG),
     *         SUMRY(1,4,IBLK), ISUMRY(1,4,IBLK),IELBEG-1)

            DO 70 IO= IOMIN, IOMAX
               WRITE (IO, 80) MAT(1,IBLK),
     *            SUMRY(1,1,IBLK), ISUMRY(1,1,IBLK),
     *            SUMRY(2,1,IBLK), ISUMRY(2,1,IBLK),
     *            SUMRY(3,1,IBLK),  SUMRY(4,1,IBLK), 'Aspect',
     *            SUMRY(1,2,IBLK), ISUMRY(1,2,IBLK),
     *            SUMRY(2,2,IBLK), ISUMRY(2,2,IBLK),
     *            SUMRY(3,2,IBLK),  SUMRY(4,2,IBLK), 'Skewness',
     *            SUMRY(1,3,IBLK), ISUMRY(1,3,IBLK),
     *            SUMRY(2,3,IBLK), ISUMRY(2,3,IBLK),
     *            SUMRY(3,3,IBLK),  SUMRY(4,3,IBLK), 'Taper',
     *            SUMRY(1,4,IBLK), ISUMRY(1,4,IBLK),
     *            SUMRY(2,4,IBLK), ISUMRY(2,4,IBLK),
     *            SUMRY(3,4,IBLK),  SUMRY(4,4,IBLK), 'Area'

               WRITE (IO, 90)
   70       CONTINUE
   80       FORMAT (I5,(T8,2(1PE15.8,I6,2X),2(1PE15.8,2X),3X,A8))
   90       FORMAT (' ------------------')
         END IF
  100 CONTINUE
      END
