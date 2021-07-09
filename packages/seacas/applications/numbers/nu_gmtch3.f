C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE GMTCH3 (COORD, DIRCOS, MASSLV, NIQSLV, TDIS,
     *    ITMP, NIQM, NIQS, DMAX, NUMNP)
      DIMENSION COORD (NUMNP,*), DIRCOS(5,*), MASSLV(2,*), NIQSLV(*),
     *    TDIS(3,*), ITMP(*)

      DO 30 IMAS = 1, NIQM
          X1 = COORD (MASSLV(1, IMAS),1)
          Y1 = COORD (MASSLV(1, IMAS),2)
          Z1 = COORD (MASSLV(1, IMAS),3)

C DIRCOS is average unit direction vector from surfaces at node

          DCS1 = DIRCOS (1, IMAS)
          DCS2 = DIRCOS (2, IMAS)
          DCS3 = DIRCOS (3, IMAS)

          DO 10 ISLV = 1, NIQS
              X0 = COORD (NIQSLV(ISLV),1)
              Y0 = COORD (NIQSLV(ISLV),2)
              Z0 = COORD (NIQSLV(ISLV),3)

              TDIS(3,ISLV) = (X0 - X1)**2 + (Y0 - Y1)**2 + (Z0 - Z1)**2
              T = -( DCS1 * (X1-X0) + DCS2 * (Y1-Y0) + DCS3 * (Z1-Z0) )
              TDIS(1,ISLV) = T
              TDIS(2,ISLV) = ABS( TDIS(3,ISLV) - T**2 )
   10     CONTINUE

          TMIN = DMAX
          NMIN = 0
          DO 20 ISLV = 1, NIQS
             IF (TDIS(2,ISLV) .LT. TMIN .AND. TDIS(3,ISLV) .LE. DMAX)
     *          THEN
                  TMIN = TDIS(2,ISLV)
                  NMIN = ISLV
              END IF
   20     CONTINUE

          IF (NMIN .NE. 0) THEN
              DIRCOS(4,IMAS) = TDIS(1,NMIN)
              DIRCOS(5,IMAS) = SQRT(TDIS(2,NMIN))
              MASSLV(2,IMAS) = NIQSLV(NMIN)
          ELSE
              MASSLV(2,IMAS) = 0
          END IF
   30 CONTINUE
      CALL CULL3 (DIRCOS, MASSLV, NIQM)
      RETURN
      END
