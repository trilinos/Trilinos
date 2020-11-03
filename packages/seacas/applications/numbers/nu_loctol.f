C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE LOCTOL (TYPE, NDIM, RV, KV)

C     This routine is used to set the tolerances and distances
C        used in the LOCATE routines.
C     If a tolerance is not entered (blank field), then
C        the tolerance is set to the entered distance value, and
C        the distance is set to 0.0
C     If a tolerance is entered, the values are returned with no
C        changes

      DIMENSION RV(*), KV(*)
      CHARACTER*(*) TYPE
      LOGICAL MATSTR

         IF (MATSTR(TYPE, 'LINE', 1)) THEN
            IF (NDIM .EQ. 3) THEN
               IF (KV(11) .EQ. -1) THEN
                  RV(11) = RV(10)
                  RV(10) = 0.0
               END IF
            ELSE
               IF (KV(9) .EQ. -1) THEN
                  RV(9) = RV(8)
                  RV(8) = 0.0
               END IF
            END IF
         ELSE IF (MATSTR(TYPE, 'PLANE', 2)) THEN
            IF (NDIM .EQ. 3) THEN
               IF (KV(11) .EQ. -1) THEN
                  RV(11) = RV(10)
                  RV(10) = 0.0
               END IF
            ELSE
               CONTINUE
            END IF
         ELSE IF (MATSTR(TYPE, 'POINT', 2)) THEN
            IF (NDIM .EQ. 3) THEN
               IF (KV(8) .EQ. -1) THEN
                  RV(8) = RV(7)
                  RV(7) = 0.0
               END IF
            ELSE
               IF (KV(7) .EQ. -1) THEN
                  RV(7) = RV(6)
                  RV(6) = 0.0
               END IF
            END IF
         END IF
      RETURN
      END
