C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE CULL3 (DIRCOS, MASSLV, NIQM)
      DIMENSION DIRCOS(5,*), MASSLV(2,*)

      J = 0
      DO 10 I=1, NIQM
          IF (MASSLV(2,I) .NE. 0) THEN
              J = J + 1
              MASSLV(1, J) = MASSLV(1, I)
              MASSLV(2, J) = MASSLV(2, I)
              DIRCOS(1, J) = DIRCOS(1, I)
              DIRCOS(2, J) = DIRCOS(2, I)
              DIRCOS(3, J) = DIRCOS(3, I)
              DIRCOS(4, J) = DIRCOS(4, I)
              DIRCOS(5, J) = DIRCOS(5, I)
          END IF
   10 CONTINUE

      NIQM = J
      RETURN
      END
