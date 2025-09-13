C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE ELCENT (ELCEN, IX, COORD, NDIM, NUMEL, NELNOD, NUMNP)
      DIMENSION ELCEN(NUMEL, *), IX(NELNOD, *), COORD(NUMNP, *)

      IF (NDIM .EQ. 2) THEN
          DO I=1, NUMEL
              ELCEN(I,1) = (COORD(IX(1,I),1) + COORD(IX(2,I),1) +
     *                      COORD(IX(3,I),1) + COORD(IX(4,I),1))/4.0
              ELCEN(I,2) = (COORD(IX(1,I),2) + COORD(IX(2,I),2) +
     *                      COORD(IX(3,I),2) + COORD(IX(4,I),2))/4.0
           end do
      ELSE
         IF (NELNOD .EQ. 8) then
          DO I=1, NUMEL
              ELCEN(I,1) = (COORD(IX(1,I),1) + COORD(IX(2,I),1) +
     *                      COORD(IX(3,I),1) + COORD(IX(4,I),1) +
     *                      COORD(IX(5,I),1) + COORD(IX(6,I),1) +
     *                      COORD(IX(7,I),1) + COORD(IX(8,I),1))/8.0

              ELCEN(I,2) = (COORD(IX(1,I),2) + COORD(IX(2,I),2) +
     *                      COORD(IX(3,I),2) + COORD(IX(4,I),2) +
     *                      COORD(IX(5,I),2) + COORD(IX(6,I),2) +
     *                      COORD(IX(7,I),2) + COORD(IX(8,I),2))/8.0

              ELCEN(I,3) = (COORD(IX(1,I),3) + COORD(IX(2,I),3) +
     *                      COORD(IX(3,I),3) + COORD(IX(4,I),3) +
     *                      COORD(IX(5,I),3) + COORD(IX(6,I),3) +
     *                      COORD(IX(7,I),3) + COORD(IX(8,I),3))/8.0
           end do
        end if
         IF (NELNOD .EQ. 4) then
          DO I=1, NUMEL
              ELCEN(I,1) = (COORD(IX(1,I),1) + COORD(IX(2,I),1) +
     *                      COORD(IX(3,I),1) + COORD(IX(4,I),1))/4.0

              ELCEN(I,2) = (COORD(IX(1,I),2) + COORD(IX(2,I),2) +
     *                      COORD(IX(3,I),2) + COORD(IX(4,I),2))/4.0

              ELCEN(I,3) = (COORD(IX(1,I),3) + COORD(IX(2,I),3) +
     *                      COORD(IX(3,I),3) + COORD(IX(4,I),3))/4.0
           end do
        end if

      END IF
      RETURN
      END
