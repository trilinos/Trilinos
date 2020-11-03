C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE ZXYZ (X, Y, Z, MAP, NUMNP, NDIM)
C=======================================================================

C -- X,Y,Z -- REAL - IN/OUT - Coordinates of nodes
C -- MAP   -- INT  - IN     - Map between new and old node numbers
C                             MAP(I) for inactive nodes = NUMNP+1
C    NOTE: MAP(I) <= I for all I unless node is deleted, then
C          MAP(I)  = NUMNP+1

      REAL    X(NUMNP+1), Y(NUMNP+1), Z(NUMNP+1)
      INTEGER MAP(NUMNP)

      DO 10 I=1, NUMNP
         X(MAP(I)) = X(I)
         Y(MAP(I)) = Y(I)
         IF (NDIM .EQ. 3) THEN
            Z(MAP(I)) = Z(I)
         END IF
   10 CONTINUE

      RETURN
      END

