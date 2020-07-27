C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE MXVECT(N,VEC,MAT,RES)
      REAL VEC(*),MAT(N,*),RES(*)

      IF (N .EQ. 4) THEN
        RES(1) = MAT(1,1)*VEC(1) + MAT(2,1)*VEC(2) + MAT(3,1)*VEC(3) +
     *    MAT(4,1)*VEC(4)

        RES(2) = MAT(1,2)*VEC(1) + MAT(2,2)*VEC(2) + MAT(3,2)*VEC(3) +
     *    MAT(4,2)*VEC(4)

        RES(3) = MAT(1,3)*VEC(1) + MAT(2,3)*VEC(2) + MAT(3,3)*VEC(3) +
     *    MAT(4,3)*VEC(4)

        RES(4) = MAT(1,4)*VEC(1) + MAT(2,4)*VEC(2) + MAT(3,4)*VEC(3) +
     *    MAT(4,4)*VEC(4)

      ELSE
        DO 2980 J = 1,N
          RES(J) = 0.0
 2980   CONTINUE
        DO 3010 I = 1,N
          DO 2990 J = 1,N
            RES(J) = RES(J) + MAT(I,J)*VEC(I)
 2990     CONTINUE
 3010   CONTINUE
      END IF
      RETURN

      END
