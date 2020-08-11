C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE MXMULT(N,MAT1,MAT2,MATR)
      REAL MAT1(N,*),MAT2(N,*),MATR(N,*)

      IF (N .EQ. 4) THEN
          MATR(1,1) = MAT1(1,1)*MAT2(1,1) + MAT1(1,2)*MAT2(2,1) +
     *                MAT1(1,3)*MAT2(3,1) + MAT1(1,4)*MAT2(4,1)
          MATR(1,2) = MAT1(1,1)*MAT2(1,2) + MAT1(1,2)*MAT2(2,2) +
     *                MAT1(1,3)*MAT2(3,2) + MAT1(1,4)*MAT2(4,2)
          MATR(1,3) = MAT1(1,1)*MAT2(1,3) + MAT1(1,2)*MAT2(2,3) +
     *                MAT1(1,3)*MAT2(3,3) + MAT1(1,4)*MAT2(4,3)
          MATR(1,4) = MAT1(1,1)*MAT2(1,4) + MAT1(1,2)*MAT2(2,4) +
     *                MAT1(1,3)*MAT2(3,4) + MAT1(1,4)*MAT2(4,4)

          MATR(2,1) = MAT1(2,1)*MAT2(1,1) + MAT1(2,2)*MAT2(2,1) +
     *                MAT1(2,3)*MAT2(3,1) + MAT1(2,4)*MAT2(4,1)
          MATR(2,2) = MAT1(2,1)*MAT2(1,2) + MAT1(2,2)*MAT2(2,2) +
     *                MAT1(2,3)*MAT2(3,2) + MAT1(2,4)*MAT2(4,2)
          MATR(2,3) = MAT1(2,1)*MAT2(1,3) + MAT1(2,2)*MAT2(2,3) +
     *                MAT1(2,3)*MAT2(3,3) + MAT1(2,4)*MAT2(4,3)
          MATR(2,4) = MAT1(2,1)*MAT2(1,4) + MAT1(2,2)*MAT2(2,4) +
     *                MAT1(2,3)*MAT2(3,4) + MAT1(2,4)*MAT2(4,4)

          MATR(3,1) = MAT1(3,1)*MAT2(1,1) + MAT1(3,2)*MAT2(2,1) +
     *                MAT1(3,3)*MAT2(3,1) + MAT1(3,4)*MAT2(4,1)
          MATR(3,2) = MAT1(3,1)*MAT2(1,2) + MAT1(3,2)*MAT2(2,2) +
     *                MAT1(3,3)*MAT2(3,2) + MAT1(3,4)*MAT2(4,2)
          MATR(3,3) = MAT1(3,1)*MAT2(1,3) + MAT1(3,2)*MAT2(2,3) +
     *                MAT1(3,3)*MAT2(3,3) + MAT1(3,4)*MAT2(4,3)
          MATR(3,4) = MAT1(3,1)*MAT2(1,4) + MAT1(3,2)*MAT2(2,4) +
     *                MAT1(3,3)*MAT2(3,4) + MAT1(3,4)*MAT2(4,4)

          MATR(4,1) = MAT1(4,1)*MAT2(1,1) + MAT1(4,2)*MAT2(2,1) +
     *                MAT1(4,3)*MAT2(3,1) + MAT1(4,4)*MAT2(4,1)
          MATR(4,2) = MAT1(4,1)*MAT2(1,2) + MAT1(4,2)*MAT2(2,2) +
     *                MAT1(4,3)*MAT2(3,2) + MAT1(4,4)*MAT2(4,2)
          MATR(4,3) = MAT1(4,1)*MAT2(1,3) + MAT1(4,2)*MAT2(2,3) +
     *                MAT1(4,3)*MAT2(3,3) + MAT1(4,4)*MAT2(4,3)
          MATR(4,4) = MAT1(4,1)*MAT2(1,4) + MAT1(4,2)*MAT2(2,4) +
     *                MAT1(4,3)*MAT2(3,4) + MAT1(4,4)*MAT2(4,4)

      ELSE

        DO 230 K = 1,N
          DO 200 J = 1,N
            MATR(K,J) = 0.0
 200      CONTINUE
          DO 220 I = 1,N
            DO 210 J = 1,N
              MATR(K,J) = MATR(K,J) + MAT1(K,I)*MAT2(I,J)
 210        CONTINUE
 220      CONTINUE
 230    CONTINUE
      END IF
      RETURN
      END
