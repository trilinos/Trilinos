C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE MXIDEN(N,MAT)
      REAL MAT(N,*)

      IF (N .EQ. 4) THEN
        MAT(1,1) = 1.0
        MAT(1,2) = 0.0
        MAT(1,3) = 0.0
        MAT(1,4) = 0.0

        MAT(2,1) = 0.0
        MAT(2,2) = 1.0
        MAT(2,3) = 0.0
        MAT(2,4) = 0.0

        MAT(3,1) = 0.0
        MAT(3,2) = 0.0
        MAT(3,3) = 1.0
        MAT(3,4) = 0.0

        MAT(4,1) = 0.0
        MAT(4,2) = 0.0
        MAT(4,3) = 0.0
        MAT(4,4) = 1.0
      ELSE
        DO 2930 I = 1,N
          DO 2950 J = 1,N
            MAT(I,J) = 0.
 2950     CONTINUE
 2930   CONTINUE
        DO 2970 I = 1,N
          MAT(I,I) = 1.
 2970   CONTINUE
      end if
      RETURN

      END
