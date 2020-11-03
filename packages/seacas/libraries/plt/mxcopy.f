C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE MXCOPY(N,MAT1,MAT2)
      REAL MAT1(N,*),MAT2(N,*)

      IF (N .EQ. 4) THEN
        DO 100 I=1, 4
          MAT2(I,1) = MAT1(I,1)
          MAT2(I,2) = MAT1(I,2)
          MAT2(I,3) = MAT1(I,3)
          MAT2(I,4) = MAT1(I,4)
 100    CONTINUE
      ELSE
        DO 2910 I = 1,N
          DO 2890 J = 1,N
            MAT2(I,J) = MAT1(I,J)
 2890     CONTINUE
 2910   CONTINUE
      end if
      RETURN

      END
