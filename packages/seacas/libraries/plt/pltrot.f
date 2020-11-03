C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE PLTROT(R,A,B)
      REAL R(3,3),A(3,3),B(3,3)

      B(1,1) = R(1,1)*A(1,1) + R(2,1)*A(2,1) + R(3,1)*A(3,1)
      B(2,1) = R(1,2)*A(1,1) + R(2,2)*A(2,1) + R(3,2)*A(3,1)
      B(3,1) = R(1,3)*A(1,1) + R(2,3)*A(2,1) + R(3,3)*A(3,1)

      B(1,2) = R(1,1)*A(1,2) + R(2,1)*A(2,2) + R(3,1)*A(3,2)
      B(2,2) = R(1,2)*A(1,2) + R(2,2)*A(2,2) + R(3,2)*A(3,2)
      B(3,2) = R(1,3)*A(1,2) + R(2,3)*A(2,2) + R(3,3)*A(3,2)

      B(1,3) = R(1,1)*A(1,3) + R(2,1)*A(2,3) + R(3,1)*A(3,3)
      B(2,3) = R(1,2)*A(1,3) + R(2,2)*A(2,3) + R(3,2)*A(3,3)
      B(3,3) = R(1,3)*A(1,3) + R(2,3)*A(2,3) + R(3,3)*A(3,3)

      RETURN

      END
