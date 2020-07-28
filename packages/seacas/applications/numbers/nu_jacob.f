C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      subroutine jacob(x1,x2,x3,x4,x5,x6,x7,x8,
     *                 y1,y2,y3,y4,y5,y6,y7,y8,
     *                 z1,z2,z3,z4,z5,z6,z7,z8, JACMN)
      REAL JACOBI, JACMN
      REAL XXI(3), XET(3), XZE(3)

      JACMN = 1.0e38

C... Jacobian at mid-point:

      xxi(1) = x2 + x3 + x6 + x7 - x1 - x4 - x5 - x8
      xet(1) = x3 + x4 + x7 + x8 - x1 - x2 - x5 - x6
      xze(1) = x5 + x6 + x7 + x8 - x1 - x2 - x3 - x4

      xxi(2) = y2 + y3 + y6 + y7 - y1 - y4 - y5 - y8
      xet(2) = y3 + y4 + y7 + y8 - y1 - y2 - y5 - y6
      xze(2) = y5 + y6 + y7 + y8 - y1 - y2 - y3 - y4

      xxi(3) = z2 + z3 + z6 + z7 - z1 - z4 - z5 - z8
      xet(3) = z3 + z4 + z7 + z8 - z1 - z2 - z5 - z6
      xze(3) = z5 + z6 + z7 + z8 - z1 - z2 - z3 - z4

      jacobi = CRSDOT(xxi,xet,xze)/64.0
      jacmn = min(jacobi, jacmn)
C... J(0,0,0):

      xxi(1) = x2 - x1
      xet(1) = x4 - x1
      xze(1) = x5 - x1

      xxi(2) = y2 - y1
      xet(2) = y4 - y1
      xze(2) = y5 - y1

      xxi(3) = z2 - z1
      xet(3) = z4 - z1
      xze(3) = z5 - z1

      jacobi = CRSDOT(xxi,xet,xze)
      jacmn = min(jacobi, jacmn)

C... J(1,0,0):

      xxi(1) = x2 - x1
      xet(1) = x3 - x2
      xze(1) = x6 - x2

      xxi(2) = y2 - y1
      xet(2) = y3 - y2
      xze(2) = y6 - y2

      xxi(3) = z2 - z1
      xet(3) = z3 - z2
      xze(3) = z6 - z2

      jacobi = CRSDOT(xxi,xet,xze)
      jacmn = min(jacobi, jacmn)

C... J(0,1,0):

      xxi(1) = x3 - x4
      xet(1) = x4 - x1
      xze(1) = x8 - x4

      xxi(2) = y3 - y4
      xet(2) = y4 - y1
      xze(2) = y8 - y4

      xxi(3) = z3 - z4
      xet(3) = z4 - z1
      xze(3) = z8 - z4

      jacobi = CRSDOT(xxi,xet,xze)
      jacmn = min(jacobi, jacmn)

C... J(0,0,1):

      xxi(1) = x6 - x5
      xet(1) = x8 - x5
      xze(1) = x5 - x1

      xxi(2) = y6 - y5
      xet(2) = y8 - y5
      xze(2) = y5 - y1

      xxi(3) = z6 - z5
      xet(3) = z8 - z5
      xze(3) = z5 - z1

      jacobi = CRSDOT(xxi,xet,xze)
      jacmn = min(jacobi, jacmn)

C... J(1,1,0):

      xxi(1) = x3 - x4
      xet(1) = x3 - x2
      xze(1) = x7 - x3

      xxi(2) = y3 - y4
      xet(2) = y3 - y2
      xze(2) = y7 - y3

      xxi(3) = z3 - z4
      xet(3) = z3 - z2
      xze(3) = z7 - z3

      jacobi = CRSDOT(xxi,xet,xze)
      jacmn = min(jacobi, jacmn)

C... J(1,0,1):

      xxi(1) = x6 - x5
      xet(1) = x7 - x6
      xze(1) = x6 - x2

      xxi(2) = y6 - y5
      xet(2) = y7 - y6
      xze(2) = y6 - y2

      xxi(3) = z6 - z5
      xet(3) = z7 - z6
      xze(3) = z6 - z2

      jacobi = CRSDOT(xxi,xet,xze)
      jacmn = min(jacobi, jacmn)

C... J(0,1,1):

      xxi(1) = x7 - x8
      xet(1) = x8 - x5
      xze(1) = x8 - x4

      xxi(2) = y7 - y8
      xet(2) = y8 - y5
      xze(2) = y8 - y4

      xxi(3) = z7 - z8
      xet(3) = z8 - z5
      xze(3) = z8 - z4

      jacobi = CRSDOT(xxi,xet,xze)
      jacmn = min(jacobi, jacmn)

C... J(1,1,1):

      xxi(1) = x7 - x8
      xet(1) = x7 - x6
      xze(1) = x7 - x3

      xxi(2) = y7 - y8
      xet(2) = y7 - y6
      xze(2) = y7 - y3

      xxi(3) = z7 - z8
      xet(3) = z7 - z6
      xze(3) = z7 - z3

      jacobi = CRSDOT(xxi,xet,xze)
      jacmn = min(jacobi, jacmn)

      return
      end

      REAL FUNCTION CRSDOT(a, b, c)
      REAL a(3), b(3), c(3)
C ... Compute d = a dot (b cross c)

      REAL xcross, ycross, zcross

      xcross = b(2) * c(3) - b(3) * c(2)
      ycross = b(3) * c(1) - b(1) * c(3)
      zcross = b(1) * c(2) - b(2) * c(1)

      crsdot = (a(1) * xcross + a(2) * ycross + a(3) * zcross)
      return
      end
