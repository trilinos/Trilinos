C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      REAL FUNCTION SPIRAL  (XA,  XK,  X,  XCEN,  YCEN,  ANGLE)
C***********************************************************************

C  FUNCTION SPIRAL = CALCULATES THE Y VALUUE GIVEN THE SPIRAL AND X

C***********************************************************************

      SPIRAL  =  XA * EXP (XK * ANGLE) * COS (ANGLE) - (X - XCEN)

      RETURN

      END
