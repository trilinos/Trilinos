C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE SKPRNT (LUNIT, LEN, RSTACK, NDIM)
C***********************************************************************

C  SUBROUTINE SKPRNT = STACK PRINTING ROUTINE

C***********************************************************************

      REAL RSTACK (NDIM)

      WRITE (LUNIT, '(I8,G12.5)') (I, RSTACK(I + 2), I = LEN, 1, -1)

      RETURN

      END
