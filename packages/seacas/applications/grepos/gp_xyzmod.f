C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE XYZMOD (NUMNP, NDIM, XOLD,YOLD,ZOLD,  XNEW,YNEW,ZNEW)
C=======================================================================
      REAL XOLD(*), YOLD(*), ZOLD(*)
      REAL XNEW(*), YNEW(*), ZNEW(*)

      DO 10 I=1, NUMNP
        XNEW(I) = XOLD(I)
        YNEW(I) = YOLD(I)
        IF (NDIM .EQ. 3) ZNEW(I) = ZOLD(I)
 10   CONTINUE
      RETURN
      END
