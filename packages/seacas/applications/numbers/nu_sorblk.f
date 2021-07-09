C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE SORBLK (IDELB, INDEX, MAT, NELBLK)
      DIMENSION IDELB(*), INDEX(*), MAT(6,*)

      CALL INDEXI (IDELB, INDEX, NELBLK, .TRUE.)
      DO 10 I=1, NELBLK
         MAT(6,I) = INDEX(I)
   10 CONTINUE
      RETURN
      END
