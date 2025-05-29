C    Copyright(C) 1999-2020, 2025 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE TRBLK (IDELB, NUMELB, NUMLNK, MAT, NELBLK, NNODES)
      DIMENSION IDELB(*), NUMELB(*), NUMLNK(*), MAT(7,*)

      IBEG = 0
      IEND = 0
      DO I=1, NELBLK
         if (nnodes .eq. 0) then
            nnodes = numlnk(i)
         else if (numlnk(i) .ne. nnodes) then
          CALL PRTERR('ERROR',
     *      'All blocks must be the same tet/hex/quad topology.')
          stop
       endif
       IBEG = IEND + 1
       IEND = IEND + NUMELB(I)
       MAT(1,I) = IDELB(I)
       MAT(2,I) = NUMELB(I)
       MAT(3,I) = IBEG
       MAT(4,I) = IEND
       MAT(5,I) = 1
       MAT(6,I) = I
       MAT(7,I) = numlnk(i)
      end do
      RETURN
      END
