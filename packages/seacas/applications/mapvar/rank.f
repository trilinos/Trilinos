C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

      SUBROUTINE RANK(N,INDX,IRANK,NDIM)

C...  create a rank array from an index array
      DIMENSION INDX(NDIM),IRANK(NDIM)

      DO J=1,N
        IRANK(INDX(J))=J
      end do

      RETURN
      END

