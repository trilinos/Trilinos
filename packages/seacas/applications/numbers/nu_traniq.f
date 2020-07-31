C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE TRANIQ (LSTSN, MAP, MASSLV, NSEG, IDIM)
      DIMENSION LSTSN(*), MAP(*), MASSLV(IDIM,*)

      DO 10 I=1,NSEG
          MASSLV(1,MAP(I)) = LSTSN(I)
   10 CONTINUE

      RETURN
      END
