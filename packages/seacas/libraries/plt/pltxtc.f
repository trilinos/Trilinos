C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE PLTXTC1(X,Y,LINE)
      CHARACTER*(*) LINE

      CALL PLTRIM(LINE,L)
      IF (L.LE.0) THEN
         RETURN

      END IF

      CALL PLTXSL(LINE(1:L),XL)
      CALL PLTXTS(X-XL/2.,Y,LINE(1:L))
      RETURN

      END
