C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE PLTIQD(ARRAY)
      DIMENSION ARRAY(*)

      DO 2040 I = 1,23
         CALL VDIQDC(I,ARRAY(I))
 2040 CONTINUE
      RETURN

      END
