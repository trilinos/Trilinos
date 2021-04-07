C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE PLTIQC(ICOLOR,R,G,B)
      DIMENSION CARRAY(3)

      CALL VDIQCO(1,ICOLOR,CARRAY,0)
      R = CARRAY(1)
      G = CARRAY(2)
      B = CARRAY(3)
      RETURN

      END
