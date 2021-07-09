C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE CFMAX (VARNP, NLNKF, LINKF1, FEMIN, FEMAX)
C=======================================================================

C   --*** CFMAX *** (DETOUR) Compute min/max values for face
C   --   Written by Amy Gilkey - revised 03/14/88
C   --
C   --CFMAX computes the minimum and maximum value of the nodal variable
C   --for the face.
C   --
C   --Parameters:
C   --   VARNP - IN - the contour function values
C   --   NLNKF - IN - the number of nodes per face
C   --   LINKF1 - IN - the connectivity for the face
C   --   FEMIN, FEMAX - OUT - minimum and maximum value

      REAL VARNP(*)
      INTEGER LINKF1(NLNKF)

C   --Compute the minimum and maximum values for the face

      FEMAX = VARNP(LINKF1(1))
      FEMIN = VARNP(LINKF1(1))
      DO 100 K = 2, NLNKF
         FEMAX = MAX (FEMAX, VARNP(LINKF1(K)))
         FEMIN = MIN (FEMIN, VARNP(LINKF1(K)))
  100 CONTINUE

      RETURN
      END
