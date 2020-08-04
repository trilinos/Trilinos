C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE SETVW (IVIEW, *)
C=======================================================================

C   --*** SETVW *** (MESH) Set up viewport for view
C   --   Written by Amy Gilkey - revised 02/12/85
C   --
C   --SETVW sets up the the viewport for the specified view.
C   --
C   --Parameters:
C   --   IVIEW - IN - the view number
C   --   * - return statement the cancel function is active
C   --
C   --Common Variables:
C   --   Uses DVIEW, WVIEW of /LAYOUT/

      PARAMETER (KLFT=1, KRGT=2, KBOT=3, KTOP=4)

      COMMON /LAYOUD/ DVIEW(KTOP,4), WVIEW(KTOP,4)

      LOGICAL MPVIEW, MPORT2, LDUM
      LOGICAL GRABRT

      IF (GRABRT ()) RETURN 1

      LDUM = MPVIEW (DVIEW(KLFT,IVIEW), DVIEW(KRGT,IVIEW),
     &   DVIEW(KBOT,IVIEW), DVIEW(KTOP,IVIEW))
      LDUM = MPORT2 (WVIEW(KLFT,IVIEW), WVIEW(KRGT,IVIEW),
     &   WVIEW(KBOT,IVIEW), WVIEW(KTOP,IVIEW))

      RETURN
      END
