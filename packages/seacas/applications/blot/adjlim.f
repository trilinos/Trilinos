C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE ADJLIM (MSHDEF,
     &   XISSYM, YISSYM, LFTSYM, BOTSYM, XAXSYM, YAXSYM,
     &   SQUARE, RDMESH, ZMMESH)
C=======================================================================

C   --*** ADJLIM *** (MESH) Adjust zoom mesh limits
C   --   Written by Amy Gilkey - revised 05/26/88
C   --
C   --ADJLIM adjusts the zoom window to make a (square) window centered on
C   --the user-requested limits.  If there are symmetric views, the center
C   --is the symmetry axis.  The zoom mesh limits may be adjusted to take
C   --account of the fact that one axis may be twice as long as the other
C   --axis in 2-view plots.
C   --
C   --On a 2-view plot, a check is made to see if the "short" limits
C   --can be made half as long.  This occurs if the user-specified window
C   --can be fit within the given "short" window.
C   --
C   --Parameters:
C   --   MSHDEF - IN - the display modes for all views (none/defined)
C   --      (as in /MSHOPT/)
C   --   XISSYM - IN - true iff X views are symmetric vs copies
C   --      (false if X not divided)
C   --   YISSYM - IN - true iff Y views are symmetric vs copies
C   --      (false if Y not divided)
C   --   LFTSYM - IN - true iff symmetry line to left of main view
C   --   BOTSYM - IN - true iff symmetry line to bottom of main view
C   --   XAXSYM - IN - the vertical symmetry axis
C   --   YAXSYM - IN - the horizontal symmetry axis
C   --   SQUARE - IN - true iff window limits must be square
C   --   RDMESH - IN - the user-requested window limits
C   --      (left, right, bottom, top)
C   --   ZMMESH - OUT - the window limits (may be RDMESH)
C   --      (left, right, bottom, top)

      PARAMETER (KLFT=1, KRGT=2, KBOT=3, KTOP=4)

      CHARACTER*(*) MSHDEF(4)
      LOGICAL XISSYM, YISSYM, LFTSYM, BOTSYM
      LOGICAL SQUARE
      REAL RDMESH(KTOP), ZMMESH(KTOP)

      IF (XISSYM) THEN
         IF (LFTSYM) THEN
            XRNG = RDMESH(KRGT) - XAXSYM
         ELSE
            XRNG = XAXSYM - RDMESH(KLFT)
         END IF
      ELSE
         XRNG = RDMESH(KRGT) - RDMESH(KLFT)
      END IF
      IF (YISSYM) THEN
         IF (BOTSYM) THEN
            YRNG = RDMESH(KTOP) - YAXSYM
         ELSE
            YRNG = YAXSYM - RDMESH(KBOT)
         END IF
      ELSE
         YRNG = RDMESH(KTOP) - RDMESH(KBOT)
      END IF

      IF (SQUARE) THEN
         RNG = MAX (XRNG, YRNG)
         IF ((MSHDEF(1) .NE. 'NONE') .AND. (MSHDEF(3) .EQ. 'NONE')) THEN
            IF (XRNG .LE. 0.51*RNG) RNG = 0.5*RNG
         END IF
      ELSE
         RNG = XRNG
      END IF
      IF (XISSYM) THEN
         IF (LFTSYM) THEN
            ZMMESH(KLFT) = XAXSYM
            ZMMESH(KRGT) = ZMMESH(KLFT) + RNG
         ELSE
            ZMMESH(KRGT) = XAXSYM
            ZMMESH(KLFT) = ZMMESH(KRGT) - RNG
         END IF
      ELSE
         ZMMESH(KLFT) = 0.5 * (RDMESH(KRGT) + RDMESH(KLFT) - RNG)
         ZMMESH(KRGT) = ZMMESH(KLFT) + RNG
      END IF

      IF (SQUARE) THEN
         RNG = MAX (XRNG, YRNG)
         IF ((MSHDEF(4) .NE. 'NONE') .AND. (MSHDEF(3) .EQ. 'NONE')) THEN
            IF (YRNG .LE. 0.51*RNG) RNG = 0.5*RNG
         END IF
      ELSE
         RNG = YRNG
      END IF
      IF (YISSYM) THEN
         IF (BOTSYM) THEN
            ZMMESH(KBOT) = YAXSYM
            ZMMESH(KTOP) = ZMMESH(KBOT) + RNG
         ELSE
            ZMMESH(KTOP) = YAXSYM
            ZMMESH(KBOT) = ZMMESH(KTOP) - RNG
         END IF
      ELSE
         ZMMESH(KBOT) = 0.5 * (RDMESH(KTOP) + RDMESH(KBOT) - RNG)
         ZMMESH(KTOP) = ZMMESH(KBOT) + RNG
      END IF

      RETURN
      END
