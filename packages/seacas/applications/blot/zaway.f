C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      LOGICAL FUNCTION ZAWAY (NLNKF, LINKF1, XN, YN, ZN, HIDENP)
C=======================================================================

C   --*** ZAWAY *** (MESH) Determine if face is hidden
C   --   Written by Amy Gilkey - revised 10/22/87
C   --              Sam Key, 03/01/85
C   --
C   --ZAWAY determines if a face on a 3D surface is hidden.
C   --It is hidden if and only if its outward normal points "into" the
C   --plotting surface.  Nodes within the face are marked as
C   --visible if the face is visible.
C   --
C   --Parameters:
C   --   NLNKF - IN - the number of nodes per face
C   --   LINKF1 - IN - the connectivity for the face
C   --   XN, YN, ZN - IN - the nodal coordinates; the Z-coordinate is
C   --      pointing towards the viewer (out of the plotting plane)
C   --   HIDENP - IN/OUT - node status (as in HIDDEN)

      PARAMETER (KNVIS=0, KNFOVR=10, KNHID=100)

      INTEGER LINKF1(NLNKF)
      REAL XN(*), YN(*), ZN(*)
      INTEGER HIDENP(*)

C   --Form X and Y components of diagonal vectors of face
C   --(calculate from midpoint of each side)

      AX = (XN(LINKF1(3)) + XN(LINKF1(4)))
     &   - (XN(LINKF1(1)) + XN(LINKF1(2)))
      BX = (XN(LINKF1(4)) + XN(LINKF1(1)))
     &   - (XN(LINKF1(2)) + XN(LINKF1(3)))
      AY = (YN(LINKF1(3)) + YN(LINKF1(4)))
     &   - (YN(LINKF1(1)) + YN(LINKF1(2)))
      BY = (YN(LINKF1(4)) + YN(LINKF1(1)))
     &   - (YN(LINKF1(2)) + YN(LINKF1(3)))

C   --Form Z component of normal vector to corner, and make node
C   --visible if normal points forward

      ZAWAY = (AX*BY .LE. BX*AY)

      IF (.NOT. ZAWAY) THEN
         DO 100 ILINK = 1, 4
            HIDENP(LINKF1(ILINK)) = KNVIS
  100    CONTINUE
      END IF

      RETURN
      END
