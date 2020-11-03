C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE FNODES (IFACE, LINKE, NODEF)
C=======================================================================

C   --*** FNODES *** (MESH) Get nodes for element's face
C   --   Written by Amy Gilkey - revised 05/23/86
C   --              Sam Key, 03/01/85
C   --
C   --FNODES returns the nodes for a given face of a hexahedral element.
C   --The 4-tuple sequence defining the face is counter-clockwise looking
C   --into the element.
C   --
C   --Parameters:
C   --   IFACE - IN - the face number
C   --   LINKE - IN - the hexahedral element connectivity
C   --   NODEF - OUT - the nodes in the extracted face

      INTEGER LINKE(8), NODEF(4)

      INTEGER KFACE(4,6)
      SAVE KFACE
C      --KFACE(*,i) - the indices of the 4 nodes in face i

      DATA KFACE / 1, 4, 3, 2,  5, 6, 7, 8,
     &   3, 4, 8, 7,  2, 3, 7, 6,  1, 2, 6, 5,  1, 5, 8, 4 /

      DO 100 K = 1, 4
         NODEF(K) = LINKE(KFACE(K,IFACE))
  100 CONTINUE

      RETURN
      END
