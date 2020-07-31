C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE TNODES (IFACE, LINKE, NODEF)
C=======================================================================

C   --*** TNODES *** (MESH) Get nodes for tet element's face
C   --
C   --TNODES returns the nodes for a given face of a tetrahedral element.
C   --The 3-tuple sequence defining the face is counter-clockwise looking
C   --into the element.
C   --
C   --Parameters:
C   --   IFACE - IN - the face number
C   --   LINKE - IN - the tetrahedral element connectivity
C   --   NODEF - OUT - the nodes in the extracted face

      INTEGER LINKE(6), NODEF(4)

      INTEGER KFACE(3,4)
      SAVE KFACE
C      --KFACE(*,i) - the indices of the 3 nodes in face i

      DATA KFACE /
     $     1, 2, 4,
     $     2, 3, 4,
     &     1, 4, 3,
     $     1, 3, 2 /

      DO 100 K = 1, 3
         NODEF(K) = LINKE(KFACE(K,IFACE))
  100 CONTINUE
      NODEF(4) = LINKE(KFACE(1,IFACE))
      RETURN
      END
