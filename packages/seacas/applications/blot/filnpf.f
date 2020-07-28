C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE FILNPF (NLNKF, LINKF1, NFACES, MAXNPF, NPFS, NOVER,
     *  NUMNP)
C=======================================================================

C   --*** FILNPF *** (MESH) Point to face from NPFS
C   --   Written by Amy Gilkey - revised 10/27/87
C   --              Sam Key, 06/01/85
C   --
C   --FILNPF puts a pointer to the given face into each of the face's node
C   --NPFS array.
C   --
C   --Parameters:
C   --   NLNKF - IN - the number of nodes per face
C   --   LINKF1 - IN - the connectivity for the face
C   --   NFACES - IN - the face number
C   --   MAXNPF - IN - the maximum length of the NPFS entry
C   --   NPFS - IN/OUT - the list of unmatched faces containing a node;
C   --      (0,i) = the length of the list
C   --   NOVER - IN/OUT - the number of overrun errors

      include 'minmax.blk'

      INTEGER LINKF1(NLNKF)
      INTEGER NPFS(NUMNP,0:MAXNPF)

      DO 100 ILINK = 1, NLNKF
         INF = LINKF1(ILINK)
         IF (NPFS(INF,0) .LT. MAXNPF) THEN
            L = NPFS(INF,0) + 1
            NPFS(INF,L) = NFACES
            NPFS(INF,0) = L
            minnod = min(minnod, inf)
            maxnod = max(maxnod, inf)
         ELSE
            NOVER = NOVER + 1
         END IF
  100 CONTINUE

      RETURN
      END
