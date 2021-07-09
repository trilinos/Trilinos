C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C========================================================================
      SUBROUTINE AVG(IGLND,INVCN,MAXLN,INVLEN,SOLEA,SOLENA,ITT,iblk)

C************************************************************************

C Subroutine AVG provides for translating nodal values of element
C variables back to the element centroids for the special case where
C too few elements can be associated with a node. Element variable
C data is simply averaged at that node.

C Called by ELTON1

C************************************************************************

C  IGLND  INT   The global node number
C  INVCN  INT   The inverse connectivity (1:maxln,1:numnda)
C  MAXLN  INT   The maximum nomber of elements connected to any node
C  INVLEN INT   The number of elements connected to this node
C  SOLEA  REAL  Element variables (1:numeba,1:nvarel)
C  SOLENA REAL  Element variables at nodes" (1:nodesa,1:nvarel)
C  NDLSTA INT   The array that identifies the local element block node
C               number with the global mesh node number (1:numnda)
C  ITT    INT   Truth table
C  iblk   INT   Block number being processed (not block ID)

C************************************************************************

      include 'aexds1.blk'
      include 'amesh.blk'
      include 'ebbyeb.blk'

      DIMENSION INVCN(MAXLN,*),SOLEA(NUMEBA,*),
     &          SOLENA(NODESA,NVAREL), ITT(NVAREL,*)

C************************************************************************

      DO 10 IVAR = 1, NVAREL
        IF (ITT(IVAR,iblk) .EQ. 0)GO TO 10
        SUM = 0.
        DO 20 J = 1, INVLEN
          SUM = SUM + SOLEA(INVCN(J,IGLND),IVAR)
   20   CONTINUE
        SOLENA(IGLND,IVAR) = SUM / INVLEN
   10 CONTINUE
      RETURN
      END
