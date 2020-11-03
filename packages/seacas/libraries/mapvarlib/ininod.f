C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

      SUBROUTINE ININOD(SOLNB,IVAR,TIMES,ISTP,IDBLK,NDLSTB,XB,YB,ZB,
     &                  SN)

C *********************************************************************

C ININOD initializes nodal variable values based on TIME, ELEMENT
C BLOCK, VARIABLE NAME, COORDINATE, etc. By default, nodal variable
C values are set  to zero. It is intended that the user rewrite this
C subroutine to provide values that are appropriate to the problem
C being solved. This is the preferred method to handle nodal variable
C assignment for recipient mesh nodes that lie outside the boundary
C of the donor mesh.

C Called by INTRPN, SINTPN

C *********************************************************************

C SOLNB   REAL  Array of nodal variable values
C               (1:nodesb,1:nvarnp,1:ntimes)
C IVAR    INT   Position of variable in SOLNB array
C TIMES   REAL  Array of times (1:ntimes)
C ISTP    INT   Position in TIMES array, time step being processed
C IDBLK   INT   The element block I. D.
C NDLSTB  INT   Array of nodes that belong to this element block
C               (1:numndb)
C XB      REAL  Array of X-coordinates (1:nodesb)
C YB      REAL  Array of Y-coordinates (1:nodesb)
C ZB      REAL  Array of Z-coordinates (1:nodesb)
C SN      REAL  Dummy array to store values from MESH-B

C *********************************************************************

      include 'exodusII.inc'

      include 'aexds1.blk'
      include 'bmesh.blk'
      include 'ebbyeb.blk'
      include 'ex2tp.blk'
      include 'inival.blk'

      DIMENSION SOLNB(NODESB,NVARNP), TIMES(*), NDLSTB(*)
      DIMENSION XB(*), YB(*), ZB(*)
      DIMENSION SN(*)

C *********************************************************************

C Code to help you find some potentially useful stuff
C The actual time (real number)
C     TIME = TIMES(ISTP)

C The pointer into VARNAM to get the variable name being processed
C     INAM = IVAR + NVARGP + NVAREL

C The name of the variable (character) being processed
C     NAME = NAMVAR(INAM)

C INOD = NDLSTB(I)

C Set value of node-NDLSTB(I); variable-IVAR to 0. by default.
C User to replace this with whatever code he wishes.

C ... Note: The exgnv step assumes that variable 'ivar' on the mesh
C           database is the same variable as 'ivar' on the results
C           database. This may or may not be the case....

      call exgvp(ntp3ex, 'N', nvar, ierr)
      if (nvar .ge. ivar) then
        CALL EXGNV(NTP3EX,ISTP,IVAR,NODESB,SN,IERR)
        DO 10 I = 1, NUMNDB
          INODE = NDLSTB(I)
          SOLNB(INODE,IVAR) = SN(INODE)
 10     CONTINUE
      else
        DO 20 I = 1, NUMNDB
          INODE = NDLSTB(I)
          SOLNB(INODE,IVAR) = VALINI
 20     CONTINUE
      end if
      RETURN
      END
