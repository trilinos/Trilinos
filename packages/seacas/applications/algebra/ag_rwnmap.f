C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE RWNMAP (NDBIN, NDBOUT, NUMNP, NUMNPO, IXNODE,
     &                   MAPND, NEWIX)
C=======================================================================
C   --*** RWNMAP *** (ALGEBRA) Read and write database node number map
C   --
C   --
C   --Parameters:
C   --   NDBIN, NDBOUT - IN - the input and output database file
C   --   NUMNP  - IN - the number of nodes
C   --   NUMNPO - IN - the number of output nodes
C   --   IXNODE - IN - the indices of the output nodes (iff NUMNPO <> NUMNP)
C   --   IOERR  - OUT - input/output error flag
C   --   MAPND - SCRATCH - the node number map
C   --   NEWIX - SCRATCH - size = NUMNP (iff NUMNPO <> NUMNP)

      INTEGER IXNODE(*)
      INTEGER MAPND(*)
      INTEGER NEWIX(*)

      call exgnnm(ndbin, mapnd, ierr)

      IF ((NUMNPO .GT. 0) .AND. (NUMNP .NE. NUMNPO)) THEN
        do ix = 1, numnpo
          newix(ix) = mapnd(ixnode(ix))
        end do
        do ix = 1, numnpo
          mapnd(ix) = newix(ix)
        end do

      END IF

      call expnnm(ndbout, mapnd, ierr)

      RETURN
      END
