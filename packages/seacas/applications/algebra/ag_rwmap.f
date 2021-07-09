C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE RWMAP (NDBIN, NDBOUT, NUMEL, NUMELO, IXELEM,
     &                  MAPEL, NEWIX)
C=======================================================================
C   --*** RWMAP *** (ALGEBRA) Read and write database element order map
C   --   Written by Amy Gilkey - revised 04/28/88
C   --   Modified for EXODUSIIV2 format 8/29/95
C   --
C   --RWMAP reads and writes the element order map to the database.
C   --Deleted elements are removed.
C   --
C   --Parameters:
C   --   NDBIN, NDBOUT - IN - the input and output database file
C   --   NUMEL  - IN - the number of elements
C   --   NUMELO - IN - the number of output elements
C   --   IXELEM - IN - the indices of the output elements (iff NUMELO <> NUMEL)
C   --   IOERR  - OUT - input/output error flag
C   --   MAPEL - SCRATCH - the element order map
C   --   NEWIX - SCRATCH - size = NUMEL (iff NUMELO <> NUMEL)
C   --
C   --Database must be positioned at start of map upon entry;
C   --upon exit at end of map.

      INTEGER IXELEM(*)
      INTEGER MAPEL(*)
      INTEGER NEWIX(*)

      call exgmap (ndbin, mapel, ierr)

      IF ((NUMELO .GT. 0) .AND. (NUMEL .NE. NUMELO)) THEN
        do ix = 1, numelo
          newix(ix) = mapel(ixelem(ix))
        end do
        do ix = 1, numelo
          mapel(ix) = newix(ix)
        end do
      END IF

      call expmap(ndbout, mapel, ierr)

      call exgenm(ndbin, mapel, ierr)

      IF ((NUMELO .GT. 0) .AND. (NUMEL .NE. NUMELO)) THEN
        do ix = 1, numelo
          newix(ix) = mapel(ixelem(ix))
        end do
        do ix = 1, numelo
          mapel(ix) = newix(ix)
        end do
      END IF

      call expenm(ndbout, mapel, ierr)

      RETURN
      END
