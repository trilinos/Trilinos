C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE ININPF (ilow, ihigh, MAXNPF, NPFS)
C=======================================================================

C   --*** ININPF *** (MESH) Initialize NPFS array
C   --   Written by Amy Gilkey - revised 10/22/87
C   --              Sam Key, 06/01/85
C   --
C   --ININPF initializes the NPFS array by setting the length for all nodes
C   --to zero.
C   --
C   --Parameters:
C   --   NUMNP - IN - the number of nodes
C   --   MAXNPF - IN - the maximum length of the NPFS entry
C   --   NPFS - OUT - the list of unmatched faces containing a node;
C   --      (0,i) = the length of the list

      INTEGER NPFS(*)
      do i=ilow, ihigh
        npfs(i) = 0
      end do

      RETURN
      END
