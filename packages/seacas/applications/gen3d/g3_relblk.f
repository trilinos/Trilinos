C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE RELBLK (IELB, NUMELB, NUMLNK, NUMATR,
     &   LINK, ATRIB, LNKTMP, *)
C=======================================================================

C   --*** RELBLK *** (GEN3D) Read database element block misc.
C   --   Written by Amy Gilkey - revised 02/27/86
C   --
C   --RELBLK reads the element block connectivity and attribute information
C   --from the database.  An error message is displayed if the end of file
C   --is read.
C   --
C   --Parameters:
C   --   IELB - IN - the element block number (for errors)
C   --   NUMELB - IN - the number of elements in the block
C   --   NUMLNK - IN - the number of nodes per element
C   --   NUMATR - IN - the number of attributes
C   --   LINK - OUT - the element connectivity for this block (4 nodes always)
C   --   ATRIB - OUT - the attributes for this block (packed)
C   --   * - return statement if end of file or read error
C   --
C   --Common Variables:
C   --   Uses NDBIN of /DBASE/

      include 'exodusII.inc'
      INCLUDE 'g3_dbase.blk'

      INTEGER LINK(4,NUMELB)
      INTEGER LNKTMP(NUMLNK, NUMELB)
      REAL ATRIB(*)

      CHARACTER*5 STRA

      if (numlnk .eq. 4) then
        call exgelc(ndbin, ielb, link, ierr)
        if (ierr .ne. 0) goto 10
      else
        call exgelc(ndbin, ielb, lnktmp, ierr)
        if (ierr .ne. 0) goto 10
        do 30 i=1, numelb
          do 20 j=1, numlnk
            link(j,i) = lnktmp(j,i)
 20       continue
 30     continue
      end if

      if (numatr .gt. 0) then
        call exgeat(ndbin, ielb, atrib, ierr)
        if (ierr .ne. 0) goto 10
      end if

      RETURN

   10 CONTINUE
      CALL INTSTR (1, 0, IELB, STRA, LSTRA)
      CALL PRTERR ('FATAL', 'Reading BLOCK ' // STRA(:LSTRA))
      RETURN 1
      END
