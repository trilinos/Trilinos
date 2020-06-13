C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C $Log: make2b.f,v $
C Revision 1.2  2009/03/25 12:36:45  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.1  1994/04/07 20:04:48  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:53:17  gdsjaar
c Added RCS Id and Log to all files
c
C=======================================================================
      SUBROUTINE MAKE2B (NELBLK, LENE, IE2ELB)
C=======================================================================

C   --*** MAKE2B *** (BLOT) Create element to element block index
C   --   Written by Amy Gilkey - revised 09/08/87
C   --
C   --MAKE2B creates a list of the element block for each element.
C   --
C   --Parameters:
C   --   NELBLK - IN - the number of element blocks to read
C   --   LENE - IN - the cumulative element counts by element block
C   --   IE2ELB - OUT - the element block for each element
C   --      filled with the element block number, not ID

      INTEGER NELBLK
      INTEGER LENE(0:NELBLK)
      INTEGER IE2ELB(*)

      DO 110 IELB = 1, NELBLK
         DO 100 IEL = LENE(IELB-1)+1, LENE(IELB)
            IE2ELB(IEL) = IELB
  100    CONTINUE
  110 CONTINUE

      RETURN
      END
