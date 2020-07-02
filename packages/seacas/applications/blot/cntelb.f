C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C $Log: cntelb.f,v $
C Revision 1.2  2009/03/25 12:36:43  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.1  1994/04/07 19:56:43  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:48:43  gdsjaar
c Added RCS Id and Log to all files
c
C=======================================================================
      SUBROUTINE CNTELB (IELBST, NELBLK, NUMON, NUMSEL)
C=======================================================================

C   --*** CNTELB *** (MESH) Counts selected element blocks
C   --   Written by Amy Gilkey - revised 01/23/87
C   --
C   --CNTELB counts the number of ON element blocks and the number of selected
C   --element blocks.
C   --
C   --Parameters:
C   --   IELBST - IN - the element block status:
C   --      -1 = OFF, 0 = ON, but not selected, 1 = selected
C   --   NELBLK - IN - the number of element blocks
C   --   NUMON - OUT - the number of ON element blocks
C   --   NUMSEL - OUT - the number of selected element blocks

      INTEGER IELBST(NELBLK)

C   --Count the number of ON element blocks

      NUMON = 0
      DO 100 IELB = 1, NELBLK
         IF (IELBST(IELB) .GE. 0) NUMON = NUMON + 1
  100 CONTINUE

C   --Count the number of selected element blocks

      NUMSEL = 0
      DO 110 IELB = 1, NELBLK
         IF (IELBST(IELB) .GT. 0) NUMSEL = NUMSEL + 1
  110 CONTINUE

      RETURN
      END
