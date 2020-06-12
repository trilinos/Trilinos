C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C $Log: spdsg1.f,v $
C Revision 1.2  2009/03/25 12:36:48  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.1  1994/04/07 20:14:10  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:58:15  gdsjaar
c Added RCS Id and Log to all files
c
C=======================================================================
      SUBROUTINE SPDSG1 (NNENUM, NENUM, IE2ELB, ISEVOK, NEL, ISEGEL)
C=======================================================================

C   --*** SPDSG1 *** (SPLOT) Find all defined element values
C   --   Written by Amy Gilkey - revised 11/05/87
C   --
C   --SPDSG1 returns a list of the defined element indices for an
C   --element variable.
C   --
C   --Parameters:
C   --   NNENUM - IN - the number of selected nodes/elements
C   --   NENUM - IN - the selected nodes/elements
C   --   IE2ELB - IN - the element block for each element
C   --   ISEVOK - IN - the element block variable truth table;
C   --      variable of block j exists iff ISEVOK(j)
C   --   NEL - OUT - the number of defined elements
C   --   ISEGEL - OUT - the NENUM indices of the defined elements;
C   --      ISEGEL(0) = the number of defined elements

      INTEGER NENUM(*)
      INTEGER IE2ELB(*)
      LOGICAL ISEVOK(*)
      INTEGER ISEGEL(0:*)

      NEL = 0
      DO 100 I = 1, NNENUM
         IF (ISEVOK(IE2ELB(NENUM(I)))) THEN
            NEL = NEL + 1
            ISEGEL(NEL) = I
         END IF
  100 CONTINUE
      ISEGEL(0) = NEL

      RETURN
      END
