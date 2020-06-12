C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C $Log: dbselb.f,v $
C Revision 1.2  2009/03/25 12:36:43  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.1  1994/04/07 19:59:11  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:49:18  gdsjaar
c Added RCS Id and Log to all files
c
C=======================================================================
      SUBROUTINE DBSELB (NELBLK, NUMEL, LENE, INELB, NLISEL, LISEL)
C=======================================================================

C   --*** DBSBEL *** (BLOT) Select elements list of element blocks
C   --   Written by Amy Gilkey - revised 01/05/88
C   --
C   --DBSBEL creates the element block selection array and the element
C   --selection array (by block) given a list of selected element blocks.
C   --
C   --Parameters:
C   --   NELBLK - IN - the number of element blocks
C   --   NUMEL - IN - the number of elements
C   --   LENE - IN - the cumulative element counts by element block
C   --   INELB - IN - the indices of the selected element blocks
C   --   NLISEL - IN/OUT - the number of selected elements for each block
C   --   LISEL - IN/OUT - the indices of the selected elements (by block)

      INTEGER LENE(0:*)
      INTEGER INELB(0:*)
      INTEGER NLISEL(0:*)
      INTEGER LISEL(0:*)

      NLISEL(0) = 0
      CALL INIINT (NELBLK, 0, NLISEL(1))
      DO 100 IX = 1, INELB(0)
         IELB = INELB(IX)
         NLISEL(IELB) = LENE(IELB) - LENE(IELB-1)
  100 CONTINUE
      NLISEL(0) = INELB(0)

      LISEL(0) = 0
      CALL INIINT (NUMEL, 0, LISEL(1))

      DO 120 IELB = 1, NELBLK
         IF (NLISEL(IELB) .GT. 0) THEN
            LISEL(0) = LISEL(0) + NLISEL(IELB)
            DO 110 IEL = LENE(IELB-1)+1, LENE(IELB)
               LISEL(IEL) = IEL
  110       CONTINUE
         END IF
  120 CONTINUE

      RETURN
      END
