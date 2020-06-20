C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE NUM2IX (NSETS, NINSET, IXSET)
C=======================================================================
C$Id: num2ix.f,v 1.2 2009/03/25 12:46:02 gdsjaar Exp $
C$Log: num2ix.f,v $
CRevision 1.2  2009/03/25 12:46:02  gdsjaar
CAdd copyright and license notice to all files.
C
CRevision 1.1.1.1  1990/08/14 16:15:45  gdsjaar
CTesting
C
c Revision 1.1  90/08/14  16:15:44  gdsjaar
c Initial revision
c
c Revision 1.1  90/08/09  13:39:39  gdsjaar
c Initial revision
c

C   --*** NUM2IX *** (ETCLIB) Change number-in-set to set indices
C   --   Written by Amy Gilkey - revised 10/23/87
C   --
C   --NUM2IX creates the set indices given the number in each set.
C   --
C   --Parameters:
C   --   NSETS - IN - the number of sets
C   --   NINSET - IN - the number in each set
C   --   IXSET - OUT - the set indices; set i = IXSET(i-1)+1 .. IXSET(i)

      INTEGER NINSET(*)
      INTEGER IXSET(0:*)

      IXSET(0) = 0
      DO 10 I = 1, NSETS
         IXSET(I) = IXSET(I-1) + NINSET(I)
   10 CONTINUE

      RETURN
      END
