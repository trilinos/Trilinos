C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE FFADDV (IRANGE, LINE)
C=======================================================================
C$Id: ffaddv.f,v 1.2 2009/03/25 12:46:02 gdsjaar Exp $
C$Log: ffaddv.f,v $
CRevision 1.2  2009/03/25 12:46:02  gdsjaar
CAdd copyright and license notice to all files.
C
CRevision 1.1.1.1  1990/08/14 16:14:18  gdsjaar
CTesting
C
c Revision 1.1  90/08/14  16:14:17  gdsjaar
c Initial revision
c
c Revision 1.1  90/08/09  13:39:23  gdsjaar
c Initial revision
c

C   --*** FFADDV *** (FFLIB) Add integer range to line
C   --   Written by Amy Gilkey - revised 11/16/87
C   --
C   --FFADDV adds an integer range (as a character string) to a line.
C   --
C   --Parameters:
C   --   IRANGE - IN - the integer range to add;
C   --      IRANGE(1) TO IRANGE(2) BY IRANGE(3)
C   --   LINE - IN/OUT - the line being built

      INTEGER IRANGE(3)
      CHARACTER*(*) LINE

      CALL FFADDI (IRANGE(1), LINE)
      IF (IRANGE(2) .NE. IRANGE(1)) THEN
         CALL FFADDC ('TO', LINE)
         CALL FFADDI (IRANGE(2), LINE)
      END IF
      IF (IRANGE(3) .NE. 1) THEN
         CALL FFADDC ('BY', LINE)
         CALL FFADDI (IRANGE(3), LINE)
      END IF

      RETURN
      END
