C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE FFADDC (CVAL, LINE)
C=======================================================================
C$Id: ffaddc.f,v 1.2 2009/03/25 12:46:01 gdsjaar Exp $
C$Log: ffaddc.f,v $
CRevision 1.2  2009/03/25 12:46:01  gdsjaar
CAdd copyright and license notice to all files.
C
CRevision 1.1.1.1  1990/08/14 16:14:07  gdsjaar
CTesting
C
c Revision 1.1  90/08/14  16:14:05  gdsjaar
c Initial revision
c
c Revision 1.1  90/08/09  13:39:22  gdsjaar
c Initial revision
c

C   --*** FFADDC *** (FFLIB) Add character string to line
C   --   Written by Amy Gilkey - revised 11/16/87
C   --
C   --FFADDC adds a character string to a line.
C   --
C   --Parameters:
C   --   CVAL - IN - the character string to add
C   --   LINE - IN/OUT - the line being built

      CHARACTER*(*) CVAL
      CHARACTER*(*) LINE

      IF (LINE .EQ. ' ') THEN
         I = 0
      ELSE
         I = LENSTR (LINE) + 1
         IF (I .LE. LEN (LINE)) LINE(I:I) = ' '
      END IF
      IF (I .LT. LEN (LINE)) THEN
         LINE(I+1:) = CVAL
      END IF

      RETURN
      END
