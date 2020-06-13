C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE MINMXL (NPTS, PTOK, PTS, VMIN, VMAX)
C=======================================================================
C$Id: minmxl.f,v 1.2 2009/03/25 12:46:02 gdsjaar Exp $
C$Log: minmxl.f,v $
CRevision 1.2  2009/03/25 12:46:02  gdsjaar
CAdd copyright and license notice to all files.
C
CRevision 1.1.1.1  1990/08/14 16:15:39  gdsjaar
CTesting
C
c Revision 1.1  90/08/14  16:15:38  gdsjaar
c Initial revision
c
c Revision 1.1  90/08/09  13:39:38  gdsjaar
c Initial revision
c

C   --*** MINMXL *** (ETCLIB) Calculate min/max value in selected list
C   --   Written by Amy Gilkey - revised 11/10/87
C   --
C   --MINMXL calculates the minimum and maximum of the data selected by
C   --a logical value.
C   --
C   --Parameters:
C   --   NPTS - IN - the number of points
C   --   ISOK - IN - use point i iff ISOK(i)
C   --   PTS - IN - the points
C   --   VMIN, VMAX - OUT - the maximum and maximum value of the points

      INTEGER NPTS
      LOGICAL PTOK(*)
      REAL PTS(*)
      REAL VMIN, VMAX

      VMIN =  1.0E+30
      VMAX = -1.0E+30
      DO 10 I = 1, NPTS
         IF (PTOK(I)) THEN
            VMIN = MIN (VMIN, PTS(I))
            VMAX = MAX (VMAX, PTS(I))
         END IF
   10 CONTINUE

      RETURN
      END
