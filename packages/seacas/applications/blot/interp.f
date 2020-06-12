C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C $Log: interp.f,v $
C Revision 1.2  2009/03/25 12:36:45  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.1  1994/04/07 20:03:56  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:52:45  gdsjaar
c Added RCS Id and Log to all files
c
C=======================================================================
      LOGICAL FUNCTION INTERP_BL (CNTR, F1, F2, PSI)
C=======================================================================

C   --*** INTERP_BL *** (DETOUR) Compute interception point
C   --   Written by Amy Gilkey - revised 11/21/85
C   --   D. P. Flanagan, 3/25/83
C   --
C   --INTERP_BL tests if a contour value falls within an interval.
C   --A degenerate interval fails regardless of the contour value. If the
C   --test is passed, the nomalized interval coordinate of the contour
C   --value is computed.  The returned function value is true only if the
C   --interval test is passed.
C   --
C   --Parameters:
C   --   CNTR - IN - the contour value
C   --   F1 - IN - the interval origin value
C   --   F2 - IN - the interval terminus value
C   --   PSI - OUT - normalized interval coordinate

C   --Test if coordinate lies within the interval

      INTERP_BL = (F2 .NE. F1) .AND.
     &   (((F1 .LE. CNTR) .AND. (CNTR .LE. F2))
     &   .OR. ((F2 .LE. CNTR) .AND. (CNTR .LE. F1)))

C   --Compute normalized interval coordinate

      IF (INTERP_BL) PSI = (CNTR - F1) / (F2 - F1)
      RETURN
      END
