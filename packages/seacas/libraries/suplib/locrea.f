C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C=======================================================================
      INTEGER FUNCTION LOCREA (VALU, NVALUS, VALUS)
C=======================================================================
C$Id: locrea.f,v 1.2 2009/03/25 12:46:02 gdsjaar Exp $
C$Log: locrea.f,v $
CRevision 1.2  2009/03/25 12:46:02  gdsjaar
CAdd copyright and license notice to all files.
C
CRevision 1.1.1.1  1990/08/14 16:15:21  gdsjaar
CTesting
C
c Revision 1.1  90/08/14  16:15:20  gdsjaar
c Initial revision
c
c Revision 1.1  90/08/09  13:39:34  gdsjaar
c Initial revision
c

C   --*** LOCREA *** (ETCLIB) Find closest real value
C   --   Written by Amy Gilkey - revised 11/03/87
C   --
C   --LOCREA returns the closest real value to the given value in a list of
C   --real values (which may not be ordered).
C   --
C   --Parameters:
C   --   VALU - the value to be searched for
C   --   NVALUS - the number of values in the list
C   --   VALUS - the list of values

      REAL VALUS(*)

      DIF = ABS (VALUS(1) - VALU)
      IX = 1
      DO 10 I = 2, NVALUS
         DIFI = ABS (VALUS(I) - VALU)
         IF (DIF .GT. DIFI) THEN
            DIF = DIFI
            IX = I
         END IF
   10 CONTINUE

      LOCREA = IX

      RETURN
      END
