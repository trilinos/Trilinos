C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C=======================================================================
      INTEGER FUNCTION IENGRX (IMAXE, IMINE)
C=======================================================================
C$Id: iengrx.f,v 1.2 2009/03/25 12:46:02 gdsjaar Exp $
C$Log: iengrx.f,v $
CRevision 1.2  2009/03/25 12:46:02  gdsjaar
CAdd copyright and license notice to all files.
C
CRevision 1.1.1.1  1990/08/14 16:14:55  gdsjaar
CTesting
C
c Revision 1.1  90/08/14  16:14:54  gdsjaar
c Initial revision
c
c Revision 1.1  90/08/09  13:39:29  gdsjaar
c Initial revision
c

C   --*** IENGX *** (STRLIB) Internal to NUMSTR
C   --   Written by Amy Gilkey - revised 02/14/86
C   --
C   --IENGRX returns the "best" engineering notation exponent for a
C   --minimum and maximum exponent.
C   --
C   --Parameters:
C   --   IMAXE - IN - the maximum exponent
C   --   IMINE - IN - the minimum exponent

      INTEGER IMAXE, IMINE

      IF (IMAXE .GT. 0) THEN
         IENGRX = INT ((IMAXE - 1) / 3) * 3
      ELSE
         IENGRX = INT ((IMAXE - 2) / 3) * 3
      END IF

      RETURN
      END
