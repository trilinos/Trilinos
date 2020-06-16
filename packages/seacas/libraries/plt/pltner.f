C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C $Id: pltner.f,v 1.1 1993/07/16 16:49:00 gdsjaar Exp $
C $Log: pltner.f,v $
C Revision 1.1  1993/07/16 16:49:00  gdsjaar
C Changed plt to library rather than single source file.
C
C=======================================================================
      LOGICAL FUNCTION PLTNER(X,Y)

      IF (ABS(X-Y).LT..001) THEN
         PLTNER = .TRUE.
         RETURN

      ELSE
         PLTNER = .FALSE.
         RETURN

      END IF

      END
