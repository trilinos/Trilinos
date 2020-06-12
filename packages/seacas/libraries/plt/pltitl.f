C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C $Id: pltitl.f,v 1.1 1993/07/16 16:48:30 gdsjaar Exp $
C $Log: pltitl.f,v $
C Revision 1.1  1993/07/16 16:48:30  gdsjaar
C Changed plt to library rather than single source file.
C
C=======================================================================
      INTEGER FUNCTION PLTITL(REALN)

      IF (PLTFRC(REALN).EQ.0.) THEN
         PLTITL = INT(REALN)
         RETURN

      END IF

      AREAL = ABS(REALN)
      PLTITL = INT(AREAL)
      IF (REALN.LT.0.) THEN
         PLTITL = -PLTITL - 1
      END IF

      RETURN

      END
