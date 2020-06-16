C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C $Id: memini.f,v 1.1 1993/07/16 16:47:00 gdsjaar Exp $
C $Log: memini.f,v $
C Revision 1.1  1993/07/16 16:47:00  gdsjaar
C Changed plt to library rather than single source file.
C
C=======================================================================
      SUBROUTINE MEMINI(MEMRY,LENGTH)
      INTEGER MEMRY(*)
      INTEGER LENGTH

      DO 2630 I = 1,LENGTH
         MEMRY(1) = LENGTH
 2630 CONTINUE
      RETURN

      END
