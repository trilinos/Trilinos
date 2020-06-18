C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C $Id: pltiqd.f,v 1.1 1993/07/16 16:48:28 gdsjaar Exp $
C $Log: pltiqd.f,v $
C Revision 1.1  1993/07/16 16:48:28  gdsjaar
C Changed plt to library rather than single source file.
C
C=======================================================================
      SUBROUTINE PLTIQD(ARRAY)
      DIMENSION ARRAY(*)

      DO 2040 I = 1,23
         CALL VDIQDC(I,ARRAY(I))
 2040 CONTINUE
      RETURN

      END
