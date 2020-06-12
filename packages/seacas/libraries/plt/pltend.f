C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C $Id: pltend.f,v 1.1 1993/07/16 16:48:06 gdsjaar Exp $
C $Log: pltend.f,v $
C Revision 1.1  1993/07/16 16:48:06  gdsjaar
C Changed plt to library rather than single source file.
C
C=======================================================================
      SUBROUTINE PLTEND

      CALL VDFRAM(1)
      CALL VDTERM
      RETURN

      END
