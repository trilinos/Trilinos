C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C $Id: cpuerr.f,v 1.1 1993/07/16 16:46:27 gdsjaar Exp $
C $Log: cpuerr.f,v $
C Revision 1.1  1993/07/16 16:46:27  gdsjaar
C Changed plt to library rather than single source file.
C
C=======================================================================
      SUBROUTINE CPUERR(STR,DISP)
      CHARACTER*(*) STR

      CALL SIORPT('CPU',STR,DISP)
      RETURN

      END
