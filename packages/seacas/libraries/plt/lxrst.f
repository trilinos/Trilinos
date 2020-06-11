C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C $Id: lxrst.f,v 1.1 1993/07/16 16:46:48 gdsjaar Exp $
C $Log: lxrst.f,v $
C Revision 1.1  1993/07/16 16:46:48  gdsjaar
C Changed plt to library rather than single source file.
C
C=======================================================================
      SUBROUTINE LXRST
      IMPLICIT INTEGER (A-Z)
      CHARACTER*504 ILINE
      COMMON /LXCOM1/ILINE
      COMMON /LXCOM2/JLINE,LXINIT

      LXINIT = 12345
      JLINE = 504
      ILINE(JLINE:JLINE) = CHAR(0)
      RETURN

      END
