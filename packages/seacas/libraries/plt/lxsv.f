C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C $Id: lxsv.f,v 1.1 1993/07/16 16:46:54 gdsjaar Exp $
C $Log: lxsv.f,v $
C Revision 1.1  1993/07/16 16:46:54  gdsjaar
C Changed plt to library rather than single source file.
C
C=======================================================================
      INTEGER FUNCTION LXSV()
      IMPLICIT INTEGER (A-Z)
      CHARACTER*504 ILINE
      COMMON /LXCOM1/ILINE
      COMMON /LXCOM2/JLINE,LXINIT

      LXSV = JLINE
      RETURN

      END
