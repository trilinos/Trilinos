C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C $Id: lxcln.f,v 1.1 1993/07/16 16:46:37 gdsjaar Exp $
C $Log: lxcln.f,v $
C Revision 1.1  1993/07/16 16:46:37  gdsjaar
C Changed plt to library rather than single source file.
C
C=======================================================================
      SUBROUTINE LXCLN
      IMPLICIT INTEGER (A-Z)
      CHARACTER*504 ILINE
      COMMON /LXCOM1/ILINE
      COMMON /LXCOM2/JLINE,LXINIT

      K = INDEX(ILINE(JLINE:),CHAR(0))
      IF (K.EQ.0) THEN
         JLINE = 504

      ELSE
         JLINE = MIN(JLINE+K,504)
      END IF

      RETURN

      END
