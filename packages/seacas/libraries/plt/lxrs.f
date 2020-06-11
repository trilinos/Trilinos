C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C $Id: lxrs.f,v 1.1 1993/07/16 16:46:47 gdsjaar Exp $
C $Log: lxrs.f,v $
C Revision 1.1  1993/07/16 16:46:47  gdsjaar
C Changed plt to library rather than single source file.
C
C=======================================================================
      SUBROUTINE LXRS(IP)
      IMPLICIT INTEGER (A-Z)
      CHARACTER*504 ILINE
      COMMON /LXCOM1/ILINE
      COMMON /LXCOM2/JLINE,LXINIT

      IF (IP.GT.0 .AND. IP.LE.504) THEN
         JLINE = IP

      ELSE
         CALL LXERR('Illegal pointer restoration',3)
      END IF

      RETURN

      END
