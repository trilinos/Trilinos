C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C $Id: lxerr.f,v 1.1 1993/07/16 16:46:38 gdsjaar Exp $
C $Log: lxerr.f,v $
C Revision 1.1  1993/07/16 16:46:38  gdsjaar
C Changed plt to library rather than single source file.
C
C=======================================================================
      SUBROUTINE LXERR(MSG,DISP)
      CHARACTER*(*) MSG
      INTEGER DISP
      CHARACTER*80 LOCMSG

      LOCMSG = MSG
      IF (DISP.GE.2) THEN
         CALL LXCLN
      END IF

      CALL SIORPT('LEX',LOCMSG,DISP)
      RETURN

      END
