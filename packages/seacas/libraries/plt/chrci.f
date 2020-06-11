C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C $Id: chrci.f,v 1.1 1993/07/16 16:46:16 gdsjaar Exp $
C $Log: chrci.f,v $
C Revision 1.1  1993/07/16 16:46:16  gdsjaar
C Changed plt to library rather than single source file.
C
C=======================================================================
      LOGICAL FUNCTION CHRCI(LINE,INTE)
      CHARACTER LINE* (*),FORM*10,CL*2

      CALL CHRTRM(LINE,LL)
      CALL CHRIC(LL,CL,NL)
      FORM = '(i'//CL(1:NL)//')'
      READ (LINE(1:LL),FORM,ERR=10) INTE
      CHRCI = .TRUE.
      RETURN

   10 CHRCI = .FALSE.
      RETURN

      END
