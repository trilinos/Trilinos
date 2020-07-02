C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C $Id: cpuqa.f,v 1.1 1993/07/16 16:46:31 gdsjaar Exp $
C $Log: cpuqa.f,v $
C Revision 1.1  1993/07/16 16:46:31  gdsjaar
C Changed plt to library rather than single source file.
C
C=======================================================================
      SUBROUTINE CPUQA(DATE1,TIME,USER,JOBID,DIV,ENUM,CASENO,CLASS)
      CHARACTER*(*) DATE1
      CHARACTER*(*) TIME
      CHARACTER*(*) USER
      CHARACTER*(*) JOBID
      CHARACTER*(*) DIV
      CHARACTER*(*) ENUM
      CHARACTER*(*) CASENO
      CHARACTER*(*) CLASS

      DATE1 = ' '
      TIME = ' '
      JOBID = ' '
      USER = ' '
      DIV = ' '
      ENUM = ' '
      CASENO = ' '
      CLASS = ' '
      RETURN

      END
