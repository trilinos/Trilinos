C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

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
