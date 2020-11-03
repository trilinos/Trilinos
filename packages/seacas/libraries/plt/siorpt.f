C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE SIORPT(MODULE,MESS,DISP)
      IMPLICIT INTEGER (A-Z)
      CHARACTER*(*) MODULE
      CHARACTER*(*) MESS
      CHARACTER*256 LOCLIN
      CHARACTER*20 ERRORT
      CHARACTER*120 ERRMSG
      CHARACTER*8 ERRMOD

      ERRMOD = MODULE
      ERRMSG = MESS
      CALL CPUDAC(ERRORT)
      CALL CHRTRM(ERRORT,LT)
      CALL CHRTRM(ERRMOD,L)
      IF (L.EQ.0) THEN
         L = 1
      END IF

      IF (DISP.EQ.1) THEN

      ELSE IF (DISP.EQ.2) THEN
         LOCLIN = '#PLT error (warning) in module '//ERRMOD(1:L)//
     *            ' at '//ERRORT(1:LT)
         WRITE (6,10) LOCLIN

   10    FORMAT (1X,A)

         LOCLIN = '#Error message: '//ERRMSG
         WRITE (6,10) LOCLIN

      ELSE IF (DISP.EQ.3) THEN
         LOCLIN = '#PLT error (traceback) in module '//ERRMOD(1:L)//
     *            ' at '//ERRORT(1:LT)
         WRITE (6,10) LOCLIN
         LOCLIN = '#Error message: '//ERRMSG
         WRITE (6,10) LOCLIN

      ELSE IF (DISP.EQ.4) THEN
         LOCLIN = '#PLT error (fatal) in module '//ERRMOD(1:L)//' at '//
     *            ERRORT(1:LT)
         WRITE (6,10) LOCLIN
         LOCLIN = '#Error message: '//ERRMSG
         WRITE (6,10) LOCLIN

      ELSE
         LOCLIN = '#PLT error (fatal) in module '//ERRMOD(1:L)//' at '//
     *            ERRORT(1:LT)
         WRITE (6,10) LOCLIN
         LOCLIN = '#Error message: '//ERRMSG
         WRITE (6,10) LOCLIN
      END IF

      RETURN

      END
