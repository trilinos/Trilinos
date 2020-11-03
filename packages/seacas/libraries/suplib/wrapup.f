C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE WRAPUP (PROG)
C=======================================================================

C   --*** WRAPUP *** (ETCLIB) Wrapup program
C   --   Written by Amy Gilkey - revised 04/23/86
C   --
C   --WRAPUP should be called at the end of any program.  It performs
C   --finishing details common to all programs.  Specifically, it
C   --   Gets and displays the CPU time used
C   --
C   --Parameters:
C   --   PROG - IN - the program name

      CHARACTER*(*) PROG

      CHARACTER*8 STR8

C   --Get and display the CPU time used

      CALL EXCPUS (CPUSEC)

      WRITE (STR8, '(F8.2)') CPUSEC
      CALL SQZSTR (STR8, LSTR)
      WRITE (*, 10000) PROG(:LENSTR(PROG)), STR8(:LSTR)
10000  FORMAT (/, 1X, A, ' used ', A, ' seconds of CPU time')

      RETURN
      END
