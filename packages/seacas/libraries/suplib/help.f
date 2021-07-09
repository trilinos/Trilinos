C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      LOGICAL FUNCTION HELP (PROG, TYPE, OPTION)
C=======================================================================
      CHARACTER*(*) PROG
      CHARACTER*(*) TYPE
      CHARACTER*(*) OPTION

      HELP = .FALSE.

      WRITE (*, 10)
     &  'Help for ', PROG(:LENSTR(PROG)), ' is not available'
      RETURN
 10   FORMAT (/, 1X, 5A)
      END

