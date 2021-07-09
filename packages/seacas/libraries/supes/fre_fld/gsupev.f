C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details
C=======================================================================
      SUBROUTINE GSUPEV (STRING)
C=======================================================================

C***********************************************************************

C     *** GSUPEV *** Returns current version number of the supes library
C                    Get SUPEs Version number
C     Parameters:
C     STRING - OUT - string containing versions of supes

C     Version Number Format:
C           n1.n2.n3

C     where n1 is the major version number
C           n2 is the minor version or change capability number
C           n3 is the bug-fix number
C***********************************************************************

      CHARACTER*32 STRING

C     Version Format: 'major.minor.bug_fixl'
      STRING = '4.4.0'

      RETURN
      END

