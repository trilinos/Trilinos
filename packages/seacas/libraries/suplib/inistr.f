C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE INISTR (LEN, IFROM, ITO)
C=======================================================================

C   --*** INISTR *** (STRLIB) Initialize all strings in list
C   --   Written by Amy Gilkey - revised 03/15/88
C   --
C   --INISTR initializes all the strings in a list to a specified value.
C   --
C   --Parameters:
C   --   LEN - IN - the number of strings in the list
C   --   IFROM - IN - the initial value
C   --   ITO - OUT - the initialized list

      INTEGER LEN
      CHARACTER*(*) IFROM
      CHARACTER*(*) ITO(*)

      DO 100 I = 1, LEN
         ITO(I) = IFROM
  100 CONTINUE

      RETURN
      END
