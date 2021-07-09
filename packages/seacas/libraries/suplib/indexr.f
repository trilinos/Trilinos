C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      INTEGER FUNCTION INDEXR (STRING, CHAR)
C=======================================================================
C   --*** INDEXR *** (STRLIB) Find last occurrence of CHAR in STRING
C   --
C   --INDEXR returns the last position of CHAR in STRING
C   --if not found, returns 0
C   --

      CHARACTER*(*) STRING
      CHARACTER*1   CHAR

      DO i=len(string),1,-1
        if (string(i:i) .eq. char) then
          indexr = i
          return
        end if
      end do

      indexr = 0
      return
      end
