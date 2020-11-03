C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details
      SUBROUTINE MYDEL (NAME1, DICT, DPOINT, LDICT, NNAMES, VOID,
     *   LVOID, NVOIDS, CHRCOL, LASTER, MYLOC, MYCLOC)

      IMPLICIT INTEGER (A-Z)
      INCLUDE 'params.inc'

C     This routine removes a name from the dictionary and returns the
C     available space to the void table.

C***********************************************************************

C     NAME1    Name to be deleted
               CHARACTER*8 NAME1
C     DICT     Dictionary name table
C     DPOINT   Dictionary pointer table
C     LDICT    Dimension of dictionary tables
C     NNAMES   Number of names in the dictionary
               CHARACTER*8 DICT(LDICT,CHRCOL)
               DIMENSION DPOINT(LDICT,CHRCOL,3), NNAMES(2)
C     VOID     Void table
C     LVOID    Dimension of void table
C     NVOIDS   Number of voids
               DIMENSION VOID(LVOID,CHRCOL,2), NVOIDS(2)
C     CHRCOL   Number of column for character names.
C     LASTER   Error return

C***********************************************************************

C     FIND NAME1 IN DICTIONARY.

      CALL MYFIND (NAME1, DICT, DPOINT, LDICT, NNAMES,
     *   CHRCOL, LASTER, ROW)
      IF (LASTER .NE. SUCESS) RETURN

      LOC = DPOINT(ROW,CHRCOL,1)
      LEN = DPOINT(ROW,CHRCOL,2)

C     DELETE DICTIONARY ENTRY.

      CALL SHFTC (DICT(1,CHRCOL), CHRCOL*LDICT, ROW+1,
     *   NNAMES(CHRCOL), 1)
      CALL SHFTI (DPOINT(1,CHRCOL,1), CHRCOL*LDICT, 3, ROW+1,
     *   NNAMES(CHRCOL), 1)
      NNAMES(CHRCOL) = NNAMES(CHRCOL) - 1
      IF (LEN .LE. 0) RETURN

C ... Using malloc/free -- let system manage void space. Return
C     memory to system via 'free'.  The value given to memret
C     is a flag which tells the system that this is a 'safe' free
C     which should actually execute. (Major Kludge...)
      LASTER = SUCESS
      memret = -999
      if (chrcol .eq. 1) then
        oldadr = loc+myloc-1
      else
        oldadr = loc+mycloc-1
      end if
      call exmemy(-len, oldadr, memret)
      if (memret .lt. 0 .or. memret .gt. len) then
        write (*,*) 'ERROR in MYDEL', memret, len
         laster = ilblk
         return
      end if

      RETURN
      END
