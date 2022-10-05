C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details
      SUBROUTINE MXFREE (OFFSET, DICT, DPOINT, LDICT, NNAMES, CHRCOL,
     *  MYLOC)

      IMPLICIT INTEGER (A-Z)
      INCLUDE 'params.inc'

C     This routine frees all memory allocated by SUPES.
C     All addresses are invalid at this point forward.

C***********************************************************************

C     OFFSET   Offset to internal reference vector
C     DICT     Dictionary name table
C     DPOINT   Dictionary pointer table
C     LDICT    Dimension of dictionary tables
C     NNAMES   Number of names in the dictionary
               CHARACTER*8 DICT(LDICT,CHRCOL)
               DIMENSION DPOINT(LDICT,CHRCOL,3), NNAMES(2)
C     CHRCOL   Number of column for character names.

C***********************************************************************

      TOFF = OFFSET
      DO ICOL = 1, CHRCOL
         IF (ICOL .EQ. 2) THEN
            TOFF = 0
         END IF
C        DICTIONARY.

         DO I = 1, NNAMES(ICOL)
           LOC = DPOINT(I,ICOL,1)
           LEN = DPOINT(I,ICOL,2)

C ... Using malloc/free -- let system manage void space. Return
C     memory to system via 'free'.  The value given to memret
C     is a flag which tells the system that this is a 'safe' free
C     which should actually execute. (Major Kludge...)
           if (len .gt. 0) then
             LASTER = SUCCESS
             memret = -999
             oldadr = loc+myloc-1
             call exmemy(-len, oldadr, memret)
             if (memret .lt. 0 .or. memret .gt. len) then
               write (*,*) 'ERROR in MXDEL', memret, len
               laster = ilblk
             end if
           end if
         end do
      end do

      RETURN
      END
