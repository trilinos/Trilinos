C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details
      SUBROUTINE MYLONG (NAME1, NEWLEN, NEWLOC, MYV, MYCHAR, MYLOC,
     *   MYCLOC, UCLOC, COFFST, OFFSET,
     *   DICT, DPOINT, LDICT, NNAMES, VOID, LVOID, NVOIDS,
     *   FILL, FDATA, CFILL, CFDATA, CHRNUM, CHRCOL, MAXSIZ, LASTER)

      IMPLICIT INTEGER (A-Z)
      INCLUDE 'params.inc'

C***********************************************************************

C     NAME1    Name of the vector which changes length
               CHARACTER*8 NAME1
C     NEWLEN   The new length of the vector (character units)
C     NEWLOC   The new location of the vector (returned)
C     MYV      Internal reference array
               DIMENSION MYV(*)
C     MYCHAR   Internal reference array.
               CHARACTER MYCHAR(*)
C     MYLOC    Address of internal array
C     MYCLOC   Address of internal character array.
C     UCLOC    Address of user's character array.
C     COFFST   Offset between internal numeric array and user's
C              character array.
C     OFFSET   Address offset from internal array to user's array
C     DICT     Dictionary name table
C     DPOINT   Dictionary pointer table
C     LDICT    Dimension of dictionary
C     NNAMES   Number of names in the dictionary
               CHARACTER DICT(LDICT,CHRCOL)
               DIMENSION DPOINT(LDICT,CHRCOL,3), NNAMES(2)
C     VOID     Void table
C     LVOID    Dimension of void table
C     NVOIDS   Number of voids
               DIMENSION VOID(LVOID,CHRCOL,2), NVOIDS(2)
C     FILL     Flag for data fill.
C     FDATA    Data for fill.
               LOGICAL FILL
C     CFILL    Flag for character data fill.
C     CFDATA   Data for fill.
               LOGICAL CFILL
               CHARACTER*1 CFDATA
C     CHRNUM   Number of characters per numeric storage unit
C     CHRCOL   Number of column for character names.
C     MAXSIZ   Dimension of static character array.
C     LASTER   Error return

C***********************************************************************

      INTLEN = (NEWLEN + CHRNUM - 1) / CHRNUM

C     Get current location and length.

      CALL MYFIND (NAME1, DICT, DPOINT, LDICT, NNAMES,
     *   CHRCOL, LASTER, ROW)
      IF (LASTER .NE. SUCESS) RETURN

C     Save the current location of the array.

      OLDLOC = DPOINT(ROW,CHRCOL,1)
      OLDLEN = DPOINT(ROW,CHRCOL,2)
      OLDCL = DPOINT(ROW,CHRCOL,3)
      NEWLOC = (OLDLOC - 1) * CHRNUM + 1 + COFFST

      memret = -998
      memlen = intlen

C ... If the old length == 0, then we don't have a valid pointer.
C      Need to call malloc instead of realloc.
      if (oldlen .eq. 0) then
         memret = 0
        call exmemy(memlen, newadr, memret)
      else

C ... Passing a size of 0 to realloc (via exmemy) is the
C     same as a 'free' which invalidates the pointer.
C     Supes wants the pointer to stay valid, so instead
C     we request a space of '1' to maintain a valid pointer.
        if (intlen .eq. 0) then
          memlen = 1
        end if

        if (chrcol .eq. 1) then
          newadr = oldloc+myloc-1
        else
          newadr = oldloc+mycloc-1
        end if

        call exmemy(-memlen, newadr, memret)
      end if

      if (memret .lt. 0 .or. memret .gt. memlen) then
        write (*,*) 'ERROR in mylong ', memret, memlen
         laster = ilblk
         return
      end if

      IF (LASTER .NE. SUCESS) RETURN

      if (chrcol .eq. 1) then
        DPOINT(ROW,CHRCOL,1) = newadr+1-myloc
      else
        DPOINT(ROW,CHRCOL,1) = newadr+1-mycloc
      end if
      NEWLOC = (DPOINT(ROW,CHRCOL,1) - 1) * CHRNUM + 1 + COFFST
      DPOINT(ROW,CHRCOL,2) = INTLEN
      DPOINT(ROW,CHRCOL,3) = NEWLEN

C     Perform data fill if appropriate.

      IF (CFILL) THEN
         I1 = NEWLOC + UCLOC - MYCLOC + OLDCL
         I2 = I1 + NEWLEN - OLDCL - 1
         DO 130 I = I1, I2
            MYCHAR(I) = CFDATA
  130    CONTINUE
      END IF

      RETURN
      END
