C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details
      SUBROUTINE MXLONG (NAME1, NEWLEN, NEWLOC, MYV, MYCHAR, MYLOC,
     *   MYCLOC, UCLOC, COFFST, OFFSET,
     *   DICT, DPOINT, LDICT, NNAMES, VOID, LVOID, NVOIDS,
     *   FILL, FDATA, CFILL, CFDATA, CHRNUM, CHRCOL, LASTER)

      IMPLICIT INTEGER (A-Z)
      INCLUDE 'params.inc'

C***********************************************************************

C     NAME1    Name of the vector which changes length
               CHARACTER*8 NAME1
C     NEWLEN   The new length of the vector
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
C     LASTER   Error return

C***********************************************************************

C     Get current location and length.

      CALL MXFIND (NAME1, DICT, DPOINT, LDICT, NNAMES,
     *   CHRCOL, LASTER, ROW)
      IF (LASTER .NE. SUCESS) RETURN

C     Save the current location of the array.

      OLDLOC = DPOINT(ROW,1,1)
      OLDLEN = DPOINT(ROW,1,2)

      LASTER = SUCESS
      oldadr = oldloc+myloc-1
      memret = -998
      memlen = newlen

C ... If the old length == 0, then we don't have a valid pointer.
C      Need to call malloc instead of realloc.
      if (oldlen .eq. 0) then
        memret = 0
        if (memlen .gt. 0) then
           call exmemy(memlen, oldadr, memret)
        end if
      else
        call exmemy(-memlen, oldadr, memret)
      end if

      if (memret .lt. 0 .or. memret .gt. memlen) then
        laster = ilblk
        write (*,*) 'ERROR in mxlong ', memret, memlen
        return
      end if

      IF (LASTER .NE. SUCESS) RETURN

      DPOINT(ROW,1,1) = oldadr+1-myloc
      NEWLOC = DPOINT(ROW,1,1) + OFFSET
      DPOINT(ROW,1,2) = NEWLEN

C     Perform data fill if appropriate.

      IF (FILL) THEN
         DO 120 I = DPOINT(ROW,1,1)+OLDLEN, DPOINT(ROW,1,1)+NEWLEN-1-7,8
            MYV(I+0) = FDATA
            MYV(I+1) = FDATA
            MYV(I+2) = FDATA
            MYV(I+3) = FDATA
            MYV(I+4) = FDATA
            MYV(I+5) = FDATA
            MYV(I+6) = FDATA
            MYV(I+7) = FDATA
  120    CONTINUE
         DO 130 J = I, DPOINT(ROW,1,1)+NEWLEN-1
            MYV(J) = FDATA
  130    CONTINUE
      END IF

      RETURN
      END
