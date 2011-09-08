C    Copyright(C) 2008 Sandia Corporation.  Under the terms of Contract
C    DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
C    certain rights in this software
C    
C    Redistribution and use in source and binary forms, with or without
C    modification, are permitted provided that the following conditions are
C    met:
C    
C    * Redistributions of source code must retain the above copyright
C       notice, this list of conditions and the following disclaimer.
C              
C    * Redistributions in binary form must reproduce the above
C      copyright notice, this list of conditions and the following
C      disclaimer in the documentation and/or other materials provided
C      with the distribution.
C                            
C    * Neither the name of Sandia Corporation nor the names of its
C      contributors may be used to endorse or promote products derived
C      from this software without specific prior written permission.
C                                                    
C    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
C    "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
C    LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
C    A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
C    OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
C    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
C    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
C    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
C    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
C    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
C    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
C    
      SUBROUTINE MYLONG (NAME1, NEWLEN, NEWLOC, MYV, MYCHAR, MYLOC,
     *   MYCLOC, UCLOC, COFFST, OFFSET,
     *   DICT, DPOINT, LDICT, NNAMES, VOID, LVOID, NVOIDS,
     *   FILL, FDATA, CFILL, CFDATA, CHRNUM, CHRCOL, MAXSIZ, LASTER)
C
      IMPLICIT INTEGER (A-Z)
      INCLUDE 'params.inc'
C
C***********************************************************************
C
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
C
C***********************************************************************
C
      INTLEN = (NEWLEN + CHRNUM - 1) / CHRNUM
C
C     Get current location and length.
C
      CALL MYFIND (NAME1, DICT, DPOINT, LDICT, NNAMES,
     *   CHRCOL, LASTER, ROW)
      IF (LASTER .NE. SUCESS) RETURN
C
C     Save the current location of the array.
C
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
C
C     Perform data fill if appropriate.
C
      IF (CFILL) THEN
         I1 = NEWLOC + UCLOC - MYCLOC + OLDCL
         I2 = I1 + NEWLEN - OLDCL - 1
         DO 130 I = I1, I2
            MYCHAR(I) = CFDATA
  130    CONTINUE
      END IF
C
      RETURN
      END
