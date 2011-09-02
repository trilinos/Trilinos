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
      SUBROUTINE MYDEL (NAME1, DICT, DPOINT, LDICT, NNAMES, VOID,
     *   LVOID, NVOIDS, CHRCOL, LASTER, MYLOC, MYCLOC)
C
      IMPLICIT INTEGER (A-Z)
      INCLUDE 'params.inc'
C
C     This routine removes a name from the dictionary and returns the
C     available space to the void table.
C
C***********************************************************************
C
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
C
C***********************************************************************
C
C     FIND NAME1 IN DICTIONARY.
C
      CALL MYFIND (NAME1, DICT, DPOINT, LDICT, NNAMES,
     *   CHRCOL, LASTER, ROW)
      IF (LASTER .NE. SUCESS) RETURN
C
      LOC = DPOINT(ROW,CHRCOL,1)
      LEN = DPOINT(ROW,CHRCOL,2)
C
C     DELETE DICTIONARY ENTRY.
C
      CALL SHFTC (DICT(1,CHRCOL), CHRCOL*LDICT, ROW+1,
     *   NNAMES(CHRCOL), 1)
      CALL SHFTI (DPOINT(1,CHRCOL,1), CHRCOL*LDICT, 3, ROW+1,
     *   NNAMES(CHRCOL), 1)
      NNAMES(CHRCOL) = NNAMES(CHRCOL) - 1
      IF (LEN .LE. 0) RETURN
C
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
