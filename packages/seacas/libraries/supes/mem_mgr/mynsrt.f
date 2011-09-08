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
      SUBROUTINE MYNSRT (NAME1, NEWLOC, NUMLEN, CLEN,
     *   DICT, DPOINT, LDICT,
     *   NNAMES, CHRCOL, LASTER)
C
      IMPLICIT INTEGER (A-Z)
      INCLUDE 'params.inc'
C
C     This routine updates the dictionary with a new name (if it is new)
C     and updates the location and length tables.  The length of the
C     dictionary is checked before the new name is added.  If LASTER is
C     not returned with a value of SUCESS, the tables are unchanged.
C
C***********************************************************************
C
C     NAME1    Name to be inserted
               CHARACTER*8 NAME1
C     NEWLOC   Location of storage
C     NUMLEN   Numeric length of storage
C     CLEN     Character length of storage
C     DICT     Dictionary name table
C     DPOINT   Dictionary pointer table
C     LDICT    Dimensioned size of dictionary
C     NNAMES   Number of entries in the dictionary
               CHARACTER*8 DICT(LDICT,CHRCOL)
               DIMENSION DPOINT(LDICT,CHRCOL,3), NNAMES(2)
C     CHRCOL   Number of column for character names.
C     LASTER   Error return
C
C***********************************************************************
C
C     IS THERE ROOM IN THE DICTIONARY?
C
      IF (NNAMES(CHRCOL) .GE. LDICT) THEN
         LASTER = DFULL
         RETURN
      END IF
C
C     FIND NAME1 IN DICTIONARY
C
      CALL MYFIND (NAME1, DICT, DPOINT, LDICT, NNAMES,
     *   CHRCOL, LASTER, ROW)
      IF (LASTER .EQ. WRTYPE) THEN
         RETURN
      ELSE IF (LASTER .EQ. SUCESS) THEN
         LASTER = BDNAME
         RETURN
      ELSE IF (LASTER .EQ. NONAME) THEN
         LASTER = SUCESS
      END IF
C
C     UPDATE DICTIONARY.
C
      CALL SHFTC (DICT(1,CHRCOL), CHRCOL*LDICT, ROW, NNAMES(CHRCOL), -1)
      CALL SHFTI (DPOINT(1,CHRCOL,1), CHRCOL*LDICT, 3, ROW,
     *   NNAMES(CHRCOL), -1)
      NNAMES(CHRCOL) = NNAMES(CHRCOL) + 1
      DICT(ROW,CHRCOL) = NAME1
      DPOINT(ROW,CHRCOL,1) = NEWLOC
      DPOINT(ROW,CHRCOL,2) = NUMLEN
      DPOINT(ROW,CHRCOL,3) = CLEN
      RETURN
      END
