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
      SUBROUTINE MYFIND (NAME1, DICT, DPOINT, LDICT, NNAMES,
     *   CHRCOL, LASTER, ROW)
C
      IMPLICIT INTEGER (A-Z)
      INCLUDE 'params.inc'
C
C***********************************************************************
C
C     NAME1    Name to be found
               CHARACTER*8 NAME1
C     DICT     Dictionary name table
C     DPOINT   Dictionary pointer table
C     LDICT
               CHARACTER*8 DICT(LDICT,CHRCOL)
               DIMENSION DPOINT(LDICT,CHRCOL,3)
C     NNAMES   Number of names in the dictionary
               DIMENSION NNAMES(2)
C     CHRCOL   Column number for character array names.
C     LASTER   Error return
C     ROW      Location of found name or place to insert new name
C
C***********************************************************************
C
      CALL SRCHC (DICT, 1, NNAMES(1), NAME1, ERR, ROW)
      IF (ERR .EQ. 1) THEN
         IF (DPOINT(ROW,1,3) .NE. -1) THEN
C
C           The found name is a name for a character array.
            LASTER = SUCESS
C
         ELSE
C
C           The names was found and is of numeric type.
            LASTER = WRTYPE
         END IF
C
      ELSE IF (CHRCOL .EQ. 1) THEN
C
C        ENTRY NOT FOUND.
C
         LASTER = NONAME
C
      ELSE
         CALL SRCHC (DICT(1,2), 1, NNAMES(2), NAME1, ERR, ROW)
         IF (ERR .EQ. 1) THEN
            LASTER = SUCESS
         ELSE
            LASTER = NONAME
         END IF
      END IF
C
      RETURN
      END
