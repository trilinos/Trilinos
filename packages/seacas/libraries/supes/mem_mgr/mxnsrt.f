C    Copyright(C) 2008-2017 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
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
C    * Neither the name of NTESS nor the names of its
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
      SUBROUTINE MXNSRT (NAME1, NEWLOC, NEWLEN,
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
C     NEWLEN   Length of storage
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
      IF (NNAMES(1) .GE. LDICT) THEN
         LASTER = DFULL
         RETURN
      END IF
C
C     FIND NAME1 IN DICTIONARY
C
      CALL MXFIND (NAME1, DICT, DPOINT, LDICT, NNAMES,
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
      CALL SHFTC (DICT, CHRCOL*LDICT, ROW, NNAMES(1), -1)
      CALL SHFTI (DPOINT, CHRCOL*LDICT, 3, ROW, NNAMES(1), -1)
      NNAMES(1) = NNAMES(1) + 1
      DICT(ROW,1) = NAME1
      DPOINT(ROW,1,1) = NEWLOC
      DPOINT(ROW,1,2) = NEWLEN
      DPOINT(ROW,1,3) = -1
      RETURN
      END
