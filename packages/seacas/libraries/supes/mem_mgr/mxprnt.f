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
      SUBROUTINE MXPRNT (NAME1, UNIT, NAME2, MYV, RMYV, OFFSET,
     *   DICT, DPOINT, LDICT, NNAMES, CHRCOL, NCOLP, LASTER)
C
      IMPLICIT INTEGER (A-Z)
      INCLUDE 'params.inc'
C
C***********************************************************************
C
C     NAME1    Name of array to be printed
               CHARACTER*8 NAME1
C     UNIT     Output unit number.
C     NAME2    Type of array to be printed
               CHARACTER*(*) NAME2
C     MYV      Internal integer array
               INTEGER MYV(*)
C     RMYV     Internal real array
               REAL RMYV(*)
C     DICT     Dictionary name table
C     DPOINT   Dictionary pointer table
C     LDICT    Dimension of dictionary
C     NNAMES   Number of names in the dictionary
               CHARACTER*8 DICT(LDICT,CHRCOL)
               DIMENSION DPOINT(LDICT,CHRCOL,2), NNAMES(2)
C     CHRCOL   Number of column for character names.
C     NCOLP    Number of print columns
C     LASTER   Error return
C
C***********************************************************************
C
C     FIND NAME1 IN DICTIONARY
C
      CALL MXFIND (NAME1, DICT, DPOINT, LDICT, NNAMES,
     *   CHRCOL, LASTER, ROW)
      IF (LASTER .NE. SUCESS) RETURN
      DELTA = DPOINT(ROW,1,1) - 1
      WRITE(UNIT,10000) DICT(ROW,1), DPOINT(ROW,1,1)+OFFSET,
     *   DPOINT(ROW,1,2)
      IF (DPOINT(ROW,1,2) .LT. 0) THEN
         LASTER = DEFRON
         WRITE (UNIT, *) 'THIS VECTOR WAS RESERVED IN THE DEFERRED '//
     *      'MODE AND IS NOT YET RESOLVED.'
         RETURN
      END IF
      IF (NAME2(1:1) .EQ. 'R') THEN
C
C        VECTOR IS REAL
C
         NCOL=(NCOLP-11)/13
         NROW=DPOINT(ROW,1,2)/NCOL+1
         NGRP=NROW/10+1
         DO 110 IGRP=1,NGRP
            WRITE(UNIT,10010)
            NPRT=(IGRP-1)*10*NCOL
            NREM=DPOINT(ROW,1,2)-NPRT
            NROW=(NREM+NCOL-1)/NCOL
            NROW=MIN0(10,NROW)
            DO 100 IROW=1,NROW
               J=NPRT+1+(IROW-1)*NCOL
               KU=MIN0(DPOINT(ROW,1,2),J+NCOL-1)
               WRITE(UNIT,10020)J,(RMYV(K),K=J+DELTA,KU+DELTA)
  100       CONTINUE
  110    CONTINUE
C
      ELSE IF (NAME2(1:1) .EQ. 'I') THEN
C
C        VECTOR IS INTEGER
C
         NCOL=(NCOLP-11)/10
         NROW=DPOINT(ROW,1,2)/NCOL+1
         NGRP=NROW/10+1
         DO 130 IGRP=1,NGRP
            WRITE(UNIT,10010)
            NPRT=(IGRP-1)*10*NCOL
            NREM=DPOINT(ROW,1,2)-NPRT
            NROW=(NREM+NCOL-1)/NCOL
            NROW=MIN0(10,NROW)
            DO 120 IROW=1,NROW
               J=NPRT+1+(IROW-1)*NCOL
               KU=MIN0(DPOINT(ROW,1,2),J+NCOL-1)
               WRITE(UNIT,10030)J,(MYV(K),K=J+DELTA,KU+DELTA)
  120       CONTINUE
  130    CONTINUE
C
      ELSE
C
C        TYPE IS UNKNOWN
C
         LASTER = BDTYPE
         RETURN
C
      END IF
      LASTER = SUCESS
      RETURN
10000 FORMAT(//' ARRAY NAME = ',A,3X,'LOCATION = ',I12,3X,
     *  'LENGTH = ',I8)
10010 FORMAT(' ')
10020 FORMAT(1X,I6,':',9(2X,1PE11.4))
10030 FORMAT(1X,I6,':',12(2X,I8))
      END
