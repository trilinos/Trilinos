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
      SUBROUTINE MYPRNT (NAME1, UNIT, MYCV, OFFSET, TOFFST,
     *   DICT, DPOINT, LDICT, NNAMES, CHRNUM,
     *   CHRCOL, NCOLP, WRDSIZ, LASTER)
C
      IMPLICIT INTEGER (A-Z)
      INCLUDE 'params.inc'
C
C***********************************************************************
C
C     NAME1    Name of array to be printed
               CHARACTER*8 NAME1
C     UNIT     Output unit number.
C     MYCV     Internal character array
               CHARACTER*1 MYCV(*)
C     OFFSET   Offset between internal reference and users reference
C              string.
C     TOFFST   Offset between internal reference and internal character
C              array.
C     DICT     Dictionary name table
C     DPOINT   Dictionary pointer table
C     LDICT    Dimension of dictionary
C     NNAMES   Number of names in the dictionary
               CHARACTER*8 DICT(LDICT,CHRCOL)
               DIMENSION DPOINT(LDICT,CHRCOL,3), NNAMES(2)
C     CHRCOL   Number of column for character names.
C     NCOLP    Number of print columns
C     WRDSIZ   Number of characters to group together in printing.
C     LASTER   Error return
C
C***********************************************************************
C
C
C     Check worklength
C
      IF (WRDSIZ .LT. 1 .OR. WRDSIZ+2+11 .GT. NCOLP) THEN
         LASTER = BADLEN
         RETURN
      END IF
C     FIND NAME1 IN DICTIONARY
C
      CALL MYFIND (NAME1, DICT, DPOINT, LDICT, NNAMES,
     *   CHRCOL, LASTER, ROW)
      IF (LASTER .NE. SUCESS) RETURN
      DELTA = (DPOINT(ROW,CHRCOL,1) - 1) * CHRNUM + 1 + OFFSET
      WRITE(UNIT,10000) DICT(ROW,CHRCOL),
     *   DELTA,
     *   DPOINT(ROW,CHRCOL,3)
      IF (DPOINT(ROW,CHRCOL,2) .LT. 0) THEN
         LASTER = DEFRON
         WRITE (UNIT, *) 'THIS VECTOR WAS RESERVED IN THE DEFERRED '//
     *      'MODE AND IS NOT YET RESOLVED.'
         RETURN
      END IF
C
      DELTA = (DPOINT(ROW,CHRCOL,1) - 1) * CHRNUM + TOFFST
      NCOL = (NCOLP - 11) / (WRDSIZ + 2)
      NROW = (DPOINT(ROW,CHRCOL,3) + WRDSIZ * NCOL - 1)
     *   / (WRDSIZ * NCOL)
      NGRP = (NROW + 9) / 10
      DO 110 IGRP = 1, NGRP
         WRITE(UNIT,10010)
         NPRT = (IGRP - 1) * 10 * NCOL * WRDSIZ
         NREM = DPOINT(ROW,CHRCOL,3) - NPRT
         NROW = (NREM + WRDSIZ * NCOL - 1) / (NCOL * WRDSIZ)
         NROW = MIN(10, NROW)
         J = NPRT + 1
         DO 100 IROW = 1, NROW
            WRITE (UNIT, 10020) J,
     *         ((MYCV(K), K=J+DELTA+WRDSIZ*(IWRD-1),
     *         MIN(DELTA+DPOINT(ROW,CHRCOL,3),
     *         J+DELTA+WRDSIZ*IWRD-1)),
     *         ' ', ' ', IWRD = 1, NCOL)
            J = J + NCOL * WRDSIZ
  100    CONTINUE
  110 CONTINUE
C
      LASTER = SUCESS
      RETURN
10000 FORMAT('0'/'0ARRAY NAME = ',A,3X,'LOCATION = ',I8,3X,
     *  'LENGTH = ',I8)
10010 FORMAT(' ')
10020 FORMAT(1X,I6,':',132A1)
      END
