C Copyright (c) 2007 Sandia Corporation. Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Governement
C retains certain rights in this software.
C 
C Redistribution and use in source and binary forms, with or without
C modification, are permitted provided that the following conditions are
C met:
C 
C     * Redistributions of source code must retain the above copyright
C       notice, this list of conditions and the following disclaimer.
C 
C     * Redistributions in binary form must reproduce the above
C       copyright notice, this list of conditions and the following
C       disclaimer in the documentation and/or other materials provided
C       with the distribution.  
C 
C     * Neither the name of Sandia Corporation nor the names of its
C       contributors may be used to endorse or promote products derived
C       from this software without specific prior written permission.
C 
C THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
C "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
C LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
C A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
C OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
C SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
C LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
C DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
C THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
C (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
C OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
C 

      SUBROUTINE GETBND(LBLK,NE,X,IND,NP,XMIN,XMAX,NDIM,ILO,IUP,
     *                  ISCR,RSCR )
C
C-----------------------------------------------------------------------
C     
C DESCRIPTION:
C
C  FIND THE ELEMENTS IN A SORTED ARRAY X WHOSE VALUES FALL IN THE
C  INTERVAL BETWEEN XMIN AND XMAX. NO ELEMENTS HAVING
C  VALUES EQUAL TO XMIN OR XMAX ARE INCLUDED. SINCE THE ARRAY IS
C  SORTED, THE ELEMENTS CAN BE SPECIFIED BY THE UPPER AND
C  LOWER ELEMENT NUMBERS IN THE RANGE.
C
C     |    X(ILO)   . . . . . . . . . . . . .   X(IUP)    |
C    XMIN                                               XMAX       X>
C
C  IT IS ASSUMED THAT THE ARRAY X HAS BEEN SORTED IN INCREASING ORDER,
C  BUT THE ELEMENTS HAVE NOT BEEN MOVED.
C  THE SORTED LIST IS DETERMINED BY THE ARRAY INDX,
C  WHICH POSITIONS THE ORIGINAL UNSORTED X ARRAY ELEMENTS
C  IN THE SORTED LIST. THUS, THE 5TH ELEMENT IN THE SORTED LIST IS
C    X(IND(5))
C
C-----------------------------------------------------------------------
C
C  INPUT:
C
C  X      -  array in unsorted order
C  IND    -  index array giving the element order in the sorted list
C  NP     -  the number of particles in the list
C  XMIN   -  the lower limit of the interval
C  XMAX   -  the upper limit of the interval
C  NDIM   -  the dimension of the arrays
C
C  OUTPUT:
C
C  ILO    -  the first element in the sorted list .gt. xmin
C  IUP    -  the last element in the sorted list .lt. xmax
C
C-----------------------------------------------------------------------
C
C
      include 'tapes.blk'
C
      DIMENSION
     *  X(NDIM),IND(NDIM),XMIN(LBLK),XMAX(LBLK),ILO(LBLK),IUP(LBLK)
      DIMENSION
     *  ISCR(5*LBLK),RSCR(LBLK)
C
C INTEGER SCRATCH SPACE
      ISP    = 1
      LIL    = ISP
      ISP    = ISP + LBLK
      LIU    = ISP
      ISP    = ISP + LBLK
      LIT    = ISP
      ISP    = ISP + LBLK
      LINDX1 = ISP
C      ISP    = ISP + LBLK
C      LINDX2 = ISP
C      ISP    = ISP + LBLK
C REAL SCRATCH SPACE
      ISP   = 1
      LXTST = ISP
C      ISP   = ISP + LBLK
C
      DO 200 J = 1, NE, LBLK
            N = MIN(LBLK,NE-J+1)
C  SEARCH TO FIND THE FIRST ELEMENT .GE. XMIN
      CALL SRCHGE(LBLK,N,X,IND,XMIN(J),1,     NP,NDIM,ILO(J),
     *     ISCR(LIL),ISCR(LIU),ISCR(LIT),ISCR(LINDX1),ISCR(LINDX1),
     *     RSCR(LXTST) )
C  SEARCH TO FIND THE FIRST ELEMENT .GT. XMAX
      CALL SRCHGT(LBLK,N,X,IND,XMAX(J),ILO(J),NP,NDIM,IUP(J),
     *     ISCR(LIL),ISCR(LIU),ISCR(LIT),ISCR(LINDX1),ISCR(LINDX1),
     *     RSCR(LXTST) )
C  THE PREVIOUS ELEMENT IS THE LAST ONE .LT. XMAX
      DO 100 JJ = 1, N
       IUP(J+JJ-1)=IUP(J+JJ-1) - 1
 100  CONTINUE
C
 200  CONTINUE
      RETURN
      END
C
