C    Copyright (c) 2014, Sandia Corporation.
C    Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
C    the U.S. Governement retains certain rights in this software.
C    
C    Redistribution and use in source and binary forms, with or without
C    modification, are permitted provided that the following conditions are
C    met:
C    
C        * Redistributions of source code must retain the above copyright
C          notice, this list of conditions and the following disclaimer.
C    
C        * Redistributions in binary form must reproduce the above
C          copyright notice, this list of conditions and the following
C          disclaimer in the documentation and/or other materials provided
C          with the distribution.
C    
C        * Neither the name of Sandia Corporation nor the names of its
C          contributors may be used to endorse or promote products derived
C          from this software without specific prior written permission.
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

C $Id: addsnp.f,v 1.1 1990/11/30 11:03:07 gdsjaar Exp $
C $Log: addsnp.f,v $
C Revision 1.1  1990/11/30 11:03:07  gdsjaar
C Initial revision
C
C
CC* FILE: [.MAIN]ADDSNP.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE ADDSNP(MSNAP,SNAPDX,NSNAP,INDEX,VALUE,ERR)
C***********************************************************************
C
C   SUBROUTINE ADDSNP = ADDS SNAP GRID LINE DEFINITIONS
C
C***********************************************************************
C
C  VARIABLES USED:
C     MSNAP  = DIMENSION OV SNAP ARRAYS
C     SNAPDX = THE SNAP GRID VALUES ARRAY (X AND Y)
C     NSNAP  = THE NUMBER OF SNAP GRID VALUES IN X AND Y
C     INDEX  = 1 FOR X VALUES, 2 FOR Y VALUES
C     VALUE  = THE GRID VALUE TO BE ADDED
C     KOUNT  = THE LOCATION OF THE SNAPDX VALUE JUST LESS THAN VALUE
C
C***********************************************************************
C
      DIMENSION SNAPDX(2,MSNAP),NSNAP(2)
C
      LOGICAL ERR
C
C  ADD THE SNAP GRID VALUE WHERE IT FITS IN NUMERICAL ORDER
C
      ERR=.FALSE.
C
      IF(NSNAP(INDEX).GT.0)THEN
         KOUNT=0
         DO 100 I=1,NSNAP(INDEX)
            IF(VALUE.LT.SNAPDX(INDEX,I))GOTO 110
            KOUNT=I
C
C  IF THE VALUE IS ALREADY THERE, THEN DON'T ADD IT AGAIN - JUST RETURN
C
            IF(VALUE.EQ.SNAPDX(INDEX,I))RETURN
C
  100    CONTINUE
  110    CONTINUE
         IF(NSNAP(INDEX).EQ.MSNAP)THEN
            CALL MESAGE('** NO MORE ROOM FOR ADDITIONAL GRID LINES **')
            WRITE(*,10000)MSNAP
            ERR=.TRUE.
            RETURN
         ENDIF
         NSNAP(INDEX)=NSNAP(INDEX)+1
         DO 120 I=NSNAP(INDEX),KOUNT+2,-1
            SNAPDX(INDEX,I)=SNAPDX(INDEX,I-1)
  120    CONTINUE
         SNAPDX(INDEX,KOUNT+1)=VALUE
C
C  JUST PUT THE FIRST VALUE WHERE IT BELONGS
C
      ELSE
         NSNAP(INDEX)=1
         SNAPDX(INDEX,1)=VALUE
      ENDIF
C
      RETURN
C
10000 FORMAT(' THE MAXIMUM NUMBER OF GRID LINES IS: ',I10)
C
      END
