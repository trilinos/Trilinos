C    Copyright(C) 2014-2017 National Technology & Engineering Solutions of
C    Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
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

C $Id: inrenm.f,v 1.2 1992/01/10 23:14:40 gdsjaar Exp $
C $Log: inrenm.f,v $
C Revision 1.2  1992/01/10 23:14:40  gdsjaar
C Fixed problem with renumbering input
C
c Revision 1.1.1.1  1990/11/30  11:10:10  gdsjaar
c FASTQ Version 2.0X
c
c Revision 1.1  90/11/30  11:10:09  gdsjaar
c Initial revision
c 
C
CC* FILE: [.MAIN]INRENM.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE INRENM (MSC, N23, CFLAG, RIN, IIN, IFOUND, NUMBER,
     &   NOROOM)
C***********************************************************************
C
C  SUBROUTINE INRENM = INPUTS A RENUMBERING CARD
C
C***********************************************************************
C
      DIMENSION NUMBER (MSC), RIN (IFOUND), IIN (IFOUND)
C
      CHARACTER * 80 NUMBER, CFLAG * 72
C
      LOGICAL NOROOM
C
      NOROOM = .TRUE.
C
      N23 = N23 + 1
      IF (N23 .GT. MSC)RETURN
      NUMBER (N23) = ' '
      NUMBER (N23) (1:5) = CFLAG (1:5)
C
C  INPUT A POINT - LINE - POINT CARD
C
      IF (CFLAG (1:5) .EQ. 'P-L-P') THEN
         IFOUND = MIN0 (IFOUND, 15)
         DO 100 IJ = 1, IFOUND
            I2 =  (IJ + 1) * 5
            I1 = I2 - 4
            WRITE (NUMBER (N23) (I1:I2), 10000)IIN (IJ)
  100    CONTINUE
C
C  INPUT AN X, Y LOCATION RENUMBERING CARD
C
      ELSEIF (CFLAG (1:3) .EQ. 'X-Y') THEN
         WRITE (NUMBER (N23) (11:20), 10010) RIN (1)
         WRITE (NUMBER (N23) (21:30), 10010) RIN (2)
C
C  INPUT A NODE UNIQUE ID RENUMBERING CARD
C
      ELSEIF (CFLAG (1:4) .EQ. 'NODE') THEN
         IFOUND = MIN0 (IFOUND, 7)
         DO 110 IJ = 1, IFOUND
            I2 =  ( (IJ + 1) * 10)
            I1 = I2 - 9
            WRITE (NUMBER (N23) (I1:I2), 10020)IIN (IJ)
  110    CONTINUE
C
C  INDICATE ERROR IN RENUMBERING FLAG
C
      ELSE
         N23 = N23 - 1
         WRITE ( * , 10030) CFLAG (1:5)
      ENDIF
C
      NOROOM = .FALSE.
      RETURN
C
10000 FORMAT (I5)
10010 FORMAT (1PE10.3)
10020 FORMAT (I10)
10030 FORMAT (' RENUMBERING KEY WORD: ', A5, ' IS NOT ALLOWABLE',  / ,
     &   ' THIS RENUMBERING LIST WILL NOT BE INPUT INTO DATABASE')
      END
