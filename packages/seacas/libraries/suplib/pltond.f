C Copyright(C) 2009 Sandia Corporation. Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
C certain rights in this software.
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


      SUBROUTINE PLTOND (MD, TITLE, NTITLE, XX, YY, XLAB, YLAB, NP,     
     1   CURVE, TPLOT)
C
C     DOUBLE PRECISION VERSION
C        OUTPUT DATA IN GRAFAID NEUTRAL FILE FORMAT FOR ONE PLOT
C
C ----------------------------------------------------------------------
C VARIABLES:
C     Name     Type     Description
C     ----     ----     ----------------------------------
C     MD       integer  Neutral file unit number - file already opened
C     NTITLE   integer  Number of lines in TITLE
C     XLAB     char*40  Label for X-Axis
C     YLAB     char*40  Label for Y-Axis
C     NP       integer  Number of Points
C     CURVE    char*16  Curve Name
C     TPLOT    logical  .TRUE.  if X-variable is time
C                       .FALSE. if X-variable is not time
C
C ARRAYS:
C     Name     Dimension   Type     Description
C     ----     ---------   ----     ----------------------------------
C     TITLE    NTITLE      char*80  Title for current curve
C     XX       NP          real     X-variable data
C     YY       NP          real     Y-variable data
C
C     SUBROUTINES AND FUNCTIONS CALLED:
C
      CHARACTER*(*) TITLE(*)
      CHARACTER*(*) XLAB,YLAB
      REAL*8 XX(1), YY(1), XMN, XMX, YMN, YMX
C
      LOGICAL MONO, TPLOT
      CHARACTER*(*) CURVE
      CHARACTER*11 BEGIN
      CHARACTER*9 ECURVE
      CHARACTER*1 COMMA,AUX
      CHARACTER*4 XTYP, AXTYP
      DATA COMMA/','/, AXTYP/'NOLO'/
      DATA BEGIN/'BEGIN CURVE'/, ECURVE/'END CURVE'/, AUX/'F'/
C
C     ...LOCATE MINIMUM AND MAXIMUM VALUES AND CHECK FOR NONMONOTONIC DATA
C
      MONO = .TRUE.
      IF (TPLOT) THEN
         XMN = XX(1)
         XMX = XX(NP)
       ELSE
         XMN = XX(1)
         XMX = XX(1)
         DO 10 I=2,NP
            XMN = MIN(XMN, XX(I))
            IF (XMX .GE. XX(I)) THEN
               MONO = .FALSE.
             ELSE
               XMX = XX(I)
             END IF
   10       CONTINUE
       END IF
      IF (MONO) THEN
         XTYP = 'MONO'
       ELSE
         XTYP = 'NONM'
       END IF
C
      YMN=YY(1)
      YMX=YY(1)
      DO 20 I=2,NP
         YMN = MIN(YMN, YY(I))
         YMX = MAX(YMX, YY(I))
   20    CONTINUE
C
C     BEGIN TO WRITE CURVE PACKET
C
      WRITE (MD, 40) BEGIN,COMMA,CURVE
      WRITE (MD, 50) NTITLE,COMMA,TITLE(1)
      DO 30 I=2,NTITLE
         WRITE (MD, 60) TITLE(I)
   30    CONTINUE
      WRITE (MD, 70) XLAB
      WRITE (MD, 70) YLAB
      WRITE (MD, 80) XMN,COMMA,XMX,COMMA,YMN,COMMA,YMX,COMMA,           
     1               NP,COMMA,AUX
      WRITE (MD, 90) AXTYP,COMMA,XTYP,COMMA
C
C     WRITE DATA PAIRS
C
      DO 35 III=1,NP
         WRITE (MD, 100) XX(III),COMMA,YY(III)
   35  CONTINUE
C
C     WRITE END OF CURVE PACKET
C
      WRITE (MD, 110) ECURVE,COMMA,CURVE
      RETURN
C
   40 FORMAT (A11,A1,A16)
   50 FORMAT (I1,A1,A)
   60 FORMAT (A)
   70 FORMAT (A)
   80 FORMAT (1PE15.7,A1,1PE15.7,A1,1PE15.7,A1,1PE15.7,A1,I5,A1,A1)
   90 FORMAT (A4,A1,A4,A1)
  100 FORMAT (1PE15.7,A1,1PE15.7)
  110 FORMAT (A9,A1,A16)
      END
