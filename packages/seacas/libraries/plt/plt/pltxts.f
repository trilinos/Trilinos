C Copyright (C) 2009 Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
C certain rights in this software
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

C $Id: pltxts.f,v 1.1 1993/07/16 16:50:02 gdsjaar Exp $ 
C $Log: pltxts.f,v $
C Revision 1.1  1993/07/16 16:50:02  gdsjaar
C Changed plt to library rather than single source file.
C 
C=======================================================================
      SUBROUTINE PLTXTS(X,Y,TEXT)
      REAL DEVCAP(23)
      REAL DEFOUT(7)
      COMMON /STATUS/DEVCAP,DEFOUT
      REAL DEVP(5)
      COMMON /DEVICE/DEVP
      REAL COLP(3)
      REAL PALETT(3,16)
      COMMON /COLOR/COLP,PALETT
      REAL TEXTP(40)
      COMMON /TEXT/TEXTP
      REAL VECTP(5)
      REAL XCUR
      REAL YCUR
      COMMON /VECTRC/VECTP,XCUR,YCUR
      INTEGER IDEX(200,2)
      INTEGER NVECT(200,2)
      REAL XSIZE(200,2)
      REAL YSIZE(200,2)
      REAL X0(2300,2)
      REAL Y0(2300,2)
      REAL X1(2300,2)
      REAL Y1(2300,2)
      COMMON /FONT/IDEX,NVECT,XSIZE,YSIZE,X0,Y0,X1,Y1
      REAL GRAPHP(100)
      COMMON /GRAPH/GRAPHP
      COMMON /MAPPAR/MAPP(11)
      REAL MAPP
      COMMON /STORAG/MEMORY(1000)
      LOGICAL CPUIFC
      CHARACTER*(*) TEXT
      CHARACTER*1 LASCHR,ESCCHR
      CHARACTER*20 ESC
      INTEGER ASCII
      LOGICAL STATUS,CHRCI
      DATA ESCCHR/'\\'/

      IFONT = 1
      CALL PLTSVV
      CALL PLTSTV(1,1.)
      CALL PLTSTV(2,TEXTP(37))
      CALL PLTRIM(TEXT,NCHAR)
      TEXTP(18) = X
      TEXTP(19) = Y
      IFLAG = 0
      YSAVE = Y
      YCHRSZ = TEXTP(1)
      LASCHR = 'M'
      I = 1
 2020 IF (.NOT. (I.LE.NCHAR)) GO TO 2040
      ASCII = ICHAR(TEXT(I:I))
      IF (ASCII.LT.1 .OR. ASCII.GT.126) THEN
         CALL CHRIC(ASCII,LASCHR,LI)
         CALL PLTFLU
         CALL SIORPT('PLTXTS','Invalid character "'//LASCHR(1:LI)//
     *               '" in text string; rest of string ignored',2)
         RETURN

      END IF

      IF (ASCII.EQ.ICHAR(ESCCHR) .AND. TEXT(I+1:I+1).EQ.ESCCHR) THEN
         I = I + 1

      ELSE IF (ASCII.EQ.ICHAR(ESCCHR)) THEN
         CALL PLTESC(TEXT,I,ESC)
         CALL CHRUP(ESC,ESC)
         IF (ESC.EQ.'^') THEN
            IF (IFLAG.NE.1) THEN
               IF (IFLAG.EQ.-1) THEN
                  CALL PLTNOR(YBUMP,YCHRSZ)
               END IF

               CALL PLTSUP(YBUMP,YCHRSZ)
               IFLAG = 1
               GO TO 2030

            END IF

         ELSE IF (ESC.EQ.'_') THEN
            IF (IFLAG.NE.-1) THEN
               IF (IFLAG.EQ.1) THEN
                  CALL PLTNOR(YBUMP,YCHRSZ)
               END IF

               CALL PLTSUB(YBUMP,YCHRSZ)
            END IF

            IFLAG = -1
            GO TO 2030

         ELSE IF (ESC.EQ.'-') THEN
            IF (IFLAG.NE.0) THEN
               CALL PLTNOR(YBUMP,YCHRSZ)
            END IF

            IFLAG = 0
            GO TO 2030

         ELSE IF (ESC.EQ.'CLO') THEN
            ASCII = 4

         ELSE IF (ESC.EQ.'CSQ') THEN
            ASCII = 5

         ELSE IF (ESC.EQ.'CDI') THEN
            ASCII = 6

         ELSE IF (ESC.EQ.'CCS') THEN
            ASCII = 7

         ELSE IF (ESC.EQ.'CX') THEN
            ASCII = 8

         ELSE IF (ESC.EQ.'CTR') THEN
            ASCII = 9

         ELSE IF (ESC.EQ.'CCI') THEN
            ASCII = 10

         ELSE IF (ESC.EQ.'CDO') THEN
            ASCII = 11

         ELSE IF (ESC.EQ.'LO') THEN
            ASCII = 12

         ELSE IF (ESC.EQ.'SQ') THEN
            ASCII = 13

         ELSE IF (ESC.EQ.'DI') THEN
            ASCII = 14

         ELSE IF (ESC.EQ.'CS') THEN
            ASCII = 15

         ELSE IF (ESC.EQ.'X') THEN
            ASCII = 16

         ELSE IF (ESC.EQ.'TR') THEN
            ASCII = 17

         ELSE IF (ESC.EQ.'CI') THEN
            ASCII = 18

         ELSE IF (ESC.EQ.'DO') THEN
            ASCII = 19

         ELSE IF (ESC.EQ.'PLUSMIN') THEN
            ASCII = 20

         ELSE IF (ESC.EQ.'LEQ') THEN
            ASCII = 21

         ELSE IF (ESC.EQ.'GEQ') THEN
            ASCII = 22

         ELSE IF (ESC.EQ.'NEQ') THEN
            ASCII = 23

         ELSE IF (ESC.EQ.'PRIME') THEN
            ASCII = 24

         ELSE IF (ESC.EQ.'NLEQ') THEN
            ASCII = 25

         ELSE IF (ESC.EQ.'NGEQ') THEN
            ASCII = 26

         ELSE IF (ESC.EQ.'LL') THEN
            ASCII = 27

         ELSE IF (ESC.EQ.'GG') THEN
            ASCII = 28

         ELSE IF (ESC.EQ.'SUM') THEN
            ASCII = 29

         ELSE IF (ESC.EQ.'NLT') THEN
            ASCII = 30

         ELSE IF (ESC.EQ.'NGT') THEN
            ASCII = 31

         ELSE IF (ESC.EQ.'APPROX') THEN
            ASCII = 127

         ELSE IF (ESC.EQ.'CR') THEN
            TEXTP(18) = X + TEXTP(30)*TEXTP(1)*TEXTP(29)
            GO TO 2030

         ELSE IF (ESC.EQ.'LF') THEN
            TEXTP(19) = YSAVE - TEXTP(30)*TEXTP(1)*TEXTP(28)
            YSAVE = TEXTP(19)
            GO TO 2030

         ELSE IF (ESC.EQ.'CL') THEN
            TEXTP(18) = X + TEXTP(30)*TEXTP(1)*TEXTP(29)
            TEXTP(19) = YSAVE - TEXTP(30)*TEXTP(1)*TEXTP(28)
            YSAVE = TEXTP(19)
            GO TO 2030

         ELSE IF (ESC.EQ.'ENG') THEN
            IFONT = 1
            GO TO 2030

         ELSE IF (ESC.EQ.'GR') THEN
            IFONT = 2
            GO TO 2030

         ELSE IF (ESC.EQ.'DDLINE') THEN
            CALL PLTSVV
            CALL PLTSTV(1,3.)
            CALL PLTDV2(TEXTP(14),1,0.,.5,4.,.5)
            TEXTP(18) = TEXTP(18) + 4.*TEXTP(28)*YCHRSZ*TEXTP(31)
            TEXTP(19) = TEXTP(19) + 4.*TEXTP(29)*YCHRSZ*TEXTP(31)
            LASCHR = 'M'
            CALL PLTREV
            GO TO 2030

         ELSE IF (ESC.EQ.'DLINE') THEN
            CALL PLTSVV
            CALL PLTSTV(1,2.)
            CALL PLTDV2(TEXTP(14),1,0.,.5,4.,.5)
            TEXTP(18) = TEXTP(18) + 4.*TEXTP(28)*YCHRSZ*TEXTP(31)
            TEXTP(19) = TEXTP(19) + 4.*TEXTP(29)*YCHRSZ*TEXTP(31)
            LASCHR = 'M'
            CALL PLTREV
            GO TO 2030

         ELSE IF (ESC.EQ.'LDLINE') THEN
            CALL PLTSVV
            CALL PLTSTV(1,5.)
            CALL PLTDV2(TEXTP(14),1,0.,.5,4.,.5)
            TEXTP(18) = TEXTP(18) + 4.*TEXTP(28)*YCHRSZ*TEXTP(31)
            TEXTP(19) = TEXTP(19) + 4.*TEXTP(29)*YCHRSZ*TEXTP(31)
            LASCHR = 'M'
            CALL PLTREV
            GO TO 2030

         ELSE IF (ESC.EQ.'MDLINE') THEN
            CALL PLTSVV
            CALL PLTSTV(1,6.)
            CALL PLTDV2(TEXTP(14),1,0.,.5,4.,.5)
            TEXTP(18) = TEXTP(18) + 4.*TEXTP(28)*YCHRSZ*TEXTP(31)
            TEXTP(19) = TEXTP(19) + 4.*TEXTP(29)*YCHRSZ*TEXTP(31)
            LASCHR = 'M'
            CALL PLTREV
            GO TO 2030

         ELSE IF (ESC.EQ.'SDLINE') THEN
            CALL PLTSVV
            CALL PLTSTV(1,4.)
            CALL PLTDV2(TEXTP(14),1,0.,.5,4.,.5)
            TEXTP(18) = TEXTP(18) + 4.*TEXTP(28)*YCHRSZ*TEXTP(31)
            TEXTP(19) = TEXTP(19) + 4.*TEXTP(29)*YCHRSZ*TEXTP(31)
            LASCHR = 'M'
            CALL PLTREV
            GO TO 2030

         ELSE IF (ESC.EQ.'SLINE') THEN
            CALL PLTSVV
            CALL PLTSTV(1,1.)
            CALL PLTDV2(TEXTP(14),1,0.,.5,4.,.5)
            TEXTP(18) = TEXTP(18) + 4.*TEXTP(28)*YCHRSZ*TEXTP(31)
            TEXTP(19) = TEXTP(19) + 4.*TEXTP(29)*YCHRSZ*TEXTP(31)
            LASCHR = 'M'
            CALL PLTREV
            GO TO 2030

         ELSE
            STATUS = CHRCI(ESC,IESC)
            IF (STATUS) THEN
               ASCII = IESC

            ELSE
               CALL PLTRIM(ESC,L)
               CALL PLTFLU
               CALL SIORPT('PLTXTS','Invalid escape sequence "'//
     *                     ESC(1:L)//'"; escape sequence ignored.',2)
               GO TO 2030

            END IF

         END IF

      END IF

      NOVECT = NVECT(ASCII,IFONT)
      J = 0
 2050 IF (.NOT. (J.LT.NOVECT)) GO TO 2060
      JN = MIN(32,NOVECT-J)
      CALL PLTDV2(TEXTP(14),JN,X0(IDEX(ASCII,IFONT)+J,IFONT),
     *            Y0(IDEX(ASCII,IFONT)+J,IFONT),
     *            X1(IDEX(ASCII,IFONT)+J,IFONT),
     *            Y1(IDEX(ASCII,IFONT)+J,IFONT))
      J = J + JN
      GO TO 2050

 2060 CONTINUE
      TEXTP(18) = TEXTP(18) + XSIZE(ASCII,IFONT)*TEXTP(28)*YCHRSZ*
     *            TEXTP(31)
      TEXTP(19) = TEXTP(19) + XSIZE(ASCII,IFONT)*TEXTP(29)*YCHRSZ*
     *            TEXTP(31)
      IF (I.LE.LEN(TEXT)) THEN
         LASCHR = TEXT(I:I)
      END IF

      IF (CPUIFC(.FALSE.)) THEN
         GO TO 2040

      END IF

 2030 I = I + 1
      GO TO 2020

 2040 CONTINUE
      CALL PLTMOV(X,Y)
      TEXTP(8) = TEXTP(18)
      TEXTP(9) = TEXTP(19)
      TEXTP(10) = X
      TEXTP(11) = YSAVE
      IF (IFLAG.NE.0) THEN
         DO 2070 I = 14,17
            TEXTP(I) = TEXTP(I)/TEXTP(32)
 2070    CONTINUE
      END IF

      CALL PLTREV
      RETURN

      END
