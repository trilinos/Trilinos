C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE PLTXSL(CHARST,LENGTH)
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
      CHARACTER*(*) CHARST
      CHARACTER ESC*20,ESCCHR*1, textip1*1
      REAL LENGTH,LENGT1
      INTEGER ASCII
      LOGICAL STATUS,CHRCI

      ESCCHR = CHAR(92)
      LENGTH = 0.
      LENGT1 = 0.
      YSZE = TEXTP(1)
      CALL PLTRIM(CHARST,NUM)
      MODE = 0
      IFONT = 1
      I = 1
 2210 IF (.NOT. (I.LE.NUM)) GO TO 2230
      ASCII = ICHAR(CHARST(I:I))
      IF (ASCII.LT.1 .OR. ASCII.GT.126) THEN
         CALL CHRIC(ASCII,ESCCHR,LI)
         CALL PLTFLU
         RETURN

      END IF

      RLINE = 0.
      if (i .lt. num) then
         textip1 = charst(i+1:i+1)
      else
         textip1 = escchr
      end if
      IF (ASCII.EQ.ICHAR(ESCCHR) .AND. textip1.EQ.ESCCHR) THEN
         I = I + 1

      ELSE IF (ASCII.EQ.ICHAR(ESCCHR)) THEN
         CALL PLTESC(CHARST,I,ESC)
         CALL CHRUP(ESC,ESC)
         IF (ESC.EQ.'^' .OR. ESC.EQ.'_') THEN
            IF (MODE.EQ.0) THEN
               YSZE = TEXTP(1)*TEXTP(32)
               MODE = 1
            END IF

            GO TO 2220

         ELSE IF (ESC.EQ.'-') THEN
            YSZE = TEXTP(1)
            MODE = 0
            GO TO 2220

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
            IF (LENGTH.GT.LENGT1) THEN
               LENGT1 = LENGTH
            END IF

            LENGTH = 0.
            GO TO 2220

         ELSE IF (ESC.EQ.'LF') THEN
            IF (LENGTH.GT.LENGT1) THEN
               LENGT1 = LENGTH
            END IF

            LENGTH = 0.
            GO TO 2220

         ELSE IF (ESC.EQ.'CL') THEN
            IF (LENGTH.GT.LENGT1) THEN
               LENGT1 = LENGTH
            END IF

            LENGTH = 0.
            GO TO 2220

         ELSE IF (ESC.EQ.'ENG') THEN
            IFONT = 1
            GO TO 2220

         ELSE IF (ESC.EQ.'GR') THEN
            IFONT = 2
            GO TO 2220

         ELSE IF (ESC.EQ.'DDLINE') THEN
            RLINE = 4.

         ELSE IF (ESC.EQ.'DLINE') THEN
            RLINE = 4.

         ELSE IF (ESC.EQ.'LDLINE') THEN
            RLINE = 4.

         ELSE IF (ESC.EQ.'MDLINE') THEN
            RLINE = 4.

         ELSE IF (ESC.EQ.'SDLINE') THEN
            RLINE = 4.

         ELSE IF (ESC.EQ.'SLINE') THEN
            RLINE = 4.

         ELSE
            STATUS = CHRCI(ESC,IESC)
            IF (STATUS) THEN
               ASCII = IESC

            ELSE
               CALL PLTRIM(ESC,L)
               CALL PLTFLU
               CALL SIORPT('PLTXTS','Invalid escape sequence "'//
     *                     ESC(1:L)//'"; escape sequence ignored.',2)
               GO TO 2220

            END IF

         END IF

      END IF

      IF (RLINE.EQ.4.) THEN
         LENGTH = LENGTH + 4.*YSZE**TEXTP(31)

      ELSE
         LENGTH = LENGTH + XSIZE(ASCII,IFONT)*YSZE*TEXTP(31)
      END IF

 2220 I = I + 1
      GO TO 2210

 2230 CONTINUE
      IF (LENGT1.GT.LENGTH) THEN
         LENGTH = LENGT1
      END IF

      IF (NUM.EQ.1 .AND. RLINE.NE.4. .AND. ASCII.GT.32) THEN
         TEMP = TEXTP(39)/TEXTP(40)
         LENGTH = LENGTH - TEMP*YSZE*TEXTP(31)
      END IF

      RETURN

      END
