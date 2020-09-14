C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

      SUBROUTINE PLTONE (MD, TITLE, NTITLE, XX, YY, XLAB, YLAB, NP,
     1   CURVE, TPLOT)

C        OUTPUT DATA IN GRAFAID NEUTRAL FILE FORMAT FOR ONE PLOT

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

C ARRAYS:
C     Name     Dimension   Type     Description
C     ----     ---------   ----     ----------------------------------
C     TITLE    NTITLE      char*80  Title for current curve
C     XX       NP          real     X-variable data
C     YY       NP          real     Y-variable data

C     SUBROUTINES AND FUNCTIONS CALLED:
C        PACKT -  Remove multiple blanks from a character string
C                 The string:     "This   is    the    title"
C                 Is returned as: "This is the title"

      CHARACTER*(*) TITLE(*)
      CHARACTER*(*) XLAB,YLAB
      DIMENSION XX(*), YY(*)

      LOGICAL MONO, TPLOT
      CHARACTER*(*) CURVE
      CHARACTER*11 BEGIN
      CHARACTER*9 ECURVE
      CHARACTER*1 COMMA,AUX
      CHARACTER*4 XTYP, AXTYP
      DATA COMMA/','/, AXTYP/'NOLO'/
      DATA BEGIN/'BEGIN CURVE'/, ECURVE/'END CURVE'/, AUX/'F'/

C     ...LOCATE MINIMUM AND MAXIMUM VALUES AND CHECK FOR NONMONOTONIC DATA

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

      YMN=YY(1)
      YMX=YY(1)
      DO 20 I=2,NP
         YMN = MIN(YMN, YY(I))
         YMX = MAX(YMX, YY(I))
   20    CONTINUE

C     BEGIN TO WRITE CURVE PACKET

      WRITE (MD, 40) BEGIN,COMMA,CURVE
      CALL PACKT (TITLE(1),80)
      WRITE (MD, 50) NTITLE,COMMA,TITLE(1)
      DO 30 I=2,NTITLE
         CALL PACKT (TITLE(I),80)
         WRITE (MD, 60) TITLE(I)
   30    CONTINUE
      CALL PACKT (XLAB,lenstr(xlab))
      CALL PACKT (YLAB,lenstr(ylab))
      WRITE (MD, 70) XLAB
      WRITE (MD, 70) YLAB
      WRITE (MD, 80) XMN,COMMA,XMX,COMMA,YMN,COMMA,YMX,COMMA,
     1               NP,COMMA,AUX
      WRITE (MD, 90) AXTYP,COMMA,XTYP,COMMA

C     WRITE DATA PAIRS

      DO 35 III=1,NP
         WRITE (MD, 100) XX(III),COMMA,YY(III)
   35  CONTINUE

C     WRITE END OF CURVE PACKET

      WRITE (MD, 110) ECURVE,COMMA,CURVE
      RETURN

   40 FORMAT (A11,A1,A16)
   50 FORMAT (I1,A1,A80)
   60 FORMAT (A)
   70 FORMAT (A)
   80 FORMAT (1PE15.7,A1,1PE15.7,A1,1PE15.7,A1,1PE15.7,A1,I5,A1,A1)
   90 FORMAT (A4,A1,A4,A1)
  100 FORMAT (1PE15.7,A1,1PE15.7)
  110 FORMAT (A9,A1,A16)
      END
