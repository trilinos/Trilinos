      CHARACTER*16 FUNCTION ENGNOT (RNUM, IPREC)
      REAL RNUM

      CHARACTER*10 TFORM
      CHARACTER*32 SCRSTR

      NSIG = IPREC
      NSIG = MIN(8, MAX(2, NSIG)) - 1

         ENGNOT = ' '
         DO 10 I=1,3
            WRITE (TFORM,  30) I, NSIG+I+7, NSIG+I
            WRITE (SCRSTR, TFORM) RNUM
            READ  (SCRSTR(NSIG+I+6:NSIG+I+7), 40) IEXP
            IF (MOD(IEXP, 3) .EQ. 0) THEN
               ENGNOT(9-(NSIG+I):16) = SCRSTR(:LENSTR(SCRSTR))
               RETURN
            END IF
   10    CONTINUE
      RETURN
   30 FORMAT ('(',I1,'PE',I2.2,'.',I2.2,')')
   40 FORMAT (I2)
      END
