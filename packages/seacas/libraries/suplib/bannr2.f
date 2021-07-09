C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE BANNR2 (NCOLS, LINEIN, IOUT)
C=======================================================================
      IMPLICIT INTEGER (A-Z)
      PARAMETER (MAXCHR = 14)
      CHARACTER*8 LETTER(7,49), SECT(7,MAXCHR)
      CHARACTER*49 MATRIX
      CHARACTER*(*) LINEIN
      CHARACTER*(MAXCHR) LINE
      CHARACTER BLANK*66
      SAVE MATRIX, LETTER, BLANK
      DATA BLANK /' '/
      DATA MATRIX(1:36)  /'ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789'/
      DATA MATRIX(37:48) /'$()*+,/.: -='/
      DATA MATRIX(49:49) /''''/

C     THE FOLLOWING CHARACTER SET IS NOT ANSI STANDARD

C      DATA MATRIX(37:68) /'!"#$%&()*+,/.:;<>?@[\]^_`{|}~ -='/
C      DATA MATRIX(69:69) /''''/

      DATA (LETTER(I, 1),I=1,7) /
     *   ' AAAAA  ',
     *   'AA   AA ',
     *   'AA   AA ',
     *   'AAAAAAA ',
     *   'AA   AA ',
     *   'AA   AA ',
     *   'AA   AA '/
      DATA (LETTER(I, 2),I=1,7) /
     *   'BBBBBB  ',
     *   'BB   BB ',
     *   'BB   BB ',
     *   'BBBBBB  ',
     *   'BB   BB ',
     *   'BB   BB ',
     *   'BBBBBB  '/
      DATA (LETTER(I, 3),I=1,7) /
     *   ' CCCCC  ',
     *   'CC   CC ',
     *   'CC      ',
     *   'CC      ',
     *   'CC      ',
     *   'CC   CC ',
     *   ' CCCCC  '/
      DATA (LETTER(I, 4),I=1,7) /
     *   'DDDDDD  ',
     *   'DD   DD ',
     *   'DD   DD ',
     *   'DD   DD ',
     *   'DD   DD ',
     *   'DD   DD ',
     *   'DDDDDD  '/
      DATA (LETTER(I, 5),I=1,7) /
     *   'EEEEEEE ',
     *   'EE      ',
     *   'EE      ',
     *   'EEEEE   ',
     *   'EE      ',
     *   'EE      ',
     *   'EEEEEEE '/
      DATA (LETTER(I, 6),I=1,7) /
     *   'FFFFFFF ',
     *   'FF      ',
     *   'FF      ',
     *   'FFFFF   ',
     *   'FF      ',
     *   'FF      ',
     *   'FF      '/
      DATA (LETTER(I, 7),I=1,7) /
     *   ' GGGGG  ',
     *   'GG   GG ',
     *   'GG      ',
     *   'GG      ',
     *   'GG  GGG ',
     *   'GG   GG ',
     *   ' GGGGG  '/
      DATA (LETTER(I, 8),I=1,7) /
     *   'HH   HH ',
     *   'HH   HH ',
     *   'HH   HH ',
     *   'HHHHHHH ',
     *   'HH   HH ',
     *   'HH   HH ',
     *   'HH   HH '/
      DATA (LETTER(I, 9),I=1,7) /
     *   '  IIII  ',
     *   '   II   ',
     *   '   II   ',
     *   '   II   ',
     *   '   II   ',
     *   '   II   ',
     *   '  IIII  '/
      DATA (LETTER(I,10),I=1,7) /
     *   '     JJ ',
     *   '     JJ ',
     *   '     JJ ',
     *   '     JJ ',
     *   '     JJ ',
     *   'JJ   JJ ',
     *   ' JJJJJ  '/
      DATA (LETTER(I,11),I=1,7) /
     *   'KK   KK ',
     *   'KK  KK  ',
     *   'KK KK   ',
     *   'KKKK    ',
     *   'KKKKK   ',
     *   'KK  KK  ',
     *   'KK   KK '/
      DATA (LETTER(I,12),I=1,7) /
     *   'LL      ',
     *   'LL      ',
     *   'LL      ',
     *   'LL      ',
     *   'LL      ',
     *   'LL      ',
     *   'LLLLLLL '/
      DATA (LETTER(I,13),I=1,7) /
     *   'M     M ',
     *   'MM   MM ',
     *   'MMM MMM ',
     *   'MM M MM ',
     *   'MM   MM ',
     *   'MM   MM ',
     *   'MM   MM '/
      DATA (LETTER(I,14),I=1,7) /
     *   'N    NN ',
     *   'NN   NN ',
     *   'NNN  NN ',
     *   'NN N NN ',
     *   'NN  NNN ',
     *   'NN   NN ',
     *   'NN    N '/
      DATA (LETTER(I,15),I=1,7) /
     *   ' OOOOO  ',
     *   'OO   OO ',
     *   'OO   OO ',
     *   'OO   OO ',
     *   'OO   OO ',
     *   'OO   OO ',
     *   ' OOOOO  '/
      DATA (LETTER(I,16),I=1,7) /
     *   'PPPPPP  ',
     *   'PP   PP ',
     *   'PP   PP ',
     *   'PPPPPP  ',
     *   'PP      ',
     *   'PP      ',
     *   'PP      '/
      DATA (LETTER(I,17),I=1,7) /
     *   ' QQQQQ  ',
     *   'QQ   QQ ',
     *   'QQ   QQ ',
     *   'QQ   QQ ',
     *   'QQ Q QQ ',
     *   'QQ  QQQ ',
     *   ' QQQQQQ '/
      DATA (LETTER(I,18),I=1,7) /
     *   'RRRRRR  ',
     *   'RR   RR ',
     *   'RR   RR ',
     *   'RRRRRR  ',
     *   'RRRRR   ',
     *   'RR  RR  ',
     *   'RR   RR '/
      DATA (LETTER(I,19),I=1,7) /
     *   ' SSSSSS ',
     *   'SS      ',
     *   'SS      ',
     *   ' SSSSS  ',
     *   '     SS ',
     *   '     SS ',
     *   'SSSSSS  '/
      DATA (LETTER(I,20),I=1,7) /
     *   'TTTTTT  ',
     *   '  TT    ',
     *   '  TT    ',
     *   '  TT    ',
     *   '  TT    ',
     *   '  TT    ',
     *   '  TT    '/
      DATA (LETTER(I,21),I=1,7) /
     *   'UU   UU ',
     *   'UU   UU ',
     *   'UU   UU ',
     *   'UU   UU ',
     *   'UU   UU ',
     *   'UU   UU ',
     *   ' UUUUU  '/
      DATA (LETTER(I,22),I=1,7) /
     *   'V      V',
     *   ' V    V ',
     *   ' VV  VV ',
     *   '  V  V  ',
     *   '  VVVV  ',
     *   '   VV   ',
     *   '   VV   '/
      DATA (LETTER(I,23),I=1,7) /
     *   'WW   WW ',
     *   'WW   WW ',
     *   'WW   WW ',
     *   'WW   WW ',
     *   'WW W WW ',
     *   'WWWWWWW ',
     *   ' WW WW  '/
      DATA (LETTER(I,24),I=1,7) /
     *   'XX   XX ',
     *   ' XX XX  ',
     *   '  XXX   ',
     *   '  XXX   ',
     *   '  XXX   ',
     *   ' XX XX  ',
     *   'XX   XX '/
      DATA (LETTER(I,25),I=1,7) /
     *   'YY   YY ',
     *   ' YY YY  ',
     *   '  YYY   ',
     *   '   YY   ',
     *   '   YY   ',
     *   '   YY   ',
     *   '   YY   '/
      DATA (LETTER(I,26),I=1,7) /
     *   'ZZZZZZZ ',
     *   '     Z  ',
     *   '    Z   ',
     *   '   Z    ',
     *   '  Z     ',
     *   ' Z      ',
     *   'ZZZZZZZ '/
      DATA (LETTER(I,27),I=1,7) /
     *   ' 000000 ',
     *   '0    00 ',
     *   '0   0 0 ',
     *   '0  0  0 ',
     *   '0 0   0 ',
     *   '00    0 ',
     *   '000000  '/
      DATA (LETTER(I,28),I=1,7) /
     *   '   1    ',
     *   '  11    ',
     *   ' 1 1    ',
     *   '   1    ',
     *   '   1    ',
     *   '   1    ',
     *   ' 11111  '/
      DATA (LETTER(I,29),I=1,7) /
     *   '  2222  ',
     *   ' 2    2 ',
     *   '     2  ',
     *   '    2   ',
     *   '   2    ',
     *   '  2     ',
     *   ' 222222 '/
      DATA (LETTER(I,30),I=1,7) /
     *   ' 33333  ',
     *   '3     3 ',
     *   '      3 ',
     *   '    33  ',
     *   '      3 ',
     *   '3     3 ',
     *   ' 33333  '/
      DATA (LETTER(I,31),I=1,7) /
     *   '    44  ',
     *   '  4  4  ',
     *   ' 4   4  ',
     *   '444444  ',
     *   '     4  ',
     *   '     4  ',
     *   '     4  '/
      DATA (LETTER(I,32),I=1,7) /
     *   '555555  ',
     *   '5       ',
     *   '5       ',
     *   '55555   ',
     *   '     5  ',
     *   '     5  ',
     *   '55555   '/
      DATA (LETTER(I,33),I=1,7) /
     *   ' 6666   ',
     *   '6       ',
     *   '6       ',
     *   '66666   ',
     *   '6    6  ',
     *   '6    6  ',
     *   ' 6666   '/
      DATA (LETTER(I,34),I=1,7) /
     *   '7777777 ',
     *   '     7  ',
     *   '    7   ',
     *   '   7    ',
     *   '  7     ',
     *   ' 7      ',
     *   '7       '/
      DATA (LETTER(I,35),I=1,7) /
     *   '  8888  ',
     *   ' 8    8 ',
     *   ' 8    8 ',
     *   '  8888  ',
     *   ' 8    8 ',
     *   ' 8    8 ',
     *   '  8888  '/
      DATA (LETTER(I,36),I=1,7) /
     *   '  9999  ',
     *   ' 9    9 ',
     *   ' 9    9 ',
     *   '  99999 ',
     *   '      9 ',
     *   '      9 ',
     *   '  9999  '/
C      DATA (LETTER(I,37),I=1,7) /
C     * '   !    ',
C     * '   !    ',
C     * '   !    ',
C     * '   !    ',
C     * '   !    ',
C     * '        ',
C     * '   !    '/
C      DATA (LETTER(I,38),I=1,7) /
C     * '  " "   ',
C     * '  " "   ',
C     * '  " "   ',
C     * '        ',
C     * '        ',
C     * '        ',
C     * '        '/
C      DATA (LETTER(I,39),I=1,7) /
C     * '  # #   ',
C     * '  # #   ',
C     * '####### ',
C     * '  # #   ',
C     * '####### ',
C     * '  # #   ',
C     * '  # #   '/
      DATA (LETTER(I,37),I=1,7) /
     *   '   $    ',
     *   ' $$$$$$ ',
     *   '$  $    ',
     *   ' $$$$$  ',
     *   '   $  $ ',
     *   '$$$$$$  ',
     *   '   $    '/
C      DATA (LETTER(I,41),I=1,7) /
C     * '%%    % ',
C     * '%%   %  ',
C     * '    %   ',
C     * '   %    ',
C     * '  %     ',
C     * ' %   %% ',
C     * '%    %% '/
C      DATA (LETTER(I,42),I=1,7) /
C     * '  &     ',
C     * ' & &    ',
C     * '  &     ',
C     * ' & &    ',
C     * '&   & & ',
C     * '&    &  ',
C     * ' &&&& & '/
      DATA (LETTER(I,38),I=1,7) /
     *   '     (  ',
     *   '    (   ',
     *   '   (    ',
     *   '   (    ',
     *   '   (    ',
     *   '    (   ',
     *   '     (  '/
      DATA (LETTER(I,39),I=1,7) /
     *   ' )      ',
     *   '  )     ',
     *   '   )    ',
     *   '   )    ',
     *   '   )    ',
     *   '  )     ',
     *   ' )      '/
      DATA (LETTER(I,40),I=1,7) /
     *   '*     * ',
     *   ' *   *  ',
     *   '  * *   ',
     *   '******* ',
     *   '  * *   ',
     *   ' *   *  ',
     *   '*     * '/
      DATA (LETTER(I,41),I=1,7) /
     *   '   +    ',
     *   '   +    ',
     *   '   +    ',
     *   '+++++++ ',
     *   '   +    ',
     *   '   +    ',
     *   '   +    '/
      DATA (LETTER(I,42),I=1,7) /
     *   '        ',
     *   '        ',
     *   '        ',
     *   '        ',
     *   '  ,,    ',
     *   '   ,    ',
     *   '  ,     '/
      DATA (LETTER(I,43),I=1,7) /
     *   '      / ',
     *   '     /  ',
     *   '    /   ',
     *   '   /    ',
     *   '  /     ',
     *   ' /      ',
     *   '/       '/
      DATA (LETTER(I,44),I=1,7) /
     *   '        ',
     *   '        ',
     *   '        ',
     *   '        ',
     *   '        ',
     *   '  ..    ',
     *   '  ..    '/
      DATA (LETTER(I,45),I=1,7) /
     *   '        ',
     *   '  ::    ',
     *   '  ::    ',
     *   '        ',
     *   '  ::    ',
     *   '  ::    ',
     *   '        '/
C      DATA (LETTER(I,51),I=1,7) /
C     * '        ',
C     * '  ;;    ',
C     * '  ;;    ',
C     * '        ',
C     * '  ;;    ',
C     * '   ;    ',
C     * '  ;     '/
C      DATA (LETTER(I,52),I=1,7) /
C     * '    <   ',
C     * '   <    ',
C     * '  <     ',
C     * ' <      ',
C     * '  <     ',
C     * '   <    ',
C     * '    <   '/
C      DATA (LETTER(I,53),I=1,7) /
C     * ' >      ',
C     * '  >     ',
C     * '   >    ',
C     * '    >   ',
C     * '   >    ',
C     * '  >     ',
C     * ' >      '/
C      DATA (LETTER(I,54),I=1,7) /
C     * '  ???   ',
C     * ' ?   ?  ',
C     * '     ?  ',
C     * '    ?   ',
C     * '   ?    ',
C     * '        ',
C     * '   ?    '/
C      DATA (LETTER(I,55),I=1,7) /
C     * ' @@@@@  ',
C     * '@     @ ',
C     * '      @ ',
C     * ' @@@  @ ',
C     * '@  @  @ ',
C     * '@  @ @  ',
C     * ' @@@@   '/
C      DATA (LETTER(I,56),I=1,7) /
C     * '  [[[   ',
C     * '  [     ',
C     * '  [     ',
C     * '  [     ',
C     * '  [     ',
C     * '  [     ',
C     * '  [[[   '/
C      DATA (LETTER(I,57),I=1,7) /
C     * '\       ',
C     * ' \      ',
C     * '  \     ',
C     * '   \    ',
C     * '    \   ',
C     * '     \  ',
C     * '      \ '/
C      DATA (LETTER(I,58),I=1,7) /
C     * ' ]]]    ',
C     * '   ]    ',
C     * '   ]    ',
C     * '   ]    ',
C     * '   ]    ',
C     * '   ]    ',
C     * ' ]]]    '/
C      DATA (LETTER(I,59),I=1,7) /
C     * '   ^    ',
C     * '  ^ ^   ',
C     * ' ^   ^  ',
C     * '^     ^ ',
C     * '        ',
C     * '        ',
C     * '        '/
C      DATA (LETTER(I,60),I=1,7) /
C     * '        ',
C     * '        ',
C     * '        ',
C     * '        ',
C     * '        ',
C     * '        ',
C     * '_______ '/
C      DATA (LETTER(I,61),I=1,7) /
C     * '   ``   ',
C     * '   `    ',
C     * '    `   ',
C     * '        ',
C     * '        ',
C     * '        ',
C     * '        '/
C      DATA (LETTER(I,62),I=1,7) /
C     * '    {   ',
C     * '   {    ',
C     * '   {    ',
C     * '  {     ',
C     * '   {    ',
C     * '   {    ',
C     * '    {   '/
C      DATA (LETTER(I,63),I=1,7) /
C     * '   |    ',
C     * '   |    ',
C     * '   |    ',
C     * '   |    ',
C     * '   |    ',
C     * '   |    ',
C     * '   |    '/
C      DATA (LETTER(I,64),I=1,7) /
C     * '  }     ',
C     * '   }    ',
C     * '   }    ',
C     * '    }   ',
C     * '   }    ',
C     * '   }    ',
C     * '  }     '/
C      DATA (LETTER(I,65),I=1,7) /
C     * ' ~~     ',
C     * '~  ~  ~ ',
C     * '    ~~  ',
C     * '        ',
C     * '        ',
C     * '        ',
C     * '        '/
      DATA (LETTER(I,46),I=1,7) /
     *   '        ',
     *   '        ',
     *   '        ',
     *   '        ',
     *   '        ',
     *   '        ',
     *   '        '/
      DATA (LETTER(I,47),I=1,7) /
     *   '        ',
     *   '        ',
     *   '        ',
     *   '------- ',
     *   '        ',
     *   '        ',
     *   '        '/
      DATA (LETTER(I,48),I=1,7) /
     *   '        ',
     *   '        ',
     *   '======= ',
     *   '        ',
     *   '======= ',
     *   '        ',
     *   '        '/
      DATA (LETTER(I,49),I=1,7) /
     *   '  ''''    ',
     *   '   ''    ',
     *   '  ''     ',
     *   '        ',
     *   '        ',
     *   '        ',
     *   '        '/

      MAXCOL = MIN(NCOLS/9, MAXCHR)

C     DELIMIT NONBLANK STRING.

      CALL STRIPB (LINEIN, ILEFT, IRIGHT)
      IF (ILEFT .LE. IRIGHT) THEN
         LINE = LINEIN (ILEFT:IRIGHT)
         LENIN = IRIGHT - ILEFT + 1
         IF (LENIN .GT. MAXCOL) THEN
            LENIN = MAXCOL
            CALL STRIPB (LINE(:LENIN), J, LENIN)
         END IF
      ELSE
         LINE = ' '
         LENIN = 0
      END IF

C     LENIN IS LAST PRINTABLE NONBLANK

C     CONVERT ALPHABET TO UPPER CASE

      DO 100 J=1,LENIN
         IF (LGE(LINE(J:J),'a') .AND. LLE(LINE(J:J),'z')) THEN
            ITEMP = ICHAR(LINE(J:J))
            LINE(J:J)=CHAR(ITEMP-(ICHAR('a')-ICHAR('A')))
         END IF
  100 CONTINUE

C     CALCULATE BLANK FILL.

      NBLANK = (NCOLS - LENIN * 9) / 2
      NBLANK = MIN (NBLANK, 66)

C     LOAD UP CHARACTERS

      DO 130 ICOL = 1, LENIN
         IPT = INDEX(MATRIX,LINE(ICOL:ICOL))
         IF (IPT .EQ. 0) THEN

C           CHARACTER NOT FOUND - REPLACE WITH A BLANK

            DO 110 IROW = 1, 7
               SECT(IROW,ICOL) = ' '
  110       CONTINUE

         ELSE

C           CHARACTER FOUND - INSERT BANNER LETTER

            DO 120 IROW = 1, 7
               SECT(IROW,ICOL) = LETTER(IROW,IPT)
  120       CONTINUE

         END IF
  130 CONTINUE

      IF ((IRIGHT - ILEFT + 1) .NE. LENIN .AND. LENIN .NE. 0) THEN

C        STRING IS TRUNCATED.

        if (iout .eq. 0) then
          WRITE (*, 5010) LINEIN(ILEFT:IRIGHT)
        else
          WRITE (IOUT, 5010) LINEIN(ILEFT:IRIGHT)
        end if

      ELSE

C        STRING IS NOT TRUNCATED OR IS NULL.

        if (iout .eq. 0) then
          WRITE (*,5000)
        else
          WRITE (IOUT,5000)
        end if

      END IF
      if (iout .eq. 0) then
        DO 140 IROW = 1, 7
          WRITE(*,5020)BLANK(:NBLANK),(SECT(IROW,J), J = 1, LENIN)
 140    CONTINUE
        WRITE (*,5000)
      else
        DO 150 IROW = 1, 7
          WRITE(IOUT,5020)BLANK(:NBLANK),(SECT(IROW,J), J = 1, LENIN)
 150    CONTINUE
        WRITE (IOUT,5000)
      end if
      RETURN
 5000 FORMAT ()
 5010 FORMAT(' WARNING, TRUNCATED BANNER STRING: ',A)
 5020 FORMAT (25(:,1X,A))
      END
