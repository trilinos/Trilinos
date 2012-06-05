C $Id: symbol.F,v 1.1 1998/04/10 04:30:25 gdsjaar Exp $
C $Log: symbol.F,v $
C Revision 1.1  1998/04/10 04:30:25  gdsjaar
C Added the ifdef needsdoubleescape to the graphics symbol
C plotting. Required moving some .f to .F files and specifying in Imakefile
C
C Revision 1.1.1.1  1990/11/30 11:17:02  gdsjaar
C FASTQ Version 2.0X
C
c Revision 1.1  90/11/30  11:17:01  gdsjaar
c Initial revision
c 
C
      SUBROUTINE SYMBOL (IDIM, X, Y, SYM)
C***********************************************************************
C
C  SUBROUTINE SYMBOL = DRAWS SYMBOLS AT X,Y ARRAY LOCATIONS
C
C***********************************************************************
C
      CHARACTER*6 SYM
C
      DIMENSION X(IDIM), Y(IDIM)
#if NeedsDoubleEscape
      IF (SYM .EQ. 'DIAMND') THEN
         CALL MPD2SY (IDIM, X, Y, '\\CDI')
      ELSEIF (SYM .EQ. 'CENTRX') THEN
         CALL MPD2SY (IDIM, X, Y, '\\CX')
      ELSEIF (SYM .EQ. 'TRINGL' )THEN
         CALL MPD2SY (IDIM, X, Y, '\\CTR')
      ELSEIF (SYM .EQ. 'CIRCLE' )THEN
         CALL MPD2SY (IDIM, X, Y, '\\CCI')
      ELSEIF (SYM .EQ. 'SQUARE' )THEN
         CALL MPD2SY (IDIM, X, Y, '\\CSQ')
      ENDIF
#else
      IF (SYM .EQ. 'DIAMND') THEN
         CALL MPD2SY (IDIM, X, Y, '\CDI')
      ELSEIF (SYM .EQ. 'CENTRX') THEN
         CALL MPD2SY (IDIM, X, Y, '\CX')
      ELSEIF (SYM .EQ. 'TRINGL' )THEN
         CALL MPD2SY (IDIM, X, Y, '\CTR')
      ELSEIF (SYM .EQ. 'CIRCLE' )THEN
         CALL MPD2SY (IDIM, X, Y, '\CCI')
      ELSEIF (SYM .EQ. 'SQUARE' )THEN
         CALL MPD2SY (IDIM, X, Y, '\CSQ')
      ENDIF
#endif
C
      RETURN
C
      END
