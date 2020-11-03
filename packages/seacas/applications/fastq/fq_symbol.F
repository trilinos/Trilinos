C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE SYMBOL (IDIM, X, Y, SYM)
C***********************************************************************

C  SUBROUTINE SYMBOL = DRAWS SYMBOLS AT X,Y ARRAY LOCATIONS

C***********************************************************************

      CHARACTER*6 SYM

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

      RETURN

      END
