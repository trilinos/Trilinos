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
