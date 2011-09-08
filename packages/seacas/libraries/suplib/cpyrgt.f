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

C=======================================================================
      SUBROUTINE CPYRGT (NOUT, YEAR)
C=======================================================================
C $Id: cpyrgt.f,v 1.4 2009/03/25 14:31:47 gdsjaar Exp $
C $Log: cpyrgt.f,v $
C Revision 1.4  2009/03/25 14:31:47  gdsjaar
C Update copyright info
C
C Revision 1.3  2009/03/25 12:46:01  gdsjaar
C Add copyright and license notice to all files.
C
C Revision 1.2  1993/07/06 21:57:53  gdsjaar
C Updated copyright output information based on latest memo from Art Silva
C
c Revision 1.1  1992/05/13  16:57:30  gdsjaar
c Added routine to output copyright notice during execution
c

C   --*** CPYRGT *** (ETCLIB) Print copyright notice
C   --   Written by Greg Sjaardema - revised 5-13-92 - 
C   --
C   --CPYRGT prints the copyright notice  at the start of any program.  
C   --The copyright notice is printed to the standard output device or 
C   --an output file.
C   --
C   --Parameters:
C   --   NOUT - IN - the output file number; 0 if standard output device
C   --   YEAR - IN - the year of the copyright

C   --Routines Called:
C   --   LENSTR - (STRLIB) Find string length

      PARAMETER (NLIN = 3)
      INTEGER NOUT
      CHARACTER*(*) YEAR

      CHARACTER*80 BANR
      CHARACTER*40 BLANK
      CHARACTER*60 TEXT(NLIN)
      SAVE BLANK

      DATA BLANK / ' ' /
      DATA TEXT  /
     *  'Under the terms of Contract',
     *  'DE-AC04-94AL85000 with Sandia Corporation, the' ,
     *  'U.S. Government retains certain rights in this software.'/

      NCEN(LEN) = MAX (1, (80 - LEN + 1) / 2)

      write (banr, 100) year(:lenstr(year))
 100  format (' Copyright ', a, ' Sandia Corporation')
      CALL SQZSTR (BANR, LBANR)
      IF (NOUT .LE. 0) THEN
         WRITE (*, 10000) BLANK(:NCEN(LBANR+8)),
     &      '+++ ', BANR(:LBANR), ' +++'
      ELSE
         WRITE (NOUT, 10000) BLANK(:NCEN(LBANR+8)),
     &      '+++ ', BANR(:LBANR), ' +++'
      END IF
      do 120 i=1, nlin
      write (banr, 110) text(i)(:lenstr(text(i)))
 110  format (A)
      CALL SQZSTR (BANR, LBANR)
      IF (NOUT .LE. 0) THEN
         WRITE (*, 10000) BLANK(:NCEN(LBANR+8)),
     &      '+++ ', BANR(:LBANR), ' +++'
      ELSE
         WRITE (NOUT, 10000) BLANK(:NCEN(LBANR+8)),
     &      '+++ ', BANR(:LBANR), ' +++'
      END IF
 120  continue

      IF (NOUT .LE. 0) THEN
         WRITE (*, *)
      ELSE
         WRITE (NOUT, *)
      END IF
10000 format (8A)
      RETURN
      END
