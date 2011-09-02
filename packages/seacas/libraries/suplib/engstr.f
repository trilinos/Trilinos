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
      SUBROUTINE ENGSTR (NNUM, NSIG, RNUM, RSTR, LSTR)
C=======================================================================
C$Id: engstr.f,v 1.2 2009/03/25 12:46:01 gdsjaar Exp $
C$Log: engstr.f,v $
CRevision 1.2  2009/03/25 12:46:01  gdsjaar
CAdd copyright and license notice to all files.
C
CRevision 1.1.1.1  1990/08/14 16:14:04  gdsjaar
CTesting
C
c Revision 1.1  90/08/14  16:14:03  gdsjaar
c Initial revision
c 
c Revision 1.1  90/08/09  13:39:21  gdsjaar
c Initial revision
c 

C   --*** ENGSTR *** (STRLIB) Convert real numbers to strings
C   --   Written by Amy Gilkey - revised 02/14/86
C   --
C   --ENGSTR converts a set of real numbers into engineering notation
C   --strings.  The real numbers are not compared in any way to get a
C   --consistent set of exponents, etc.
C   --
C   --Parameters:
C   --   NNUM - IN - the number of numbers in the array set
C   --   NSIG - IN - the maximum number of significant digits, max of 8
C   --   RNUM - IN - the array of real numbers to be converted
C   --   RSTR - OUT - the set of real number strings
C   --   LSTR - OUT - the maximum length of the number strings

C   --Routines Called:
C   --   IENGRX - (included) Get engineering notation exponent

      INTEGER NNUM
      INTEGER NSIG
      REAL RNUM(*)
      CHARACTER*(*) RSTR(*)
      INTEGER LSTR

      CHARACTER*10 SCRFMT
      CHARACTER*20 SCRSTR
      CHARACTER*15 FFMT

C   --Do not try to use a common exponent, but use engineering notation

      WRITE (SCRFMT, 10, IOSTAT=IDUM) NSIG+7, NSIG
   10 FORMAT ('(0PE', I2.2, '.', I2.2, ')')

      LSTR = 0
      DO 40 I = 1, NNUM
         WRITE (SCRSTR(1:NSIG+7), SCRFMT, IOSTAT=IDUM) RNUM(I)
         READ (SCRSTR(NSIG+5:NSIG+7), '(I3)', IOSTAT=IDUM) IE
         ISIGN = 0
         IF (RNUM(I) .LT. 0.0) ISIGN = 1

         NEWEXP = IENGRX (IE, IE)

         EXPDIV = 10.0 ** NEWEXP

         NWHOLE = MAX (0, IE - NEWEXP)
         NFRAC = MAX (0, MIN (NEWEXP - IE + NSIG,
     &      NSIG - (IE - NEWEXP)))
         NTOTAL = ISIGN + NWHOLE + 1 + NFRAC

         WRITE (FFMT, 20, IOSTAT=IDUM) NTOTAL, NFRAC
   20    FORMAT ('(F', I2.2, '.', I2.2, ')')
         WRITE (FFMT(8:15), 30, IOSTAT=IDUM) NEWEXP
   30    FORMAT (',''E', SP, I3.2, ''')')
         LSTR = MAX (LSTR, NTOTAL+4)

         WRITE (RSTR(I), FFMT, IOSTAT=IDUM) RNUM(I)/EXPDIV
   40 CONTINUE

      RETURN
      END
