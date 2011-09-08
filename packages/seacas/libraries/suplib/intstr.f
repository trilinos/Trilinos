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
      SUBROUTINE INTSTR (NNUM, LNGSTR, INUM, ISTR, LSTR)
C=======================================================================
C$Id: intstr.f,v 1.2 2009/03/25 12:46:02 gdsjaar Exp $
C$Log: intstr.f,v $
CRevision 1.2  2009/03/25 12:46:02  gdsjaar
CAdd copyright and license notice to all files.
C
CRevision 1.1.1.1  1990/08/14 16:15:14  gdsjaar
CTesting
C
c Revision 1.1  90/08/14  16:15:12  gdsjaar
c Initial revision
c 
c Revision 1.1  90/08/09  13:39:32  gdsjaar
c Initial revision
c 

C   --*** INTSTR *** (STRLIB) Convert integer numbers to strings
C   --   Written by Amy Gilkey - revised 03/14/88
C   --
C   --INTSTR converts a set of integer numbers into a consistent set of
C   --strings of the same length (right justified).
C   --
C   --Parameters:
C   --   NNUM - IN - the number of integers numbers in the set
C   --   LNGSTR - IN - the number of digits in the string; <= 0 if length
C   --      of longest integer
C   --   INUM - IN - the array of integers to be converted
C   --   ISTR - OUT - the set of integer strings
C   --   LSTR - OUT - the maximum length of the number strings

      INTEGER NNUM
      INTEGER LNGSTR
      INTEGER INUM(*)
      CHARACTER*(*) ISTR(*)
      INTEGER LSTR

      CHARACTER*20 SCISTR
      CHARACTER*5 FFMT

      IF ((NNUM .EQ. 1) .AND. (LNGSTR .LE. 0)) THEN

C      --Handle special case of single number, smallest length

         WRITE (FFMT, 10000, IOSTAT=IDUM) LEN (ISTR(1))

         WRITE (ISTR(1), FFMT, IOSTAT=IDUM) INUM(1)
         CALL SQZSTR (ISTR(1), LSTR)

      ELSE

C      --Determine smallest length of largest number string

         IF (LNGSTR .LE. 0) THEN
            MINE = 0
            MAXE = 0
            DO 100 I = 1, NNUM
               MINE = MIN (MINE, INUM(I))
               MAXE = MAX (MAXE, INUM(I))
  100       CONTINUE

            MM = MAX (IABS (MAXE), IABS (MINE))
            IF (MM .EQ. IABS (MINE)) MM = -MM
            WRITE (SCISTR, '(I12)', IOSTAT=IDUM) MM
            CALL SQZSTR (SCISTR, LSTR)

         ELSE
            LSTR = LNGSTR
         END IF

         LSTR = MIN (LSTR, LEN (ISTR(1)))

C      --Convert all integers to string of given length

         WRITE (FFMT, 10000, IOSTAT=IDUM) LSTR
10000     FORMAT ('(I', I2.2, ')')

         DO 110 I = 1, NNUM
            WRITE (ISTR(I), FFMT, IOSTAT=IDUM) INUM(I)
  110    CONTINUE
      END IF

      RETURN
      END
