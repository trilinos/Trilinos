C    Copyright(C) 2008-2017 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
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
      SUBROUTINE STRIPB( STRING,ILEFT,IRIGHT )
      CHARACTER*(*) STRING
C
************************************************************************
C
C     FREFLD INPUT SYSTEM - ANSI FORTRAN - UTILITY ROUTINE
C
C     DESCRIPTION:
C     This routine strips leading and trailing blanks from a string. It
C     does not modify or copy the string, but simply returns the
C     location of the first and last non-blank characters.  If the
C     string is completely blank, ILEFT=LEN(STRING)+1 and IRIGHT=0 will
C     be returned.
C
C     FORMAL PARAMETERS:
C     STRING  CHARACTER  Any character string.
C     ILEFT   INTEGER    Position of first non-blank character.
C     IRIGHT  INTEGER    Position of last non-blank character.
C
************************************************************************
C
C ... Needed for 64-bit Solaris compile. Arg to CHAR must be
C     an integer of correct size.

      INTEGER TABC
      TABC = 9

C Get length of the string -
      LS = LEN( STRING )
C
C Find the first non-blank character -
      DO 10 N = 1 , LS
         IF ( STRING(N:N) .NE. ' '
     $        .AND. STRING(N:N) .NE. CHAR(TABC) ) GO TO 15
   10 CONTINUE
   15 CONTINUE
      ILEFT = N
C
C Find the last non-blank character -
      DO 20 N = LS , 1 , -1
         IF ( STRING(N:N) .NE. ' '
     $        .AND. STRING(N:N) .NE. CHAR(TABC) ) GO TO 25
   20 CONTINUE
   25 CONTINUE
      IRIGHT = N
C
      END
