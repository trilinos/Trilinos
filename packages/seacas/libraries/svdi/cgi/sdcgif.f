C Copyright (C) 2009-2017 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
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
C
C     * Neither the name of NTESS nor the names of its
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
C

c sdcgif - FORTRAN shell for handling strings


C  CESC - Escape
      SUBROUTINE CESC (FUNCID, LDR, DATA)
      INTEGER FUNCID
      INTEGER LDR
      CHARACTER*(*) DATA

      CALL CESC1( FUNCID, LDR, DATA, LEN(DATA) )
      RETURN
      END


C  CTX - Text
      SUBROUTINE CTX(X, Y, FLAG, TEXT )
      REAL X,Y
      INTEGER FLAG
      CHARACTER*(*) TEXT

      CALL CTX1( X, Y, FLAG, TEXT, LEN(TEXT) )
      RETURN
      END


C  CGTXX - Get Text Extent
      SUBROUTINE CGTXX( X, Y, STRING, VSTAT, VCONC, XCONC, YCONC,
     1                  X1, Y1, X2, Y2, X3, Y3, X4, Y4)
      REAL X,Y
      CHARACTER*(*) STRING
      INTEGER VSTAT, VCONC
      REAL XCONC, YCONC
      REAL X1, Y1, X2, Y2, X3, Y3, X4, Y4

      CALL CGTXX1( X, Y, STRING, VSTAT, VCONC, XCONC, YCONC,
     1             X1, Y1, X2, Y2, X3, Y3, X4, Y4, LEN(STRING) )
      RETURN
      END


C  CQCHH - Inquire List of Available Character Heights
      SUBROUTINE CQCHH( FONT, TXP, NREQ, FIRST, VSTAT, NTOTAL,
     1                  NLIST, CHHIT )
      CHARACTER*(*)FONT
      INTEGER TXP, NREQ, FIRST, VSTAT, NTOTAL, NLIST
      INTEGER CHHIT(*)

      CALL CQCHH1( FONT, TXP, NREQ, FIRST, VSTAT, NTOTAL,
     1             NLIST, CHHIT, LEN(FONT) )
      RETURN
      END
