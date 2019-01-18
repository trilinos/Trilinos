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

C $Id: chrrvc.f,v 1.1 1993/07/16 16:46:21 gdsjaar Exp $ 
C $Log: chrrvc.f,v $
C Revision 1.1  1993/07/16 16:46:21  gdsjaar
C Changed plt to library rather than single source file.
C 
C=======================================================================
      SUBROUTINE CHRRVC(BUFF,TXT,L)
      CHARACTER*16 LOCTXT
      CHARACTER*(*) TXT

      LT = LEN(TXT)
      WRITE (LOCTXT,'(1pg13.6)') BUFF
      DO 2110 I = 1,LT
         IF (LOCTXT(I:I).NE.' ') THEN
            L1 = I
            GO TO 2120

         END IF

 2110 CONTINUE
 2120 CONTINUE
      DO 2130 I = L1,LT
         IF (LOCTXT(I:I).EQ.' ') THEN
            L2 = I
            GO TO 2140

         END IF

 2130 CONTINUE
 2140 CONTINUE
      TXT = LOCTXT(L1:L2)
      L = L2 - L1
      RETURN

      END
