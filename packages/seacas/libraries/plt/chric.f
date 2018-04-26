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

C $Id: chric.f,v 1.1 1993/07/16 16:46:19 gdsjaar Exp $ 
C $Log: chric.f,v $
C Revision 1.1  1993/07/16 16:46:19  gdsjaar
C Changed plt to library rather than single source file.
C 
C=======================================================================
      SUBROUTINE CHRIC(JIN,DS,ND)
      CHARACTER*(*) DS
      CHARACTER SDG*16

      DS = ' '
      J = ABS(JIN)
      ND = 0
 2030 CONTINUE
      ND = ND + 1
      KD = MOD(J,10)
      J = J/10
      SDG(ND:ND) = CHAR(ICHAR('0')+KD)

      IF (.NOT. (J.EQ.0)) GO TO 2030
      IJ = 1
      IF (JIN.LT.0) THEN
         IJ = 2
         DS(1:1) = '-'
         ND = ND + 1
      END IF

      I = IJ

 2060 IF (.NOT. (I.LE.ND)) GO TO 2080
      DS(I:I) = SDG(ND-I+1:ND-I+1)
      I = I + 1
      GO TO 2060

 2080 CONTINUE
      RETURN

      END
