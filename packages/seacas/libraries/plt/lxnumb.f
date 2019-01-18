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

C $Id: lxnumb.f,v 1.2 1993/07/16 18:28:50 gdsjaar Exp $ 
C $Log: lxnumb.f,v $
C Revision 1.2  1993/07/16 18:28:50  gdsjaar
C Changed real*8 to double precision
C
c Revision 1.1  1993/07/16  16:46:44  gdsjaar
c Changed plt to library rather than single source file.
c 
C=======================================================================
      LOGICAL FUNCTION LXNUMB(N,ND,CH)
      CHARACTER*(*) CH
      DOUBLE PRECISION N
      LOGICAL LXSET

      N = 0.
      ND = 0
 2580 IF (.NOT. (LXSET('0123456789',CH))) GO TO 2590
      N = N*10. + FLOAT(ICHAR(CH)-ICHAR('0'))
      ND = ND + 1
      GO TO 2580

 2590 CONTINUE
      IF (ND.EQ.0) THEN
         LXNUMB = .FALSE.
         RETURN

      END IF

      LXNUMB = .TRUE.
      RETURN

      END
