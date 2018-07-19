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

C $Id: memfre.f,v 1.1 1993/07/16 16:46:59 gdsjaar Exp $ 
C $Log: memfre.f,v $
C Revision 1.1  1993/07/16 16:46:59  gdsjaar
C Changed plt to library rather than single source file.
C 
C=======================================================================
      LOGICAL FUNCTION MEMFRE(IB,MEMRY)
      INTEGER IB
      INTEGER MEMRY(*)

      IPR = IB - 2
      MEMFRE = .FALSE.
      IP = 3
      IQ = 3
      IPN = MEMRY(IP)
 2710 IF (.NOT. (IPN.NE.0)) GO TO 2720
 2730 IF (.NOT. (IPN.LT.0)) GO TO 2740
      IPM = MEMRY(-IPN)
      IF (IPM.GT.0) THEN
         IQ = IP
         IP = -IPN
         IPN = IPM

      ELSE IF (IPM.EQ.0) THEN
         MEMRY(IP) = 0
         MEMRY(2) = IP
         RETURN

      ELSE IF (IPM.LT.0) THEN
         IPN = IPM
         MEMRY(IP) = IPN
      END IF

      GO TO 2730

 2740 CONTINUE
      IF (IP.EQ.IPR) THEN
         IPN = -IPN
         IF (MEMRY(IQ).LT.0) THEN
            IP = IQ
         END IF

         MEMRY(IP) = IPN
         MEMFRE = .TRUE.
         IB = 0

      ELSE
         IQ = IP
         IP = ABS(IPN)
         IPN = MEMRY(IP)
      END IF

      GO TO 2710

 2720 CONTINUE
      RETURN

      END
