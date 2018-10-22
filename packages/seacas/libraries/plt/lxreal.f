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

C $Id: lxreal.f,v 1.2 1993/07/16 18:28:51 gdsjaar Exp $ 
C $Log: lxreal.f,v $
C Revision 1.2  1993/07/16 18:28:51  gdsjaar
C Changed real*8 to double precision
C
c Revision 1.1  1993/07/16  16:46:45  gdsjaar
c Changed plt to library rather than single source file.
c 
C=======================================================================
      LOGICAL FUNCTION LXREAL(VALUE,CH)
      LOGICAL LXSET,LXNUMB
      LOGICAL PREF,POSF
      DOUBLE PRECISION VALUE,PREDIG,ESIGN
      DOUBLE PRECISION FN
      CHARACTER CH
      INTEGER ND

      LXREAL = .FALSE.
      ESIGN = 1.
      PREDIG = 0.
      ISIGN = 1
      ISAVE = LXSV()
      IF (LXSET('+-',CH)) THEN
         IF (CH.EQ.'-') THEN
            ISIGN = -1
         END IF

      END IF

      PREF = LXNUMB(PREDIG,ND,CH)
      POSF = .FALSE.
      IF (LXSET('.',CH)) THEN
         IF (LXNUMB(FN,ND,CH)) THEN
            POSF = .TRUE.
            PREDIG = PREDIG + FN*10.** (FLOAT(-ND))
         END IF

      END IF

      PREDIG = PREDIG*ISIGN
      IF (.NOT. (PREF.OR.POSF)) THEN
         CALL LXRS(ISAVE)
         RETURN

      END IF

      IF (LXSET('EeDdQq',CH)) THEN
         IF (LXSET('+-',CH)) THEN
            IF (CH.EQ.'-') THEN
               ESIGN = -1.
            END IF

         END IF

         IF (LXNUMB(FN,ND,CH)) THEN
            PREDIG = PREDIG*10.** (ESIGN*FN)

         ELSE
            CALL LXRS(ISAVE)
            RETURN

         END IF

      END IF

      VALUE = PREDIG
      LXREAL = .TRUE.
      RETURN

      END
