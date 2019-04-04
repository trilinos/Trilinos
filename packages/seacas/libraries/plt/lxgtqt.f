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

C $Id: lxgtqt.f,v 1.1 1993/07/16 16:46:41 gdsjaar Exp $
C $Log: lxgtqt.f,v $
C Revision 1.1  1993/07/16 16:46:41  gdsjaar
C Changed plt to library rather than single source file.
C
C=======================================================================
      LOGICAL FUNCTION LXGTQT(STR,NS,CH)
      IMPLICIT INTEGER (A-Z)
      CHARACTER*504 ILINE
      COMMON /LXCOM1/ILINE
      COMMON /LXCOM2/JLINE,LXINIT
      CHARACTER*(*) STR
      CHARACTER CH

      CH = ILINE(JLINE:JLINE)
      LXGTQT = (CH.EQ.CHAR(39)) .OR. (CH.EQ.CHAR(34))
      IF (.NOT.LXGTQT) THEN
         RETURN

      END IF

      NS = 0
      J = JLINE
      LP = 1
 2470 IF (.NOT. (LP.GT.0)) GO TO 2480
      J = J + 1
      CH = ILINE(J:J)
      IF (CH.EQ.CHAR(39)) THEN
         LP = LP - 1

      ELSE IF (CH.EQ.CHAR(34)) THEN
         LP = LP - 1

      ELSE IF (CH.EQ.CHAR(0)) THEN
         LXGTQT = .FALSE.
         RETURN

      END IF

      NS = NS + 1
      STR(NS:NS) = CH
      GO TO 2470

 2480 CONTINUE
      NS = NS - 1
      JLINE = J + 1
      RETURN

      END
