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

C $Id: lxscnp.f,v 1.1 1993/07/16 16:46:50 gdsjaar Exp $
C $Log: lxscnp.f,v $
C Revision 1.1  1993/07/16 16:46:50  gdsjaar
C Changed plt to library rather than single source file.
C
C=======================================================================
      LOGICAL FUNCTION LXSCNP(DELIM,STR,NS,CH)
      IMPLICIT INTEGER (A-Z)
      CHARACTER*504 ILINE
      COMMON /LXCOM1/ILINE
      COMMON /LXCOM2/JLINE,LXINIT
      CHARACTER*(*) DELIM
      CHARACTER*(*) STR
      INTEGER NS
      CHARACTER CH
      LOGICAL QFLAG
      CHARACTER QCH

      NS = 0
      J = JLINE
      LP = 0
      QFLAG = .FALSE.
 2600 CONTINUE
      CH = ILINE(J:J)
      IF (LP.EQ.0) THEN
         ID = INDEX(DELIM,CH)
         IF (ID.GT.0) THEN
            GO TO 2620

         END IF

      END IF

      IF (.NOT.QFLAG) THEN
         IF (LP.EQ.0) THEN
            ID = INDEX(DELIM,CH)
            IF (ID.GT.0) THEN
               GO TO 2620

            END IF

         END IF

         IF (CH.EQ.CHAR(0)) THEN
            GO TO 2620

         ELSE IF (CH.EQ.'(' .OR. CH.EQ.'[' .OR. CH.EQ.'{') THEN
            LP = LP + 1

         ELSE IF (CH.EQ.')' .OR. CH.EQ.']' .OR. CH.EQ.'}') THEN
            LP = LP - 1

         ELSE IF (CH.EQ.'''' .OR. CH.EQ.'"' .OR. CH.EQ.CHAR(96)) THEN
            QFLAG = .TRUE.
            QCH = CH
         END IF

         IF (LP.LT.0) THEN
            GO TO 2620

         END IF

         NS = NS + 1
         STR(NS:NS) = CH

      ELSE
         NS = NS + 1
         STR(NS:NS) = CH
         IF (CH.EQ.CHAR(0)) THEN
            GO TO 2620

         ELSE IF (CH.EQ.QCH) THEN
            IF (CH.EQ.ILINE(J+1:J+1)) THEN
               NS = NS + 1
               STR(NS:NS) = CH
               J = J + 2
               GO TO 2610

            ELSE
               QFLAG = .FALSE.
            END IF

         END IF

      END IF

      J = J + 1
 2610 GO TO 2600

 2620 CONTINUE
      JLINE = J
      LXSCNP = ((LP.EQ.0) .AND. .NOT.QFLAG)
      RETURN

      END
