C Copyright(C) 2011-2017 National Technology & Engineering Solutions of
C Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C Redistribution and use in source and binary forms, with or without
C modification, are permitted provided that the following conditions are
C met:
C
C * Redistributions of source code must retain the above copyright
C    notice, this list of conditions and the following disclaimer.
C
C * Redistributions in binary form must reproduce the above
C   copyright notice, this list of conditions and the following
C   disclaimer in the documentation and/or other materials provided
C   with the distribution.
C
C * Neither the name of NTESS nor the names of its
C   contributors may be used to endorse or promote products derived
C   from this software without specific prior written permission.
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
      SUBROUTINE INCID (TYPE, IDLST, NUMID, IDINC)
C=======================================================================
C
      CHARACTER*(*) TYPE
      INTEGER       IDLST(*)
      CHARACTER*80  STRING
      CHARACTER*8   STRA

C ... This routine increments all IDS by IDINC.

      IF (TYPE(:1) .EQ. 'M') THEN
         STRA = 'Material'
      ELSE IF (TYPE(:1) .EQ. 'S') THEN
         STRA = 'Sideset'
      ELSE IF (TYPE(:1) .EQ. 'N') THEN
         STRA = 'Nodeset'
      ELSE
         CALL PRTERR ('PROGRAM', 'unrecognized id type in INCID')
         RETURN
      END IF

C ... Just increment all ids by IDINC.  Note that we assume IDs valid
C     at entrance to routine, therefore, list will be valid at exit.

      IF (NUMID .LE. 0) RETURN
      DO 20 ID = 1, NUMID
         WRITE (STRING, 10) STRA, IDLST(ID), STRA, IDLST(ID)+IDINC
   10    FORMAT (A,1X,I10,' changed to ',A,1X,I10)
         CALL SQZSTR (STRING, LSTR)
         CALL PRTERR ('CMDSPEC', STRING(:LSTR))
         IDLST(ID) = IDLST(ID) + IDINC
   20 CONTINUE

      RETURN
      END
