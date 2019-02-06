C    Copyright(C) 2008-2017 National Technology & Engineering Solutions of
C    Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
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

C=======================================================================
      SUBROUTINE DELCMD (INLINE, INTYP, CFIELD, *)
C=======================================================================

C   --*** DELCMD *** (ALGEBRA) Perform DELETE command
C   --   Written by Amy Gilkey - revised 11/23/87
C   --
C   --DELCMD processes the input DELETE command.  It sets the ISTVAR
C   --flag in the /VAR../ arrays for all deleted variables.
C   --
C   --Parameters:
C   --   INLINE - IN/OUT - the parsed input lines for the log file
C   --   INTYP - IN - the field types
C   --   CFIELD - IN - the character fields
C   --   * - return statement if command not executed
C   --
C   --Common Variables:
C   --   Uses IXLHS, NAMVAR of /VAR../
C   --   Sets ISTVAR of /VAR../

      PARAMETER (ICURTM = 1, ILSTTM = 2, IONETM = 3)
      include 'ag_namlen.blk'
      include 'ag_var.blk'
C     database type, num_of qa and info records
      include 'ag_dbnumq.blk'

      CHARACTER*(*) INLINE
      INTEGER INTYP(*)
      CHARACTER*(*) CFIELD(*)

      LOGICAL FFEXST, MATSTR
      CHARACTER*(maxnam) NAME

      NUMLHS = MAXVAR - IXLHS + 1

      IFLD = 1
  100 CONTINUE
      IF (FFEXST (IFLD, INTYP)) THEN

         CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', NAME)

         IF (MATSTR (NAME, 'QAINFO', 6)) THEN
           NQAREC = 0
           NINFO  = 0
         ELSE IF (MATSTR (NAME, 'QA', 2)) THEN
           NQAREC = 0
         ELSE IF (MATSTR (NAME, 'INFORMATION', 11)) THEN
           NINFO  = 0
         ELSE
           IVAR = LOCSTR (NAME, NUMLHS, NAMVAR(IXLHS))
           IF (IVAR .EQ. 0) THEN
             CALL PRTERR ('CMDERR',
     &         '"' // NAME(:LENSTR(NAME)) // '" not defined, ignored')
             GOTO 110
           END IF
           ISTVAR(ICURTM,IVAR+IXLHS-1) = -1
         END IF
         CALL FFADDC (NAME, INLINE)


  110    CONTINUE
         GOTO 100
      END IF

      RETURN
      END
