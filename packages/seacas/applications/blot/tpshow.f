C Copyright(C) 2009-2017 National Technology & Engineering Solutions of
C Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
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

C=======================================================================
      SUBROUTINE TPSHOW (SHOTYP, NAMES, MAPEL, MAPND)
C=======================================================================

C   --*** TPSHOW *** (TPLOT) Display TPLOT parameter information
C   --   Written by Amy Gilkey - revised 12/17/87
C   --
C   --TPSHOW displays the TPLOT plot parameters.
C   --
C   --The SHOW options with the items they display are:
C   --   TYPLOT   - the curves to be plotted for the plot set
C   --   XYPLOT   -
C   --   PLOT     -
C   --   HARDCOPY -
C   --   NEUTRAL  -
C   --
C   --Parameters:
C   --   SHOTYP - IN - the expanded SHOW option string
C   --   NAMES - IN - the variable names
C   --
C   --Common Variables:
C   --   Uses NTPVAR, TIMPLT, ITVID, ITVNE of /TPVARS/

      include 'params.blk'
      include 'tpvars.blk'

      CHARACTER*(*) SHOTYP
      CHARACTER*(*) NAMES(*)
      INTEGER MAPEL(*), MAPND(*)

      LOGICAL ISABRT
      CHARACTER*(1024) PV1, PV2
      CHARACTER*2 STR2

      IF ((SHOTYP .EQ. 'TYPLOT') .OR. (SHOTYP .EQ. 'XYPLOT')
     &   .OR. (SHOTYP .EQ. 'PLOT') .OR. (SHOTYP .EQ. 'HARDCOPY')
     &   .OR. (SHOTYP .EQ. 'NEUTRAL')) THEN
         N = 1
         DO 100 NP = 1, NTPCRV
            IF (ISABRT ()) RETURN
            IF (TIMPLT) THEN
               PV1 = 'TIME'
            ELSE
               CALL TPLABV (-1, ITVID(N), NAMES(ITVID(N)), ITVNE(N),
     &            PV1, MAPEL, MAPND)
               N = N + 1
            END IF
            CALL TPLABV (-1, ITVID(N), NAMES(ITVID(N)), ITVNE(N), PV2,
     *        MAPEL, MAPND)
            N = N + 1
            WRITE (STR2, '(I2)', IOSTAT=IDUM) NP
            WRITE (*, 10000) 'Curve ', STR2, ' :  ', PV2(:LENSTR(PV2)),
     &         '  -vs-  ', PV1(:LENSTR(PV1))
  100    CONTINUE

      END IF

      RETURN

10000  FORMAT (1X, 10A)
      END
