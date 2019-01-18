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
      SUBROUTINE SPSHOW (SHOTYP, NAMES, NENUM, LIDSP)
C=======================================================================

C   --*** SPSHOW *** (SPLOT) Display SPLOT parameter information
C   --   Modified by John Glick - 11/9/88
C   --   Written by Amy Gilkey - revised 12/17/87
C   --
C   --SPSHOW displays the SPLOT plot parameters.
C   --
C   --The SHOW options with the items they display are:
C   --   SYPLOT   - the curves to be plotted for the plot set
C   --   PLOT     -
C   --   HARDCOPY -
C   --   NEUTRAL  -
C   --   NODES    - the number of selected nodes or elements, and the
C   --                 selected node/element numbers (if lower-case)
C   --   ELEMENTS -
C   --   DISPVAR  - select history variables, global variables, and/or
C   --              TIME whose values will be displayed on the plot legend.
C   --
C   --Parameters:
C   --   SHOTYP - IN - the expanded SHOW option string
C   --   NAMES - IN - the variable names
C   --   NENUM - IN - the selected node/element numbers
C   --   LIDSP(0:*)  - IN/OUT - the indices of the selected variables
C   --          whose values will be displayed on the plot legend.
C   --          ABS(LIDSP(0)) = the number of variables in the list.
C   --          SIGN(LIDSP(0)) specifies whether the variables in the
C   --                   list should have their values displayed on
C   --                   the plot legend.  If >0, they should;
C   --                   If <=0, they should not.
C   --          LIDSP(i) identifies the ith variable in the list.
C   --          If LIDSP(i) > 0, LIDSP(i) is the id of a history variable.
C   --          If LIDSP(i) < 0, -LIDSP(i) is the id of a global variable.
C   --          If LIDSP(i) = 0, TIME is to be displayed on the plot legend.
C   --
C   --Common Variables:
C   --   Uses NODVAR, NNENUM of /SELNE/
C   --   Uses NSPVAR, ISVID of /SPVARS/

      include 'params.blk'
      include 'selne.blk'
      include 'spvars.blk'

      CHARACTER*(*) SHOTYP
      CHARACTER*(*) NAMES(*)
      INTEGER NENUM(*)
      INTEGER LIDSP(0:*)

      LOGICAL ISABRT
      CHARACTER*2 STR2
      CHARACTER*5 STR5
      CHARACTER*(MXNAME) NAM
      CHARACTER*(MXNAME) TNAME(5)

      IF ((SHOTYP .EQ. 'SYPLOT')
     &   .OR. (SHOTYP .EQ. 'PLOT') .OR. (SHOTYP .EQ. 'HARDCOPY')
     &   .OR. (SHOTYP .EQ. 'NEUTRAL')) THEN
         DO 100 N = 1, NSPVAR
            IF (ISABRT ()) RETURN
            WRITE (STR2, '(I2)', IOSTAT=IDUM) N
            NAM = NAMES(ISVID(N))
            WRITE (*, 10010) 'Curve ', STR2, ' :  ',
     &          NAM(:LENSTR(NAM)), ' -vs-  DISTANCE'
  100    CONTINUE

      ELSE IF ((SHOTYP .EQ. 'NODES') .OR. (SHOTYP .EQ. 'ELEMENTS')
     &   .OR. (SHOTYP .EQ. 'nodes') .OR. (SHOTYP .EQ. 'elements')) THEN
         CALL INTSTR (1, 0, NNENUM, STR5, LSTR)
         IF (NODVAR) THEN
            WRITE (*, 10010)
     &         'Number of selected node numbers = ', STR5(:LSTR)
         ELSE
            WRITE (*, 10010)
     &         'Number of selected element numbers = ', STR5(:LSTR)
         END IF

         IF ((SHOTYP .EQ. 'NODES') .OR. (SHOTYP .EQ. 'ELEMENTS')) THEN
            WRITE (*, 10020) (NENUM(I), I=1,NNENUM)
         END IF

      ELSE IF (SHOTYP .EQ. 'DISPVAR') THEN
         IF (LIDSP(0) .GT. 0) THEN
            WRITE (*, 10010) 'Values of display variables will',
     &         ' appear on plot legend'
         ELSE
            WRITE (*, 10010) 'Values of display variables will not',
     &         ' appear on plot legend'
         ENDIF
         WRITE (*, 10010)
         WRITE (*, 10010) 'Display variables - '
         NUMDSP = ABS(LIDSP(0))
         NUM = 0
         NLAST = 0
  110    CONTINUE
         IF (NUM .LT. NUMDSP) THEN
            NFIRST = NLAST + 1
            NLAST = NFIRST + 4
            IF (NLAST .GT. NUMDSP) NLAST = NUMDSP
            N = NLAST - NFIRST + 1
            NID = NFIRST
            DO 120 I = 1, N
               IF (LIDSP(NID) .EQ. 0) THEN
                  TNAME(I) = 'TIME'
               ELSE IF (LIDSP(NID) .GT. 0) THEN
                  CALL DBVIX_BL ('H', LIDSP(NID), IDVAR)
                  TNAME(I) = NAMES(IDVAR)
               ELSE IF (LIDSP(NID) .LT. 0) THEN
                  CALL DBVIX_BL ('G', -LIDSP(NID), IDVAR)
                  TNAME(I) = NAMES(IDVAR)
               ENDIF
               NID = NID + 1
  120       CONTINUE
            WRITE (*, 10000)(TNAME(I), I=1,N)
10000        FORMAT (5X, 5(A32,2X))
            NUM = NUM + N
            GO TO 110
         ENDIF

      END IF

      RETURN

10010  FORMAT (1X, 10A)
10020  FORMAT (1X, 12I9)
      END
