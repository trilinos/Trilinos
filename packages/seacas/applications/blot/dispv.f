C Copyright(C) 2009 Sandia Corporation. Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
C certain rights in this software.
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
C     * Neither the name of Sandia Corporation nor the names of its
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

C============================================================================
      SUBROUTINE DISPV (INIT, INLINE, IFLD, INTYP, CFIELD,
     &   NAMES, LIDSP, NAMLEN)
C============================================================================

C   --*** DISPV ***  (BLOT) Process DISPVAR command
C   --   Written by John Glick - 11/1/88
C   --
C   --Parameters:
C   --   INIT  - IN - .TRUE. iff the purpose of the call is to
C   --            initialize the list of display variables and
C   --            not to modify it.
C   --   INLINE - IN/OUT - the parsed input lines for the log file
C   --   IFLD, INTYP, CFIELD, - IN/OUT - the free-field reader index
C   --          and charcter field.
C   --   NAMES - IN - the variable names
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
C
      include 'params.blk'
      LOGICAL INIT
      CHARACTER*(*) INLINE(*)
      INTEGER IFLD, INTYP(*)
      CHARACTER*(*) CFIELD(*)
      CHARACTER*(NAMLEN) NAMES(*)
      INTEGER LIDSP(0:*)

      include 'dbnums.blk'
      include 'dbnumgq.blk'

      LOGICAL FFMATC

      CHARACTER*(MXNAME) TNAME(4096)
      SAVE TNAME
      INTEGER LLIST
      SAVE LLIST
      LOGICAL FIRST
      SAVE FIRST

      DATA FIRST /.TRUE./

C *****************************************************************

C        Check that database is in the EXODUS format

      CALL CKEXOD (EXODUS, *140)

      IF (FIRST) THEN

C        Construct array containing names of history and global
C        variables and TIME.

         if (nvarhi+nvargl+1 .gt. 4096) then
            call prterr('PROGRAM',
     $           'Too many names; increase TNAME size')
         end if
         TNAME(1) = 'TIME'
         CALL DBVIX ('H', 1, IXHV)
         CALL CPYSTR (NVARHI, NAMES(IXHV), TNAME(2))
         CALL DBVIX ('G', 1, IXGL)
         CALL CPYSTR (NVARGL, NAMES(IXGL), TNAME(NVARHI+2))
         LLIST = 1 + NVARHI + NVARGL
         FIRST = .FALSE.
      ENDIF

      IF (INIT) THEN

         LIDSP(0) = 1
         LIDSP(1) = 0

      ELSE

C              If there are no fields on the DISPV command line,
C              toggle the ON/OFF flag.

         IF (INTYP(IFLD) .LE. -1) THEN
            LIDSP(0) = -LIDSP(0)

C              Check for ON flag.

         ELSE IF (FFMATC (IFLD, INTYP, CFIELD, 'ON', 2)) THEN
            LIDSP(0) = ABS(LIDSP(0))
            CALL FFADDC ('ON', INLINE(1))

C              Check for OFF flag.

         ELSE IF (FFMATC (IFLD, INTYP, CFIELD, 'OFF', 3)) THEN
            LIDSP(0) = - ABS(LIDSP(0))
            CALL FFADDC ('OFF', INLINE(1))

C              Check for ALL parameter

         ELSE IF (FFMATC (IFLD, INTYP, CFIELD, 'ALL', 3)) THEN
            LIDSP(0) = LLIST
            LIDSP(1) = 0
            DO 100 I = 1, NVARHI
               LIDSP(I+1) = I
  100       CONTINUE
            DO 110 I = 1, NVARGL
               LIDSP(I+NVARHI+1) = -I
  110       CONTINUE
            CALL FFADDC ('ALL', INLINE(1))

         ELSE

C              Temporarily convert existing list to a list
C              containing all positive integers.

            LIDSP(0) = ABS(LIDSP(0))
            DO 120 I = 1, LIDSP(0)
               IF (LIDSP(I) .EQ. 0) THEN
                  LIDSP(I) = 1
               ELSE IF (LIDSP(I) .LT. 0) THEN
                  LIDSP(I) = 1 + NVARHI - LIDSP(I)
               ENDIF
  120       CONTINUE

C              Get list of selected variables

            CALL RIXWRD (INLINE(1), IFLD, INTYP, CFIELD,
     &         'display variable name', LLIST, TNAME,
     &         LIDSP(0), LIDSP(1), *140)

C              Convert the list back to the proper format.

            DO 130 I = 1, LIDSP(0)
               IF (LIDSP(I) .EQ. 1) THEN
                  LIDSP(I) = 0
               ELSE IF (LIDSP(I) .LE. (NVARHI+1)) THEN
                  LIDSP(I) = LIDSP(I) - 1
               ELSE
                  LIDSP(I) = -(LIDSP(I) - NVARHI - 1)
               ENDIF
  130       CONTINUE

         ENDIF

      ENDIF
      GO TO 150

  140 CONTINUE
      INLINE(1) = ' '

  150 CONTINUE
      RETURN
      END
