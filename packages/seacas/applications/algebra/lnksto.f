C    Copyright(C) 2008 Sandia Corporation.  Under the terms of Contract
C    DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
C    certain rights in this software
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
C    * Neither the name of Sandia Corporation nor the names of its
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
      SUBROUTINE LNKSTO (NAMEGV, NUMSTO, LTMENT, MERR)
C=======================================================================

C   --*** LNKSTO *** (ALGEBRA) Assign storage for variables
C   --   Written by Amy Gilkey - revised 11/30/87
C   --
C   --LNKSTO sets up the storage locations for all the variables.
C   --Input and output variables of the same name share storage (unless
C   --one is a history/global and the other is not).  Time and history/global
C   --variables are all in the first storage location: time is in slot 1,
C   --followed by the input history variables (if any), the input global
C   --variables (if any), then the output only history/global variables.
C   --IDVAR contains the slot index of the time or history/global variables.
C   --
C   --Parameters:
C   --   NAMEGV - IN - the global variable names
C   --   NUMSTO - OUT - the number of variable storage locations needed
C   --   LTMENT - OUT - the number of slots in the time /globals entry
C   --   * - return statement if an error is found; message is printed
C   --
C   --Common Variables:
C   --   Sets ISTVAR, IDVAR of /VAR../
C   --   Uses NUMINP, IXLHS, NAMVAR, TYPVAR of /VAR../
C   --   Uses NVARGL of /DBNUMS/
C   --   Uses IHVBEG, IHVEND, IGVBEG, IGVEND of /DBXVAR/

      PARAMETER (ICURTM = 1, ILSTTM = 2, IONETM = 3)
      include 'namlen.blk'
      include 'var.blk'
      include 'dbnums.blk'
      include 'dbxvar.blk'

      INTEGER MERR
      CHARACTER*(namlen) NAMEGV(*)

      INTEGER ISTTHG(3)

      MERR = 0

ccSave entry 1 for time and history/global variables
C     Save entry 1 for time and global variables
      NUMSTO = 1

C   --Assign current, last and first location for all time global
C   --variables

      ISTTHG(ICURTM) = NUMSTO
      ISTTHG(ILSTTM) = 0
      ISTTHG(IONETM) = 0

      DO 110 IVAR = 1, NUMINP
         IF ((TYPVAR(IVAR) .EQ. 'T')
     &      .OR. (TYPVAR(IVAR) .EQ. 'G')) THEN
            DO 100 ITM = 1, 3
               IF (ISTVAR(ITM,IVAR) .NE. 0) THEN
                  IF (ISTTHG(ITM) .EQ. 0) THEN
                     NUMSTO = NUMSTO + 1
                     ISTTHG(ITM) = NUMSTO
                  END IF
               END IF
  100       CONTINUE
         END IF
  110 CONTINUE

      DO 130 IVAR = IXLHS, MAXVAR
         IF ((TYPVAR(IVAR) .EQ. 'T')
     &      .OR. (TYPVAR(IVAR) .EQ. 'G')) THEN
            DO 120 ITM = 1, 3
               IF (ISTVAR(ITM,IVAR) .NE. 0) THEN
                  IF (ISTTHG(ITM) .EQ. 0) THEN
                     NUMSTO = NUMSTO + 1
                     ISTTHG(ITM) = NUMSTO
                  END IF
               END IF
  120       CONTINUE
         END IF
  130 CONTINUE

C   --Check if any history or global variables must be input, and assign
C   --index in time entry

c      IXHV = 0
      IXGV = 0
      DO 140 IVAR = 1, NUMINP
         IF (TYPVAR(IVAR) .EQ. 'G') THEN
            IXGV = 1
         END IF
  140 CONTINUE

      LTMENT = 1
      IF (IXGV .GT. 0) THEN
         IXGV = LTMENT + 1
         LTMENT = LTMENT + NVARGL
      END IF

C   --Reserve storage for all input variables

      DO 160 IVAR = 1, NUMINP
C      --If previous time is needed, current time must be read
         IF (ISTVAR(ILSTTM,IVAR) .NE. 0) ISTVAR(ICURTM,IVAR) = 1

         IF ((TYPVAR(IVAR) .EQ. 'T')
     &      .OR. (TYPVAR(IVAR) .EQ. 'G')) THEN

C         --Assign current and last for all time /global variables so
C         --current to last move is done properly

            IF (ISTVAR(ICURTM,IVAR) .NE. 0) THEN
               ISTVAR(ICURTM,IVAR) = ISTTHG(ICURTM)
               IF (ISTTHG(ILSTTM) .GT. 0)
     &            ISTVAR(ILSTTM,IVAR) = ISTTHG(ILSTTM)
            END IF
            IF (ISTVAR(IONETM,IVAR) .NE. 0)
     &         ISTVAR(IONETM,IVAR) = ISTTHG(IONETM)

            IF (TYPVAR(IVAR) .EQ. 'T') THEN
C            --TIME is in slot 1 (output time shares storage)
               IDVAR(IVAR) = 1

            ELSE IF (TYPVAR(IVAR) .EQ. 'G') THEN

C            --Input globals start at assigned slot
               IDVAR(IVAR) = IXGV
            END IF

         ELSE

C         --Reserve entry for array variable for current, last, and
C         --first time (as needed)

            DO 150 ITM = 1, 3
               IF (ISTVAR(ITM,IVAR) .NE. 0) THEN
                  NUMSTO = NUMSTO + 1
                  ISTVAR(ITM,IVAR) = NUMSTO
               END IF
  150       CONTINUE
         END IF
  160 CONTINUE

C   --Reserve storage for output variables or link to input variable

      DO 190 IVAR = IXLHS, MAXVAR

         IF ((TYPVAR(IVAR) .EQ. 'T')
     &      .OR. (TYPVAR(IVAR) .EQ. 'G')) THEN

C         --Assign current and last for all time / history/global variables
C         --so current to last move is done properly; output history/globals
C         --are in same entry as input history/globals

            IF (ISTVAR(ICURTM,IVAR) .NE. 0) THEN
               ISTVAR(ICURTM,IVAR) = ISTTHG(ICURTM)
               IF (ISTTHG(ILSTTM) .GT. 0)
     &            ISTVAR(ILSTTM,IVAR) = ISTTHG(ILSTTM)
            END IF
            IF (ISTVAR(IONETM,IVAR) .NE. 0)
     &         ISTVAR(IONETM,IVAR) = ISTTHG(IONETM)

            IF (TYPVAR(IVAR) .EQ. 'T') THEN

C            --Output time shares storage with input time
               IDVAR(IVAR) = 1

            ELSE IF (TYPVAR(IVAR) .EQ. 'G') THEN

C            --Share slot with input global variable of the same name (if any)
C            --or reserve new slot

               IF (IXGV .GT. 0) THEN
                  IINP = LOCSTR (NAMVAR(IVAR), NVARGL, NAMEGV)
               ELSE
                  IINP = 0
               END IF
               IF (IINP .GT. 0) THEN
                  IDVAR(IVAR) = IINP + IXGV - 1
               ELSE
                  LTMENT = LTMENT + 1
                  IDVAR(IVAR) = LTMENT
               END IF
            END IF

         ELSE

C         --Share array storage with input variable of the same name (if any)
C         --or reserve entry for array variable

            IINP = LOCSTR (NAMVAR(IVAR), NUMINP, NAMVAR)
            IF (IINP .GT. 0) THEN
               IF ((TYPVAR(IINP) .EQ. 'H')
     &            .OR. (TYPVAR(IINP) .EQ. 'G')) IINP = 0
            END IF

            IF (IINP .LE. 0) THEN

C            --Reserve entry for array variable for current, last, and
C            --first time (as needed)

               DO 170 ITM = 1, 3
                  IF (ISTVAR(ITM,IVAR) .NE. 0) THEN
                     NUMSTO = NUMSTO + 1
                     ISTVAR(ITM,IVAR) = NUMSTO
                  END IF
  170          CONTINUE

            ELSE

C            --Share array storage with input variable; if either has
C            --current/last pair, both must have current and last so
C            --current to last move is done properly

               DO 180 ITM = 1, 3
                  IF ((ISTVAR(ITM,IVAR) .NE. 0)
     &               .OR. (ISTVAR(ITM,IINP) .NE. 0)) THEN
                     IF (ISTVAR(ITM,IINP) .EQ. 0) THEN
                        NUMSTO = NUMSTO + 1
                        ISTVAR(ITM,IVAR) = NUMSTO
                        IF (ITM .NE. IONETM) ISTVAR(ITM,IINP) = NUMSTO
                     ELSE
                        ISTVAR(ITM,IVAR) = ISTVAR(ITM,IINP)
                     END IF
                  END IF
  180          CONTINUE
            END IF
         END IF
  190 CONTINUE

C   --Fix up variables which request the first time but not the current time
C   --by changing the current storage location to the first storage location
C   --negated

      DO 200 IVAR = 1, NUMINP
         IF ((ISTVAR(ICURTM,IVAR) .EQ. 0)
     &      .AND. (ISTVAR(IONETM,IVAR) .NE. 0)) THEN
            ISTVAR(ICURTM,IVAR) = - ISTVAR(IONETM,IVAR)
         END IF
  200 CONTINUE

      RETURN
      END
