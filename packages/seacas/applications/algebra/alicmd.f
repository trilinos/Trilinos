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
      SUBROUTINE ALICMD (INLINE, INTYP, CFIELD, IFIELD, NAMES, *)
C=======================================================================

C   --*** ALICMD *** (ALGEBRA) Perform ALIAS command
C   --   Written by Amy Gilkey - revised 12/02/87
C   --
C   --ALICMD processes the input ALIAS command.  It adds the alias and
C   --the equivalent variables to the /ALIAS../ arrays.
C   --
C   --Parameters:
C   --   INLINE - IN/OUT - the parsed input lines for the log file
C   --   INTYP - IN - the field types
C   --   CFIELD - IN - the character fields
C   --   IFIELD - IN - the integer fields
C   --   NAMES - IN - the names of the global, nodal, and element variables
C   --   * - return statement if command not executed
C   --
C   --Common Variables:
C   --   Sets NUMALI, NAMALI, NIXALI, IXALI of /ALIAS../
C   --   Uses NVARHI, NVARGL, NVARNP, NVAREL of /DBNUMS/


      include 'params.blk'
      include 'namlen.blk'

      include 'alias.blk'
      include 'dbnums.blk'

      CHARACTER*(*) INLINE
      INTEGER INTYP(*)
      CHARACTER*(*) CFIELD(*)
      INTEGER IFIELD(*)
      CHARACTER*(namlen) NAMES(*)

      LOGICAL FFEXST, FFNUMB
      CHARACTER*(maxnam) NAME
      CHARACTER TYPNAM, TYP

c      DATA NUMALI / 0 /
      NUMALI = 0 

      IF (.NOT. FFEXST (1, INTYP)) THEN
         CALL PRTERR ('CMDERR', 'No options on ALIAS command')
         GOTO 130
      END IF

      IF (NUMALI .GE. MAXALI) THEN
         CALL PRTERR ('CMDERR', 'ALIAS table full, command ignored')
         GOTO 130
      END IF

      NVNAMS = NVARHI + NVARGL + NVARNP + NVAREL

C   --Get the tensor name

      IFLD = 1
      CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', NAME)

      IVAR = LOCSTR (NAME, NVNAMS, NAMES)
      IF (IVAR .GT. 0) THEN
         CALL PRTERR ('CMDERR', 'ALIAS name must be unique')
         GOTO 130
      END IF

      IVAR = LOCSTR (NAME, NUMALI, NAMALI)
      IF (IVAR .GT. 0) THEN
         CALL PRTERR ('CMDERR', 'ALIAS name must be unique')
         GOTO 130
      END IF

      IALI = NUMALI + 1
      NAMALI(IALI) = NAME
      CALL FFADDC (NAME, INLINE)

C   --Get and find the first variable name

      CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', NAME)
      CALL FFADDC (NAME, INLINE)

      IVAR = LOCSTR (NAME, NVNAMS, NAMES)
      IF (IVAR .GT. 0) THEN
         NIXALI(IALI) = 1
         IXALI(1,IALI) = IVAR
         CALL DBVTYP (IVAR, TYPNAM, IDUM)
      ELSE
         CALL PRTERR ('CMDERR',
     &      'Invalid variable name "' // NAME(:LENSTR(NAME)) // '"')
         GOTO 130
      END IF

C   --If a number follows, assign the next n variables

      IF (FFNUMB (IFLD, INTYP)) THEN
         CALL FFINTG (IFLD, INTYP, IFIELD,
     &      'number of variables', 1, NIXALI(IALI), *130)
         CALL FFADDI (NIXALI(IALI), INLINE)

         IF (NIXALI(IALI) .GT. MAXALN) THEN
            CALL PRTERR ('CMDERR', 'Too many variables to alias')
            GOTO 130
         END IF

         I = IVAR + NIXALI(IALI) - 1
         CALL DBVTYP (I, TYP, IDUM)
         IF (TYPNAM .NE. TYP) THEN
            CALL PRTERR ('CMDERR', 'Variables do not exist')
            GOTO 130
         END IF

         DO 100 I = 1, NIXALI(IALI)
            IXALI(I,IALI) = IVAR + I - 1
  100    CONTINUE

C   --Otherwise assign the named variables

      ELSE

  110    CONTINUE
         IF (FFEXST (IFLD, INTYP)) THEN

            CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', NAME)

            IVAR = LOCSTR (NAME, NVNAMS, NAMES)
            IF (IVAR .GT. 0) THEN
               NIXALI(IALI) = NIXALI(IALI) + 1
               IF (NIXALI(IALI) .GT. MAXALN) THEN
                  CALL PRTERR ('CMDERR', 'Too many variables to alias')
                  GOTO 130
               END IF

               IXALI(NIXALI(IALI),IALI) = IVAR

               CALL DBVTYP (IVAR, TYP, IDUM)
               IF (TYPNAM .NE. TYP) THEN
                  CALL PRTERR ('CMDERR',
     &               'Variables must be of the same type')
                  GOTO 120
               END IF

            ELSE
               CALL PRTERR ('CMDERR',
     &            'Invalid variable name' // NAME(:LENSTR(NAME)))
               GOTO 120
            END IF

            CALL FFADDC (NAME, INLINE)
  120       CONTINUE
            GOTO 110
         END IF
      END IF

      NUMALI = IALI

      RETURN

  130 CONTINUE
      RETURN 1
      END
