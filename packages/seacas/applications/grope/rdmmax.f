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
      SUBROUTINE RDMMAX (IFLD, INTYP, CFIELD,
     &     NAMIGV, NAMINV, NAMIEV,
     &     NAMOGV, NAMONV, NAMOEV,
     &     NSTEP, MMSTEP, MMNAME, MMTYP, MMVAR, MMNUM, *)
C=======================================================================

C     --*** RDMMAX *** (GROPE) Parse the MINMAX command parameters
C     --
C     --RDMMAX parses the MINMAX command parameters.  It first reads
C     --the name, if any, and determines its type and variable number.
C     --It then reads the 'ALL' versus 'THIS' parameter.
C     --
C     --Parameters:
C     --   IFLD - IN/OUT - the number of the next entry to scan, incremented
C     --   INTYP - IN - the free-format field types
C     --      -1 = none, 0 = name, 1 = real, 2 = integer
C     --   CFIELD - IN - the input character fields
C     --   NAMIGV - IN - the names of the global variables as input
C     --   NAMINV - IN - the names of the nodal variables as input
C     --   NAMIEV - IN - the names of the element variables as input
C     --   NAMOGV - IN - the names of the global variables for comparison
C     --   NAMONV - IN - the names of the nodal variables for comparison
C     --   NAMOEV - IN - the names of the element variables for comparison
C     --   NSTEP - IN - the current step number
C     --   MMSTEP - IN/OUT - the requested step number, <=0 for all
C     --   MMNAME - IN/OUT - min/max variable name
C     --   MMTYP - IN/OUT - min/max variable type:
C     --      'G'lobal, 'N'odal, 'E'lement
C     --   MMVAR - IN/OUT - min/max variable number
C     --   MMNUM - IN/OUT - number of sequential min/max requests for this
C     --      variable
C     --   * - return statement if invalid range of integers specified;
C     --      message printed
C     --
C     --Common Variables:
C     --   Uses NVARNP, NVAREL, NVARGL of /DBNUMS/

      include 'params.blk'
      INCLUDE 'dbnums.blk'

      INTEGER INTYP(*)
      CHARACTER*(*) CFIELD(*)
      CHARACTER*(*) NAMIGV(*), NAMINV(*), NAMIEV(*)
      CHARACTER*(*) NAMOGV(*), NAMONV(*), NAMOEV(*)
      CHARACTER*(*) MMNAME
      CHARACTER MMTYP

      CHARACTER*(256) NAME, WORD
      CHARACTER CH

      CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', NAME)

      IF (NAME .NE. ' ') THEN

C     --Get the variable name

         IV = LOCSTR (NAME, NVARGL, NAMOGV)
         IF (IV .GT. 0) THEN
            CH = 'G'
            NAME = NAMIGV(IV)
         ELSE
            IV = LOCSTR (NAME, NVARNP, NAMONV)
            IF (IV .GT. 0) THEN
               CH = 'N'
               NAME = NAMINV(IV)
            ELSE
               IV = LOCSTR (NAME, NVAREL, NAMOEV)
               IF (IV .GT. 0) THEN
                  CH = 'E'
                  NAME = NAMIEV(IV)
               END IF
            END IF
         END IF
         IF (IV .LE. 0) THEN
            CALL PRTERR ('CMDERR', 'Expected variable name')
            GOTO 100
         END IF

C     --Get the ALL (versus this time step) field

         CALL FFCHAR (IFLD, INTYP, CFIELD, 'THIS', WORD)
         IF (WORD .EQ. 'ALL') THEN
            MMSTEP = 0
         ELSE IF (WORD .EQ. 'THIS') THEN
            MMSTEP = NSTEP
         ELSE
            CALL PRTERR ('CMDERR', 'Expected "ALL" or "THIS"')
            GOTO 100
         END IF

         MMNAME = NAME
         MMVAR = IV
         MMTYP = CH
         MMNUM = 1
      ELSE

C     --Get next min/max for last variable

         IF (MMVAR .EQ. 0) THEN
            CALL PRTERR ('CMDERR', 'Expected variable name')
            GOTO 100
         END IF
      END IF

      RETURN

 100  CONTINUE
      RETURN 1
      END
