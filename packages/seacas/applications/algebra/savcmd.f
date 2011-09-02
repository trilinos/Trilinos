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
      SUBROUTINE SAVCMD (INLINE, INTYP, CFIELD, NAMES, *)
C=======================================================================

C   --*** SAVCMD *** (ALGEBRA) Perform SAVE command
C   --   Written by Amy Gilkey - revised 02/22/88
C   --
C   --SAVCMD processes the input SAVE command.  It adds all the variables
C   --to be saved to the /VAR../ arrays.
C   --
C   --Parameters:
C   --   INLINE - IN/OUT - the parsed input lines for the log file
C   --   INTYP - IN - the field types
C   --   CFIELD - IN - the character fields
C   --   NAMES - IN - the global, nodal, and element variable names
C   --   * - return statement if command not executed
C   --
C   --Common Variables:
C   --   Sets NUMINP, IXLHS of /VAR../
C   --   Uses NVARHI, NVARGL, NVARNP, NVAREL of /DBNUMS/

      include 'params.blk'
      include 'namlen.blk'
      include 'var.blk'
      include 'dbnums.blk'

      CHARACTER*(*) INLINE
      INTEGER INTYP(*)
      CHARACTER*(*) CFIELD(*)
      CHARACTER*(*) NAMES(*)

      LOGICAL FFEXST
      CHARACTER*(maxnam) WORD, NAME
      CHARACTER TYPABR
      CHARACTER*5 STRA, STRB

      CHARACTER*(MXSTLN) TYPTBL(5)
      SAVE TYPTBL
      DATA TYPTBL /
     &  'GLOBALS                         ',
     *  'NODALS                          ',
     *  'ELEMENTS                        ',
     *  'ALL                             ',
     &  '                                ' /
c      DATA TYPTBL /
c     &   'HISTORY ', 'GLOBALS ', 'NODALS  ', 'ELEMENTS', 'ALL     ',
c     &   '        ' /

C   --Save the /VAR../ indices so they can be restored in case of error
      NINP = NUMINP
      ILHS = IXLHS

      IF (.NOT. FFEXST (1, INTYP)) THEN
         CALL PRTERR ('CMDERR', 'No options on SAVE command')
         GOTO 120
      END IF

      IFLD = 1
  100 CONTINUE
      IF (FFEXST (IFLD, INTYP)) THEN

         CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', NAME)

         IX = LOCSTR (NAME, NVARGL+NVARNP+NVAREL, NAMES)
         IF (IX .GT. 0) THEN
            CALL FFADDC (NAME, INLINE)
            CALL DBVTYP (IX, TYPABR, IVAR)
         ELSE
            IVAR = 0
            if ((name .eq. 'NODE') .or. (name .eq. 'NODES'))
     &         name = 'NODAL'
            CALL ABRSTR (WORD, NAME, TYPTBL)
            TYPABR = WORD(1:1)
            IF (TYPABR .EQ. ' ') THEN
               CALL PRTERR ('CMDWARN', 'Invalid SAVE option "'
     &            // NAME(:LENSTR(NAME)) // '", ignored')
               GOTO 110
            END IF
            CALL FFADDC (WORD, INLINE)
         END IF

         IF ((TYPABR .EQ. 'G') .OR. (TYPABR .EQ. 'A')) THEN
            IF (IVAR .EQ. 0) THEN
               CALL DBVIX ('G', 1, IGV)
               CALL ADDVAR (NAMES(IGV), NVARGL, 'G', 1,
     &            NINP, ILHS)
            ELSE
               CALL ADDVAR (NAME, 1, 'G', 1, NINP, ILHS)
            END IF
         END IF

         IF ((TYPABR .EQ. 'N') .OR. (TYPABR .EQ. 'A')) THEN
            IF (IVAR .EQ. 0) THEN
               CALL DBVIX ('N', 1, INV)
               CALL ADDVAR (NAMES(INV), NVARNP, 'N', 1,
     &            NINP, ILHS)
            ELSE
               CALL ADDVAR (NAME, 1, 'N', IVAR, NINP, ILHS)
            END IF
         END IF

         IF ((TYPABR .EQ. 'E') .OR. (TYPABR .EQ. 'A')) THEN
            IF (IVAR .EQ. 0) THEN
               CALL DBVIX ('E', 1, IEV)
               CALL ADDVAR (NAMES(IEV), NVAREL, 'E', 1,
     &            NINP, ILHS)
            ELSE
               CALL ADDVAR (NAME, 1, 'E', IVAR, NINP, ILHS)
            END IF
         END IF

  110    CONTINUE
         GOTO 100
      END IF

      IF (NUMINP .GE. IXLHS) THEN
         N = NUMINP + (MAXVAR - IXLHS + 1)
         CALL INTSTR (1, 0, N, STRA, LSTRA)
         CALL INTSTR (1, 0, MAXVAR, STRB, LSTRB)
         CALL PRTERR ('CMDSPEC',
     &      'Too many variable names to store, '
     &      // STRA(:LSTRA) // ' > ' // STRB(:LSTRB))
         CALL PRTERR ('CMDERR', 'SAVE command ignored')

         NUMINP = NINP
         IXLHS = ILHS

         GOTO 120
      END IF

      RETURN

  120 CONTINUE
      RETURN 1
      END
