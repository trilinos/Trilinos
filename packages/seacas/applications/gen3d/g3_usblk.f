C Copyright(C) 2011-2017 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
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
      SUBROUTINE USBLK (IFLD, INTYP, CFIELD, IFIELD,
     &   NEWTYP, NELBLK, IDELB, BLKTYP, *)
C=======================================================================

C   --*** USBLK *** (GEN3D) Read list of element block IDs
C   --   Written by Amy Gilkey - revised 05/21/86
C   --
C   --USBLK processes a list of element block IDs.  If an element block
C   --is in the list, the block type is changed to the new block type.
C   --
C   --Parameters:
C   --   IFLD - IN/OUT - the free-field index
C   --   INTYP - IN - the free-field type
C   --   CFIELD - IN - the free-field characters
C   --   IFIELD - IN - the free-field integers
C   --   NEWTYP - IN - the element block type to be set
C   --   NELBLK - IN - the number of element blocks
C   --   IDELB - IN - the ids for each block
C   --   BLKTYP - IN/OUT - the element block type
C   --   * - return statement iff serious error

      PARAMETER (MAXSET=10)

      INTEGER INTYP(*)
      CHARACTER*8 CFIELD(*)
      INTEGER IFIELD(*)
      CHARACTER NEWTYP
      INTEGER IDELB(NELBLK)
      CHARACTER BLKTYP(NELBLK)

      LOGICAL FFEXST
      CHARACTER*5 ISTR

      IF (.NOT. FFEXST (IFLD, INTYP)) THEN
         CALL INISTR (NELBLK, NEWTYP, BLKTYP)
      END IF

   10 CONTINUE
      IF (FFEXST (IFLD, INTYP)) THEN
         CALL FFINTG (IFLD, INTYP, IFIELD,
     &      'block id', 0, ID, *20)
         IELB = LOCINT (ID, NELBLK, IDELB)
         IF (IELB .LE. 0) THEN
            CALL INTSTR (1, 0, ID, ISTR, LSTR)
            CALL PRTERR ('CMDERR',
     &         'Invalid block id ' // ISTR(:LSTR) // ', ignored')
            GOTO 20
         END IF
         BLKTYP(IELB) = NEWTYP
   20    CONTINUE
         GOTO 10
      END IF

      RETURN
      END
