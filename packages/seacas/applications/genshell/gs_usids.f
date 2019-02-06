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
      SUBROUTINE USIDS (IFLD, INTYP, CFIELD, IFIELD,
     &   LOLD1, IDOLD1, LOLD2, IDOLD2, LENNEW, IDNEW, *)
C=======================================================================

C   $Id: usids.f,v 1.1 1990/08/20 12:23:20 gdsjaar Exp $
C   $Log: usids.f,v $
C   Revision 1.1  1990/08/20 12:23:20  gdsjaar
C   Initial revision
C

C   --*** USIDS *** (GEN3D) Read list of IDs
C   --   Written by Amy Gilkey - revised 05/21/86
C   --
C   --USIDS processes a list of IDs.  The ID list is checked for
C   --repetitions with the existing list and with itself.  Repetitions
C   --are flagged with an error message and ignored.
C   --
C   --Parameters:
C   --   IFLD - IN/OUT - the free-field index
C   --   INTYP - IN - the free-field type
C   --   CFIELD - IN - the free-field characters
C   --   IFIELD - IN - the free-field integers
C   --   LOLD1, LOLD2 - IN - the length of the existing IDs
C   --   IDOLD1, IDOLD2 - IN - the existing IDs
C   --   LENNEW - OUT - the length of the returned IDs
C   --   IDNEW - OUT - the returned IDs
C   --   * - return statement iff serious error

      PARAMETER (MAXSET=10)

      INTEGER INTYP(*)
      CHARACTER*8 CFIELD(*)
      INTEGER IFIELD(*)
      INTEGER IDOLD1(*), IDOLD2(*)
      INTEGER IDNEW(*)

      CHARACTER*5 STRA
      LOGICAL DUPID

      CALL INIINT (MAXSET, 0, IDNEW)
      LENNEW = 0

   10 CONTINUE
      IF (INTYP(IFLD) .GE. -1) THEN

         CALL FFINTG (IFLD, INTYP, IFIELD,
     &      'set id', 0, ID, *50)

         IF (LENNEW .GE. MAXSET) THEN
            CALL INTSTR (1, 0, MAXSET, STRA, LSTRA)
            CALL PRTERR ('CMDERR',
     &         'Number of IDs must be less than ' // STRA(:LSTRA))
            GOTO 60
         END IF

         DUPID = .FALSE.
         DO 20 I = 1, LOLD1
            IF (ID .EQ. IDOLD1(I)) DUPID = .TRUE.
   20    CONTINUE
         DO 30 I = 1, LOLD2
            IF (ID .EQ. IDOLD2(I)) DUPID = .TRUE.
   30    CONTINUE
         DO 40 I = 1, LENNEW
            IF (ID .EQ. IDNEW(I)) DUPID = .TRUE.
   40    CONTINUE
         IF (DUPID) THEN
            CALL INTSTR (1, 0, ID, STRA, LSTRA)
            CALL PRTERR ('CMDWARN',
     &         'Duplicate ID ' // STRA(:LSTRA) // ' ignored')
            GOTO 50
         END IF

         LENNEW = LENNEW + 1
         IDNEW(LENNEW) = ID
   50    CONTINUE
         GOTO 10
      END IF

   60 CONTINUE
      RETURN
      END
