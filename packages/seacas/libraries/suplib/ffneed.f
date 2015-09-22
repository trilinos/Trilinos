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

C=======================================================================
      SUBROUTINE FFNEED (IFLD, INTYP, FTYPE, NFLD, EXPECT, *)
C=======================================================================
C$Id: ffneed.f,v 1.3 2009/03/25 12:46:02 gdsjaar Exp $
C$Log: ffneed.f,v $
CRevision 1.3  2009/03/25 12:46:02  gdsjaar
CAdd copyright and license notice to all files.
C
CRevision 1.2  1993/08/19 18:40:54  gdsjaar
CFixed incorrect itype for integers
C
c Revision 1.1.1.1  1990/08/14  16:14:31  gdsjaar
c Testing
c
c Revision 1.1  90/08/14  16:14:29  gdsjaar
c Initial revision
c 
c Revision 1.1  90/08/09  13:39:25  gdsjaar
c Initial revision
c 

C   --*** FFNEED *** (FFLIB) Check free-field fields for type
C   --   Written by Amy Gilkey - revised 10/21/86
C   --
C   --FFNEED checks that the next free-format fields exist and are of the
C   --appropriate type.
C   --
C   --Parameters:
C   --   IFLD - IN - the index of the current field number, NOT incremented
C   --   INTYP - IN - the input types from the free-field reader
C   --   FTYPE - IN - the expected field type:
C   --      C for character, R for real, I for integer, other for character
C   --   NFLD - IN - the number of expected fields
C   --   EXPECT - IN - the value to expect string, for error
C   --   * - return statement if the fields do not exist or are not of the
C   --      expected type; message is printed

      INTEGER IFLD
      INTEGER INTYP(*)
      CHARACTER*(*) FTYPE
      INTEGER NFLD
      CHARACTER*(*) EXPECT

      CHARACTER*80 ERRMSG

      IF (FTYPE(1:1) .EQ. 'R') THEN
         ITYPE = 1
      ELSE IF (FTYPE(1:1) .EQ. 'I') THEN
         ITYPE = 2
      ELSE
         ITYPE = 0
      END IF

      DO 100 ICHK = IFLD, IFLD + NFLD - 1
         IF (INTYP(ICHK) .NE. ITYPE) THEN
            IF ((ITYPE .NE. 1) .OR. (INTYP(ICHK) .NE. 2)) THEN
               ERRMSG = 'Expected ' // EXPECT
               CALL PRTERR ('CMDERR', ERRMSG(:LENSTR(ERRMSG)))
               GOTO 110
            END IF
         END IF
  100 CONTINUE

      RETURN

  110 CONTINUE
      RETURN 1
      END
