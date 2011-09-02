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
      SUBROUTINE FFSFLD (TYPE, WORD, INUM, RNUM,
     &   IFLD, INTYP, CFIELD, IFIELD, RFIELD)
C=======================================================================
C$Id: ffsfld.f,v 1.2 2009/03/25 12:46:02 gdsjaar Exp $
C$Log: ffsfld.f,v $
CRevision 1.2  2009/03/25 12:46:02  gdsjaar
CAdd copyright and license notice to all files.
C
CRevision 1.1.1.1  1990/08/14 16:14:45  gdsjaar
CTesting
C
c Revision 1.1  90/08/14  16:14:43  gdsjaar
c Initial revision
c 
c Revision 1.1  90/08/09  13:39:27  gdsjaar
c Initial revision
c 

C   --*** FFSFLD *** (FFLIB) Set free-field field
C   --   Written by Amy Gilkey - revised 05/07/87
C   --
C   --FFSFLD sets a free-field reader field.
C   --
C   --Parameters:
C   --   TYPE - IN - the field type:
C   --      '*' = end of fields field
C   --      ' ' = empty field
C   --      'C' = character field
C   --      'R' = real number field
C   --      'I' = integer number field
C   --   WORD - IN - the character field (if type 'C')
C   --   INUM - IN - the integer field (if type 'I')
C   --   RNUM - IN - the real field (if type 'R')
C   --   IFLD - IN - the index of the field to be set
C   --   INTYP - IN/OUT - the free-field reader type array
C   --   CFIELD - IN/OUT - the free-field reader string array
C   --   IFIELD - IN/OUT - the free-field reader integer array
C   --   RFIELD - IN/OUT - the free-field reader real array

      CHARACTER TYPE
      CHARACTER*(*) WORD
      INTEGER INUM
      REAL RNUM
      INTEGER IFLD
      INTEGER INTYP(*)
      CHARACTER*(*) CFIELD(*)
      INTEGER IFIELD(*)
      REAL RFIELD(*)

      CHARACTER*32 NUMBCH

      IF (TYPE .EQ. '*') THEN
         INTYP(IFLD) = -999
      ELSE IF (TYPE .EQ. ' ') THEN
         INTYP(IFLD) = -1
         CFIELD(IFLD) = ' '
         IFIELD(IFLD) = 0
         RFIELD(IFLD) = 0
      ELSE IF (TYPE .EQ. 'C') THEN
         INTYP(IFLD) = 0
         CFIELD(IFLD) = WORD
         IFIELD(IFLD) = 0
         RFIELD(IFLD) = 0
      ELSE IF (TYPE .EQ. 'I') THEN
         INTYP(IFLD) = 2
         WRITE (NUMBCH, '(I32)', IOSTAT=K) INUM
         CALL SQZSTR (NUMBCH, L)
         CFIELD(IFLD) = NUMBCH
         IFIELD(IFLD) = INUM
         RFIELD(IFLD) = INUM
      ELSE IF (TYPE .EQ. 'R') THEN
         INTYP(IFLD) = 1
         WRITE (NUMBCH, '(E8.3E1)', IOSTAT=K) RNUM
         CALL SQZSTR (NUMBCH, L)
         CFIELD(IFLD) = NUMBCH
         IFIELD(IFLD) = 0
         RFIELD(IFLD) = RNUM
      END IF

      RETURN
      END
