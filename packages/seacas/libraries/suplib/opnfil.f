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
      SUBROUTINE OPNFIL (IUNIT, INOUT, FFORM, IDAPAR, IERR)
C=======================================================================
C$Id: opnfil.f,v 1.4 2009/03/25 12:46:02 gdsjaar Exp $
C   --*** OPNFIL *** (ETCLIB) Open a file
C   --   Written by Amy Gilkey - revised 04/20/88
C   --
C   --OPNFIL opens an input or output file.  The file name is retrieved
C   --with a call to EXNAME.  The file is rewound after open.
C   --
C   --Parameters:
C   --   IUNIT - IN - the file unit number
C   --   INOUT - IN -
C   --      'I' if input file (status = "old")
C   --      'O' if output file (status = "new")
C   --      'S' if scratch file (status = "scratch")
C   --      'U' if opened with status = "unknown"
C   --   FFORM - IN - the file type:
C   --      'F' if formatted file
C   --      'U' if unformatted file
C   --      'L' if formatted, no carriage control file
C   --   IDAPAR - IN - direct access parameters:
C   --      (1) = number of words in record; 0 if not direct access
C   --      (2) = number of records
C   --   IERR - OUT - the open error status (0 = no error)

C   --Routines Called:
C   --   EXNAME - (SUPES) Get filename from unit
C   --   LENSTR - (STRLIB) Get string length

      CHARACTER INOUT
      CHARACTER FFORM
      INTEGER IDAPAR(2)

      CHARACTER*1024 FILNAM
      CHARACTER*11 FORM
      CHARACTER*7 STAT
      CHARACTER*8 CDUM

      IERR = 999
      CALL EXNAME (IUNIT, FILNAM, LNAM)

      IF (FFORM .EQ. 'U') THEN
         FORM = 'UNFORMATTED' 
      ELSE
         FORM = 'FORMATTED' 
      END IF
      LFORM = LENSTR (FORM)

      IF (INOUT .EQ. 'I') THEN
         STAT = 'OLD' 
      ELSE IF (INOUT .EQ. 'O') THEN
         STAT = 'NEW' 
      ELSE IF (INOUT .EQ. 'U') THEN
         STAT = 'UNKNOWN'
      ELSE IF (INOUT .EQ. 'S') THEN
         STAT = 'SCRATCH' 
      ELSE
         IERR = -999
         GOTO 10
      END IF
      LSTAT = LENSTR (STAT)

      NINREC = IDAPAR(1)
      IF (NINREC .EQ. 0) THEN
         OPEN (UNIT=IUNIT, FILE=FILNAM(:LNAM), FORM=FORM(:LFORM),
     &      STATUS=STAT(:LSTAT), IOSTAT = IERR)
         IF (IERR .NE. 0 .AND. INOUT .EQ. 'I') THEN
            call screrr(iunit, filnam, lnam, '...', '...')
         END IF
      ELSE
         NRECS = IDAPAR(2)
         CALL EXPARM (CDUM, CDUM, IDUM, KSCU, KNSU, IDAU)
         IF (IDAU .EQ. 0) THEN
            NICE = 512
         ELSE
            NICE = (512 * KNSU - KNSU + 1) / KSCU
         END IF
         NBLKS = (NINREC + NICE-1) / NICE
         NINREC = NBLKS * NICE
         IF (IDAU .EQ. 0) THEN
            LREC = (NINREC * KSCU + KNSU-1) / KNSU
         ELSE
            LREC = NINREC
         END IF

         IF (INOUT .EQ. 'S') THEN
            STAT = 'unknown'
            LSTAT = LENSTR (STAT)
         END IF
         NSIZE = NINREC * NRECS
         OPEN (UNIT=IUNIT, FILE=FILNAM(:LNAM), FORM=FORM(:LFORM),
     &      ACCESS='DIRECT', RECL=LREC,
     &      STATUS=STAT(:LSTAT), ERR=10)

      END IF
 30   CONTINUE

      IERR = 0
 10   CONTINUE
      RETURN
      END
