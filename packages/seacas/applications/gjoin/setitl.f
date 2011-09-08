C Copyright (c) 2008 Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
C certain rights in this software
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
C 
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
C 

C -*- Mode: fortran -*-
C=======================================================================
C $Id: setitl.f,v 1.1 1999/01/18 19:21:26 gdsjaar Exp $
C $Log: setitl.f,v $
C Revision 1.1  1999/01/18 19:21:26  gdsjaar
C ExodusII version of gjoin, needs testing and syncing with exodus 1 version, but is being committed to permit easier testing and modifications.  This was created by Dave Fry at Goodyear
C
c Revision 1.1.1.1  1998/11/05  16:23:28  a294617
c Initial import == gjoin 1.36
c
C Revision 1.2  1997/04/04 20:06:43  gdsjaar
C Better command input in setitl (TITLE submenu). It was using getinp
C instead of FREFLD, but that caused problems if the command contained
C leading spaces. Replaced to use FREFLD.
C
C Revision 1.1.1.1  1990/11/12 14:36:05  gdsjaar
C GJOIN - X1.00.40 - 7/17/90
C
c Revision 1.1  90/11/12  14:36:03  gdsjaar
c Initial revision
c 
      SUBROUTINE SETITL (TWODB)
C=======================================================================

C   --*** SETITL *** (GJOIN) Select database title
C   --   Written by Amy Gilkey - revised 02/23/88
C   --
C   --SETITL selects the database title from the existing titles and
C   --user-supplied instructions.
C   --
C   --Parameters:
C   --   TWODB - IN - true iff two databases

      include 'exodusII.inc'
      include 'params.blk'
      include 'titles.blk'
      INCLUDE 'filnum.blk'

      PARAMETER (MAXFLD=1)
      LOGICAL TWODB

      CHARACTER*8 WORD, VERB
      INTEGER INTYP(MAXFLD+1)
      CHARACTER*8 CFIELD(MAXFLD)
      INTEGER IFIELD(MAXFLD)
      REAL RFIELD(MAXFLD)

      CHARACTER*8 CMDTBL(8)
      SAVE CMDTBL
C      --CMDTBL - the valid commands table

C   --Command table follows.  Remember to change the dimensioned size when
C   --changing the table.
      DATA CMDTBL /
     1   '1       ', '2       ', 'CHANGE  ',
     2   'LIST    ', 'HELP    ', 'UP      ', 'EXIT    ',
     3   '        ' /

C   --Print the input database titles and the output database title.

   50 CONTINUE

      WRITE (*, *)
      IF ((.NOT. TWODB) .OR. (TITLE1 .EQ. TITLE2)) THEN
         WRITE (*, 55) 'Database title:'
         WRITE (*, 55) TITLE1(:LENSTR(TITLE1))
      ELSE
         WRITE (*, 55) 'Database titles:'
         WRITE (*, 55) TITLE1(:LENSTR(TITLE1))
         WRITE (*, 55) TITLE2(:LENSTR(TITLE2))
      END IF
      WRITE (*, 55) 'Output database title:'
      WRITE (*, 55) TITLE(:LENSTR(TITLE))
   55 FORMAT (1X, 5A)

  100 CONTINUE

      WRITE (*, *)
      CALL FREFLD (0, 0, 'TITLE> ', MAXFLD,
     &   IOSTAT, NUMFLD, INTYP, CFIELD, IFIELD, RFIELD)

      IF (IOSTAT .LT. 0) GOTO 110
      IF (NUMFLD .EQ. 0) GOTO 100
      INTYP(MIN(MAXFLD,NUMFLD)+1) = -999

      IFLD = 1
      CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)
      CALL ABRSTR (VERB, WORD, CMDTBL)
      IF (VERB .EQ. ' ') VERB = WORD
      CFIELD(IFLD-1) = VERB
      INTYP(IFLD-1)  = 0
      CALL OUTLOG (KLOG, MIN(MAXFLD,NUMFLD), INTYP,
     $     CFIELD, IFIELD, RFIELD)

      IF (VERB .EQ. '1') THEN
         TITLE = TITLE1

         GOTO 50

      ELSE IF (VERB .EQ. '2') THEN
         IF (TWODB) THEN
            TITLE = TITLE2
         ELSE
            TITLE = TITLE1
         END IF

         GOTO 50

      ELSE IF (VERB .EQ. 'CHANGE') THEN
         CALL GETINP (0, 0, 'New title> ', TITLE, IOSTAT)
         CALL OUTLOG (KLOG, 1, 0, TITLE, IDUM, RDUM)

         GOTO 50

      ELSE IF (VERB .EQ. 'LIST') THEN
         GOTO 50

      ELSE IF (VERB .EQ. 'HELP') THEN
         WRITE (*, 10000)
10000     FORMAT (
     &      /,1X,'Valid Commands:',
     &      /,4X,'1  -  copy title from first database',
     &      /,4X,'2  -  copy title from second database (if any)',
     &      /,4X,'CHANGE  -  change title to user-specified title'
     &      /,4X,'LIST  -  list database titles'
     &      /,4X,'UP  -  go up a command level'
     &      )

      ELSE IF (VERB .EQ. 'UP' .OR. VERB .EQ. 'EXIT') THEN
         GOTO 110

      ELSE
         CALL PRTERR ('CMDERR', '"' // VERB(:LENSTR(VERB))
     &      // '" is an invalid command')
      END IF

      GOTO 100

  110 CONTINUE
      RETURN
      END
