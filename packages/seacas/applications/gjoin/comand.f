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
C $Id: comand.f,v 1.3 2001/06/26 17:38:54 gdsjaar Exp $

      SUBROUTINE COMAND (A, IELBST, IDELB, NUMELB, NUMLNK, NUMATR,
     &  NAMELB, INPSST, IDNPS, NNNPS, IESSST, IDESS, NEESS, DONE, *)
C=======================================================================

C   --*** COMAND *** (GJOIN) Command input and execution
C   --   Written by Amy Gilkey - revised 02/23/88
C   --
C   --COMAND inputs and executes an user command.
C   --
C   --Parameters:
C   --   A - IN - the dynamic memory base array
C   --   IELBST - IN/OUT - the element block status
C   --   IDELB - IN/OUT - the element block ID for each block
C   --   NUMELB - IN - the number of elements for each block
C   --   NUMLNK - IN - the number of nodes per element for each block
C   --   NUMATR - IN - the number of attributes for each block
C   --   INPSST - IN/OUT - the nodal point set status
C   --   IDNPS - IN/OUT - the nodal point set ID for each set
C   --   NNNPS - IN - the number of nodes for each set
C   --   IESSST - IN/OUT - the element side set status
C   --   IDESS - IN/OUT - the element side set ID for each set
C   --   NEESS - IN - the number of elements for each set
C   --   DONE - OUT - true iff no more input files
C   --   * - return statement if quit
C   --
C   --Common Variables:
C   --   Sets and uses /TITLES/
C   --   Uses /DBVARS/

      PARAMETER (MAXFLD=80)

      include 'exodusII.inc'

      include 'params.blk'
      INCLUDE 'filnum.blk'
      INCLUDE 'titles.blk'
      INCLUDE 'dbvars.blk'

      DIMENSION A(*)
      INTEGER IELBST(*)
      INTEGER IDELB(*)
      INTEGER NUMELB(*)
      INTEGER NUMLNK(*)
      INTEGER NUMATR(*)
      CHARACTER*(MXSTLN) NAMELB(*)
      INTEGER INPSST(*)
      INTEGER IDNPS(*)
      INTEGER NNNPS(*)
      INTEGER IESSST(*)
      INTEGER IDESS(*)
      INTEGER NEESS(*)
      LOGICAL DONE

      CHARACTER*8 WORD, VERB, STRA, STRB
      INTEGER     INTYP(MAXFLD+1)
      CHARACTER*8 CFIELD(MAXFLD)
      INTEGER     IFIELD(MAXFLD)
      REAL        RFIELD(MAXFLD)
      LOGICAL ISON

      CHARACTER*8 CMDTBL(14)
      SAVE CMDTBL
C      --CMDTBL - the valid commands table

C   --Command table follows.  Remember to change the dimensioned size when
C   --changing the table.
      DATA CMDTBL /
     $   'TITLE   ', 'BLOCKS  ', 'MATERIAL', 'NSETS   ', 'SSETS   ',
     $   'HELP    ', 'FINISH  ', 'EXIT    ', 'END     ', 'QUIT    ',
     $   'NODESETS', 'SIDESETS', 'ADD     ', '        ' /

      WRITE (*, *)

C   --Initialize the database title

      TITLE = TITLE1

C   --Initialize the element blocks

      DO 100 IELB = 1, NEWELB
         IELBST(IELB) = 0
  100 CONTINUE
      ISON = .FALSE.
      DO 110 IELB = 1, NEWELB
         I = LOCINT (IDELB(IELB), NEWELB, IDELB)
         IF (I .LT. IELB) THEN
            ISON = .TRUE.
            IELBST(IELB) = I
            IELBST(I) = I
C ... Check for compatible element blocks if combining...
            IF ( (NUMLNK(I) .NE. NUMLNK(IELB)) .OR.
     &           (NUMATR(I) .NE. NUMATR(IELB)) ) THEN
               CALL INTSTR (1, -1, I,    STRA, LSTRA)
               CALL INTSTR (1, -1, IELB, STRB, LSTRB)
               CALL PRTERR ('CMDWARN',
     &              'Element blocks '//STRA(:LSTRA)//' and '
     &              //STRB(:LSTRB)//
     &              ' have duplicate IDS but are not compatible')
             ELSE IF (NAMELB(I) .NE. NAMELB(IELB)) THEN
               CALL INTSTR (1, -1, I,    STRA, LSTRA)
               CALL INTSTR (1, -1, IELB, STRB, LSTRB)
               CALL PRTERR ('CMDWARN',
     &              'Element blocks '//STRA(:LSTRA)//' and '
     &              //STRB(:LSTRB)//
     &              ' have duplicate IDS but different element types')
            END IF
         END IF
  110 CONTINUE
      IF (ISON)
     &   CALL PRTERR ('CMDWARN', 'Duplicate IDs in element blocks'
     &   // ' - combined unless changed')

C   --Initialize the nodal point sets

      DO 120 INPS = 1, NEWNPS
         INPSST(INPS) = 0
  120 CONTINUE
      ISON = .FALSE.
      DO 130 INPS = 1, NEWNPS
         I = LOCINT (IDNPS(INPS), NEWNPS, IDNPS)
         IF (I .LT. INPS) THEN
            ISON = .TRUE.
            INPSST(INPS) = I
            INPSST(I) = I
         END IF
  130 CONTINUE
      IF (ISON)
     &   CALL PRTERR ('CMDWARN', 'Duplicate IDs in nodal point sets'
     &   // ' - combined unless changed')

C   --Initialize the element side sets

      DO 140 IESS = 1, NEWESS
         IESSST(IESS) = 0
  140 CONTINUE
      ISON = .FALSE.
      DO 150 IESS = 1, NEWESS
         I = LOCINT (IDESS(IESS), NEWESS, IDESS)
         IF (I .LT. IESS) THEN
            ISON = .TRUE.
            IESSST(IESS) = I
            IESSST(I) = I
         END IF
  150 CONTINUE
      IF (ISON)
     &   CALL PRTERR ('CMDWARN', 'Duplicate IDs in element side sets'
     &   // ' - combined unless changed')

  160 CONTINUE

      DONE = .FALSE.

C   --Read command line

      WRITE (*, *)
      CALL FREFLD (0, 0, 'GJOIN> ', MAXFLD,
     &   IOSTAT, NUMFLD, INTYP, CFIELD, IFIELD, RFIELD)
      IF (IOSTAT .LT. 0) GOTO 210
      IF (NUMFLD .EQ. 0) GOTO 160

      INTYP(MIN(MAXFLD,NUMFLD)+1) = -999

      IFLD = 1
      CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)
      CALL ABRSTR (VERB, WORD, CMDTBL)
      IF (VERB .EQ. ' ') VERB = WORD
      CFIELD(IFLD-1) = VERB
      CALL OUTLOG (KLOG, NUMFLD, INTYP, CFIELD, IFIELD, RFIELD)

      IF (VERB .EQ. 'TITLE') THEN

C      --Allow user to select the database title

  170    CONTINUE
         CALL SETITL (TWODB)

      ELSE IF ((VERB .EQ. 'BLOCKS') .OR. (VERB .EQ. 'MATERIAL')) THEN

C      --Allow user to change element blocks

         CALL CKNONE (NEWELB, .FALSE., 'element blocks', *210)

         CALL MDRSRV ('ISCR', KISCR, NEWELB)
         CALL MDSTAT (NERR, MEM)
         IF (NERR .GT. 0) GOTO 210

  180    CONTINUE
         CALL PRTELB (IELBST, NELBL1, NELBL2,
     &      IDELB, NUMELB, NUMLNK, NUMATR, A(KISCR))
         CALL SETSTA ('BLOCKS> ', 'element block',
     &      NELBL1, NELBL2, IDELB, IELBST, *180)

         CALL MDDEL ('ISCR')

      ELSE IF (VERB .EQ. 'NSETS' .or. VERB .EQ. 'NODESETS') THEN

C      --Allow user to change nodal point sets

         CALL CKNONE (NEWNPS, .FALSE., 'nodal point sets', *210)

         CALL MDRSRV ('ISCR', KISCR, NEWNPS)
         CALL MDSTAT (NERR, MEM)
         IF (NERR .GT. 0) GOTO 210

  190    CONTINUE
         CALL PRTNPS (INPSST, NNPS1, NNPS2, IDNPS, NNNPS,
     &      A(KISCR))
         CALL SETSTA ('NSETS> ', 'nodal point set',
     &      NNPS1, NNPS2, IDNPS, INPSST, *190)

         CALL MDDEL ('ISCR')

      ELSE IF (VERB .EQ. 'SSETS' .or. VERB .eq. 'SIDESETS') THEN

C      --Allow user to change element side sets

         CALL CKNONE (NEWESS, .FALSE., 'element side sets', *210)

         CALL MDRSRV ('ISCR', KISCR, NEWESS)
         CALL MDSTAT (NERR, MEM)
         IF (NERR .GT. 0) GOTO 210

  200    CONTINUE
         CALL PRTESS (IESSST, NESS1, NESS2,
     &      IDESS, NEESS, A(KISCR))
         CALL SETSTA ('SSETS> ', 'element side set',
     &      NESS1, NESS2, IDESS, IESSST, *200)

         CALL MDDEL ('ISCR')

      ELSE IF (VERB .EQ. 'HELP') THEN
         WRITE (*, 10000)
10000     FORMAT (
     &      /,1X,'Valid Commands:',
     &      /,4X,'TITLE  -  change the database title'
     &      /,4X,'BLOCKS  -  manipulate the element blocks'
     &      /,4X,'NSETS  -  manipulate the nodal point sets'
     &      /,4X,'SSETS  -  manipulate the element side sets'
     &      /,4X,'EXIT  -  end command input, start processing'
     $      /,4X,'FINISH - end command input, write output file'
     &      /,4X,'QUIT  -  abort processing'
     &      )

      ELSE IF ((VERB .EQ. 'EXIT') .OR. (VERB .EQ. 'END')) THEN
         CALL PRTERR ('CMDWARN', 'Please use "ADD" or "FINISH"')
         IF (TWODB) THEN
            WRITE (*, *)
            CALL FREFLD (0, 0, 'Is there another database? ', MAXFLD,
     &           IOSTAT, NUMFLD, INTYP, CFIELD, IFIELD, RFIELD)
            CALL OUTLOG (KLOG, NUMFLD, INTYP, CFIELD, IFIELD, RFIELD)
            DONE = (CFIELD(1)(1:1) .NE. 'Y')
         ELSE
            DONE = .TRUE.
         END IF
         CALL SCNEOF
         GOTO 210

      ELSE IF (VERB .EQ. 'FINISH') THEN
         DONE = .TRUE.
         CALL SCNEOF
         GOTO 210

      ELSE IF (VERB .EQ. 'ADD') THEN
         DONE = .FALSE.
         CALL SCNEOF
         GOTO 210

      ELSE IF (VERB .EQ. 'QUIT') THEN
         CALL SCNEOF
         RETURN 1

      ELSE
         CALL PRTERR ('CMDERR', '"' // VERB(:LENSTR(VERB))
     &      // '" is an invalid command')
      END IF

      GOTO 160

  210 CONTINUE
      RETURN
      END
