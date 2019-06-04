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

C     -*- Mode: fortran -*-
      SUBROUTINE GETSPL(A)
      DIMENSION A(*)

C $Id: getspl.f,v 1.6 1997/11/26 16:28:38 gdsjaar Exp $
C $Log: getspl.f,v $
C Revision 1.6  1997/11/26 16:28:38  gdsjaar
C Added NOSCALE option to the spline sweeping. If specified, then the
C input 2D mesh is not stretched to match the spline extent. The output
C X,Y coordinates are the same as the input X,Y coordinates.
C
C Revision 1.5  1994/08/25 13:45:24  gdsjaar
C Fixed spline input to skip back or bottom. Added number of points
C output for easier debugging.
C
c Revision 1.4  1993/05/27  22:15:32  gdsjaar
c Added xsweep and ysweep options to spline transformation
c
c Revision 1.3  1991/03/28  15:39:11  gdsjaar
c Modified spline from GEN3D to work with GENSHELL
c
c Revision 1.2  1990/11/09  15:44:08  gdsjaar
c Added help option
c

      INCLUDE 'gs_splxyz.blk'

      PARAMETER (BINGO = 1.0E38)
      PARAMETER (MAXFLD = 10)
      PARAMETER (NDEFS = 64)
      CHARACTER*8 WORD, VERB
      INTEGER     INTYP(MAXFLD)
      CHARACTER*8 CFIELD(MAXFLD)
      INTEGER     IFIELD(MAXFLD)
      REAL        RFIELD(MAXFLD)

      LOGICAL HELP, ISHELP

      CHARACTER*8 CMDTBL(17)
      SAVE CMDTBL
C     --CMDTBL - the valid commands table

C     --Command table follows.  Remember to change the dimensioned size when
C     --changing the table.
      DATA CMDTBL /
     1     'LEFT    ', 'RIGHT   ', 'ANGULAR ', 'TOP     ', 'FRONT   ',
     $     'BOTTOM  ', 'BACK    ', 'END     ', 'EXIT    ', 'LIST    ',
     3     'HELP    ', 'SPHERICA', 'XSWEEP  ', 'YSWEEP  ',
     *     'SCALE   ', 'NOSCALE ', '        ' /

      CALL SHOCMD ('COMMANDS', CMDTBL)

C     ... Initialize default values

      RDTHET = .FALSE.
      IPTA = 0
      IPTB = 0
      SLLFT(1) = BINGO
      SLRGT(1) = BINGO
      SLLFT(2) = BINGO
      SLRGT(2) = BINGO
      ISPL = 1
      SWEEP = SPHERI

C     ... Allocate arrays, guess on amount and increase if more points entered.

      CALL MDRSRV ('RSPLA',  KRSPLA,  NDEFS)
      CALL MDRSRV ('ZSPLA',  KZSPLA,  NDEFS)
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) THEN
         CALL MEMERR
         STOP 'Memory Allocation in GETSPL'
      END IF

      NTOTA = NDEFS
   10 CONTINUE

C     --Read command line

      WRITE (*, *)
      CALL FREFLD (0, 0, '   Spline Option > ', MAXFLD,
     &     IOSTAT, NUMFLD, INTYP, CFIELD, IFIELD, RFIELD)
      IF (IOSTAT .LT. 0) GOTO 20
      IF (NUMFLD .EQ. 0) GOTO 10
      INTYP(MIN(MAXFLD,NUMFLD)+1) = -999

      IFLD = 1
      CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)
      CALL ABRSTR (VERB, WORD, CMDTBL)
      IF (VERB .EQ. ' ') VERB = WORD

      IF (VERB .EQ. 'EXIT' .OR. VERB .EQ. 'END') GO TO 20

      IF (INTYP(1) .EQ. 0) THEN
         IF ( VERB .EQ. 'LEFT') THEN
            CALL FFREAL (IFLD, INTYP, RFIELD,
     &           'slope at left end', 0.0, SLLFT(ISPL), *10)
         ELSE IF (VERB .EQ. 'RIGHT') THEN
           CALL FFREAL (IFLD, INTYP, RFIELD,
     &       'slope at right end', 0.0, SLRGT(ISPL), *10)
         ELSE IF (VERB .EQ. 'SCALE') THEN
           NOSCAL = .FALSE.
         ELSE IF (VERB .EQ. 'NOSCALE') THEN
           NOSCAL = .TRUE.
           IF (RDTHET) THEN
             CALL PRTERR('CMDERR',
     *         'NOSCALE option cannot be used with ANGULAR spline')
           END IF
         ELSE IF (VERB .EQ. 'ANGULAR') THEN
            RDTHET = .TRUE.
         ELSE IF (VERB .EQ. 'SPHERICAL') THEN
            SWEEP = SPHERI
         ELSE IF (VERB .EQ. 'XSWEEP') THEN
            SWEEP = XSWEEP
         ELSE IF (VERB .EQ. 'YSWEEP') THEN
            SWEEP = YSWEEP
         ELSE IF (VERB .EQ. 'TOP' .OR. VERB .EQ. 'FRONT') THEN
            ISPL = 1
         ELSE IF (VERB .EQ. 'BOTTOM*' .OR. VERB .EQ. 'BACK*') THEN
            ISPL = 2
         ELSE IF (VERB .EQ. 'HELP') THEN
            ISHELP = HELP ('GEN3D', 'COMMANDS', CFIELD(IFLD))
            IF (.NOT. ISHELP) CALL SHOCMD ('COMMANDS', CMDTBL)
            VERB = ' '

         ELSE IF (VERB .EQ. 'LIST') THEN
            CALL SHOCMD ('COMMANDS', CMDTBL)
         END IF
      ELSE IF (NUMFLD .GE. 2) THEN
         IF (ISPL .EQ. 1) THEN
            IPTA = IPTA + 1
            IF (IPTA .GT. NTOTA) THEN
               NTOTA = NTOTA + NDEFS
               CALL MDLONG ('RSPLA', KRSPLA, NTOTA)
               CALL MDLONG ('ZSPLA', KZSPLA, NTOTA)
               CALL MDSTAT (NERR, MEM)
               IF (NERR .GT. 0) THEN
                  CALL MEMERR
                  STOP 'Memory Allocation in GETSPL'
               END IF
            END IF
            A(KRSPLA + IPTA - 1) = RFIELD(1)
            A(KZSPLA + IPTA - 1) = RFIELD(2)
         ELSE
C ... NOTE: ISPL = 2 Retained so that spline files generated for GEN3D
C           can also be used for GENSHELL
            CONTINUE
         END IF
      END IF
      GO TO 10

C     ... Done reading data, compress arrays and allocate remaining arrays

   20 CONTINUE
      CALL MDLONG ('RSPLA',  KRSPLA, IPTA)
      CALL MDLONG ('ZSPLA',  KZSPLA, IPTA)
      CALL MDRSRV ('ZSPL2A', KSPL2A, IPTA)
      CALL MDRSRV ('SCRA',   KSCRA,  IPTA)
      CALL MDRSRV ('DISTA',  KDISTA, IPTA)

      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) THEN
         CALL MEMERR
         STOP 'Memory Allocation in GETSPL'
      END IF

      NSPL(1) = IPTA

      RETURN
      END
