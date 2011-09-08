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
C     $Id: ffistr.f,v 1.6 2008/03/14 13:22:40 gdsjaar Exp $
      SUBROUTINE FFISTR( LINE,MFIELD,IDCONT,NFIELD,KVALUE,CVALUE,IVALUE,
     *                   RVALUE )
      CHARACTER*(*) LINE,CVALUE(MFIELD)
      CHARACTER*32 CFIELD
      INTEGER KVALUE(MFIELD),IVALUE(MFIELD)
      REAL RVALUE(MFIELD)

      INTEGER TABC
      TABC = 9
C
************************************************************************
C
C     FREFLD INPUT SYSTEM - ANSI FORTRAN - USER INTERFACE ROUTINE
C
C     DESCRIPTION:
C     This routine is the main parsing routine of the SUPES Free Field
C     Input system. It parses a CHARACTER string into data fields, and
C     returns the CHARACTER, REAL, and INTEGER value for each field. A
C     value which indicates whether the character and numeric values
C     were explicitly defined by a valid string or simply set to a
C     default blank or zero is also returned.
C
C     FORMAL PARAMETERS:
C     LINE    CHARACTER  Input string.
C     MFIELD  INTEGER    Maximum number of data fields to be returned.
C     IDCONT  INTEGER    Continuation flag ( 0 = NO continuation )
C     NFIELD  INTEGER    Number of data fields found.
C     KVALUE  INTEGER    Translation states of the data fields:
C                           -1 = This is a null field.
C                            0 = This is a non-numeric field.
C                            1 = This is a REAL numeric field.
C                            2 = This is an INTEGER numeric field.
C     CVALUE  CHARACTER  Character values of the data fields.
C     RVALUE  REAL       Floating-point values of the data fields.
C     IVALUE  INTEGER    Integer values of the data fields.
C
C
C     ROUTINES CALLED:
C     STRIPB            Strip leading/trailing blanks from a string.
C     EXUPCS            Convert a string to ANSI FORTRAN character set.
C     QUOTED            Process a quoted string.
C
************************************************************************
C
C     Initialize output arrays to their default values and zero field
C     counter, unless IDCONT indicates that this is a continuation record.
C
      IF ( IDCONT .EQ. 0 ) THEN
         DO 300 I = 1 , MFIELD
            KVALUE(I) = -1
            CVALUE(I) = ' '
            RVALUE(I) = 0.
            IVALUE(I) = 0
  300    CONTINUE
         NFIELD = 0
      END IF
      IDCONT = 0
C

************************************************************************
C
C     Isolate the effective portion of the input line. At the end
C     of this phase LINE(ILEFT:ISTOP) will represent this portion.
C     The continuation flag IDCONT will indicate whether or not a
C     continuation line is to follow this record. Exit at any
C     point where the effective portion of the line becomes null.
C
      ILEFT = 1
      ISTOP = LEN ( LINE )
C
C     Now start processing fields.
C     Upper range of loop is a dummy maximum.
C
c  Had to fix more VAX FORTRAN specific stuff.  We'll try this only
c  for a short while.  JRR.
c
      DO 1 IFLD = 1, ISTOP
C
         CALL STRIPB( LINE(ILEFT:ISTOP), IL, ISTOP )
         ISTOP = ISTOP + ILEFT - 1
         ILEFT = IL + ILEFT - 1
         IF ( ILEFT .GT. ISTOP ) THEN
C
C           Remainder of line is null.
C
            RETURN
         ELSE IF ( LINE(ILEFT:ILEFT) .EQ. '$' ) THEN
C
C           Rest is comment.
C
            RETURN
         ELSE IF ( LINE(ILEFT:ILEFT) .EQ. '*' ) THEN
C
C           Continuation.
C
            IDCONT = 1
            RETURN
         ELSE IF ( LINE(ILEFT:ILEFT) .EQ. '''' ) THEN
C
C           This is the beginning of a quoted string.  Call a special handler.
C
           CALL QUOTED ( LINE(ILEFT:ISTOP), IL, IRIGHT )
           IF ( IRIGHT .NE. 0 ) IRIGHT = IRIGHT + ILEFT - 1
           ILEFT = IL+ ILEFT - 1
         ELSE IF ( INDEX ( ',=', LINE(ILEFT:ILEFT) ) .NE. 0 ) THEN
C
C           This is a null field.
C
            IRIGHT = 0
         ELSE
C
C           Find the end of this token.  
C           Valid delimiters, are ' ', '*', ',', '=', '$'.
C
            IBLNK = INDEX ( LINE(ILEFT:ISTOP), ' ' ) + ILEFT - 2
            IAST  = INDEX ( LINE(ILEFT:ISTOP), '*' ) + ILEFT - 2
            ICOMA = INDEX ( LINE(ILEFT:ISTOP), ',' ) + ILEFT - 2
            IEQLS = INDEX ( LINE(ILEFT:ISTOP), '=' ) + ILEFT - 2
            IDOLR = INDEX ( LINE(ILEFT:ISTOP), '$' ) + ILEFT - 2
            
            ITAB  = INDEX ( LINE(ILEFT:ISTOP), CHAR(TABC)) + ILEFT - 2
            IF ( IBLNK .LT. ILEFT ) IBLNK = ISTOP + 1
            IF ( IAST  .LT. ILEFT ) IAST  = ISTOP + 1
            IF ( ICOMA .LT. ILEFT ) ICOMA = ISTOP + 1
            IF ( IEQLS .LT. ILEFT ) IEQLS = ISTOP + 1
            IF ( IDOLR .LT. ILEFT ) IDOLR = ISTOP + 1
            IF ( ITAB  .LT. ILEFT ) ITAB = ISTOP + 1
            IRIGHT = MIN ( IBLNK, IAST, ICOMA, IEQLS, IDOLR,
     $           ITAB, ISTOP )
C
C           Convert data to standard character set -
C
            CALL EXUPCS( LINE(ILEFT:IRIGHT) )
         END IF
C
C        Process this field.
C        Don't process this field unless there is room in the data arrays -
C
         NFIELD = NFIELD + 1
         IF ( NFIELD .LE. MFIELD ) THEN
C
C           Calculate the effective length of this field -
C
            LFIELD = IRIGHT - ILEFT + 1
            IF ( LFIELD .LE. 0 ) THEN
C
C              This is a null field; skip it -
C
            ELSE IF ( LFIELD .GT. 32 ) THEN
C
C              This field exceeds the maximum allowable numeric 
C              field size; define only the character value -
C
               CVALUE(NFIELD) = LINE(ILEFT:IRIGHT)
               KVALUE(NFIELD) = 0
            ELSE
C
C              Define the character value for this field, 
C              then right-justify and attempt numeric translations -
C
               CVALUE(NFIELD) = LINE(ILEFT:IRIGHT)
               KVALUE(NFIELD) = 0
               CFIELD = ' '
               IJUST = 32 - LFIELD + 1
               CFIELD(IJUST:32) = LINE(ILEFT:IRIGHT)
C
C              See if a digit is present in this field.
C              If there is no digit present, then do not accept
C              this token as a valid real or integer value.
C              This is needed, since some systems accept a token like
C              "E", "Q", "INF", ... as a valid real number even though
C              no digit appears in the token.
               IF ( INDEX(CFIELD(IJUST:32),'0').NE.0 .OR.
     1              INDEX(CFIELD(IJUST:32),'1').NE.0 .OR.
     2              INDEX(CFIELD(IJUST:32),'2').NE.0 .OR.
     3              INDEX(CFIELD(IJUST:32),'3').NE.0 .OR.
     4              INDEX(CFIELD(IJUST:32),'4').NE.0 .OR.
     5              INDEX(CFIELD(IJUST:32),'5').NE.0 .OR.
     6              INDEX(CFIELD(IJUST:32),'6').NE.0 .OR.
     7              INDEX(CFIELD(IJUST:32),'7').NE.0 .OR.
     8              INDEX(CFIELD(IJUST:32),'8').NE.0 .OR.
     9              INDEX(CFIELD(IJUST:32),'9').NE.0 ) THEN
C ... One more check. If the field contains a digit, but starts with
C     'D' or 'E', several systems will interpret this as a valid
C      Integer and/or real number. This is not the desired behavior
                  IF (CFIELD(IJUST:IJUST) .NE. 'D' .AND.
     &                CFIELD(IJUST:IJUST) .NE. 'E') THEN 
                     IDIG = 1
                  ELSE
                     IDIG = 0
                  END IF
               ELSE
                  IDIG = 0
               END IF
C
C ... It should not be necessary to initialize ITRANS, but the gcc-4.0.0 gfortran
C     Does not correctly set ITRANS after the first execution if it is non-zero.
C     It seems to work correctly if initialized to zero.
               ITRANS = 0
               READ( CFIELD,3000,IOSTAT=ITRANS ) RFIELD
               IF ( IDIG .EQ. 1 .AND. ITRANS .EQ. 0 ) THEN
C
C                 This field has a valid floating-point value -
C
                  RVALUE(NFIELD) = RFIELD
                  KVALUE(NFIELD) = 1
               END IF
               READ( CFIELD,4000,IOSTAT=ITRANS ) IFIELD
               IF ( IDIG .EQ. 1 .AND. ITRANS .EQ. 0 ) THEN
C
C                 This field has a valid integer value -
C
                  IVALUE(NFIELD) = IFIELD
                  KVALUE(NFIELD) = 2
               ELSE IF ( KVALUE(NFIELD) .EQ. 1 .AND. 
     *            ABS ( RVALUE(NFIELD) ) .LE. 1.E9 ) THEN
C
C                 This field has a valid real that did not automatically
C                 Translate to an integer.  Try to convert the real to an 
C                 integer.
C
                  IFIELD = RVALUE(NFIELD)
                  RFIELD = IFIELD
                  IF ( RFIELD .EQ. RVALUE(NFIELD) ) THEN
C
C                 Successful conversion of real to integer.
C
                     IVALUE(NFIELD) = IFIELD
                     KVALUE(NFIELD) = 2
                  END IF
               END IF
            END IF
         END IF
C
C        Remove any trailing delimiters before looping.
C
         IF ( IRIGHT .GT. 0 ) ILEFT = MAX ( ILEFT, IRIGHT ) + 1
         IF ( ILEFT .GT. ISTOP ) RETURN
         CALL STRIPB( LINE(ILEFT:ISTOP), IL, ISTOP )
         ISTOP = ISTOP + ILEFT - 1
         ILEFT = IL + ILEFT - 1
         IF ( ILEFT .GT. ISTOP ) RETURN
         IF ( INDEX ( ',=', LINE(ILEFT:ILEFT) ) .NE. 0 ) 
     *      ILEFT = ILEFT + 1
         IF ( ILEFT .GT. ISTOP ) RETURN
 1       continue
c
c  The end of the VAX FORTRAN DO Loop that I commented out.  jrr.
c
c      END DO
 3000 FORMAT( E32.0 )
 4000 FORMAT( I32 )
C
      END
