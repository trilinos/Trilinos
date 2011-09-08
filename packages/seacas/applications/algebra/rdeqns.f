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
      SUBROUTINE RDEQNS (A, C, NAMECO, BLKTYP, NAMES, QAREC, INFREC,
     &  TIMES, IPTIMS, SELELB, IDELB, VISELB, KIEVOK, MAXSTK, NLOG,
     *  MERR)
C=======================================================================

C   --*** RDEQNS *** (ALGEBRA) Read and check the equations
C   --   Written by Amy Gilkey - revised 05/23/88
C   --   Modified for ExodusIIV2 format - 9/5/95
C   --
C   --RDEQNS processes the input equations as follows:
C   --   o Reads the line and parses it into fields.
C   --   o Checks the equation for syntax.
C   --   o Stores the variable names from the equations.
C   --   o Puts the equation in postfix form.
C   --   o Adds the equation variables to a global list.
C   --
C   --If any errors are found, that equation is ignored.
C   --After all the equations have been read in, all the variables are
C   --gathered into the /VAR../ arrays.  The assigned variables are
C   --after the expression variables.
C   --
C   --Parameters:
C   --   A      - IN/OUT - the dynamic memory base array
C   --   NAMECO     - IN - the coordinate names
C   --   BLKTYP     - IN - the element block names
C   --   NAMES      - IN - the database variable names
C   --   QAREC      - IN - the QA records containing:
C   --               (1) - the analysis code name
C   --               (2) - the analysis code QA descriptor
C   --               (3) - the analysis date
C   --               (4) - the analysis time
C   --   INFREC     - IN - the information records
C   --   TIMES      - IN - the database time steps
C   --   IPTIMS    - OUT - the selected times steps
C   --   SELELB    - OUT - the selected element blocks?
C   --   IDELB     - IN  - the element block IDs
C   --   VISELB(i) - OUT - true iff element block i is to be written
C   --   KIEVOK - IN/OUT - the dynamic memory index of ISEVOK;
C   --                     input variables on input, temporary and
C   --                     output added on output
C   --   MAXSTK    - OUT - the maximum stack size needed for any equation
C   --   MERR      - OUT - error flag
C   --
C   --Common Variables:
C   --   Sets and uses NUMEQN, NUMENT of /ENT../
C   --   Uses NAMENT of /ENT../
C   --   Sets NUMINP, IXLHS, NAMVAR, TYPVAR, IDVAR, ISTVAR, IEVVAR of /VAR../
C   --   Uses NSTEPS of /DBNUMS/
C   --   Sets NPTIMS, TMIN, TMAX, DELT, NINTV, WHONLY of /TIMES/
C   --   Sets ISZOOM of /ZOOM/
C   --   Sets EQNLIN of /EQNLNS/

      include 'params.blk'
      include 'namlen.blk'
      include 'numeqn.blk'
      include 'ent.blk'
      include 'var.blk'
      include 'times.blk'
      include 'dbnums.blk'
      include 'zoom.blk'
      include 'filter.blk'
      include 'remove.blk'
      include 'eqnlns.blk'

      PARAMETER (ICURTM = 1, ILSTTM = 2, IONETM = 3)

      DIMENSION    A(*)
      CHARACTER*1  C(1)
      CHARACTER*(namlen)  NAMECO(*)
      CHARACTER*(MXSTLN)  BLKTYP(*)
      CHARACTER*(namlen)  NAMES(*)
      CHARACTER*(MXSTLN)  QAREC(4,*)
      CHARACTER*(MXLNLN)  INFREC(*)
      REAL    TIMES(*)
      INTEGER IPTIMS(*)
      LOGICAL SELELB(NELBLK)
      INTEGER IDELB(*)
      LOGICAL VISELB(NELBLK)
      INTEGER KIEVOK
      INTEGER MAXSTK
      INTEGER MERR

      INTEGER      IENTYP(MAXENT+1)
      CHARACTER*(maxnam) CENTRY(MAXENT)
      INTEGER      IENTRY(MAXENT)
      REAL         RENTRY(MAXENT)

      CHARACTER*132 LINE
      CHARACTER*256 INLINE(5)
      CHARACTER*(MXSTLN) RETVRB
      CHARACTER*5 STRA
      LOGICAL PRTEQN
      LOGICAL SAVLOG
      LOGICAL ISEQN

C     Open the log file - temporary file unless the user decides to save it
      SAVLOG = .FALSE.

      MERR = 0

C      WHONLY = .FALSE.

C   --Initialize the selected times (0 = select all times)
      CALL INITIM (0, NSTEPS, TIMES, TMIN, TMAX,
     &             DELT, NINTV, NPTIMS, IPTIMS)

      ISZOOM   = .FALSE.
      ISFILTER = .FALSE.
      ISREMOVE = .FALSE.
      IRMCNT   = 0
      do i=1, 1024
        idsrem(i) = 0
      end do

C     Initializes all logicals in a list to a specified value
      CALL INILOG (NELBLK, .TRUE., VISELB)

      IDEFEV = 0


      CALL INILOG (NELBLK, .TRUE., SELELB)

C   --TIME is always an input variable

      NAMVAR(1) = 'TIME'
      TYPVAR(1) = 'T'
      IDVAR(1) = -999
      ISTVAR(ICURTM,1) = 1
      ISTVAR(ILSTTM,1) = 0
      ISTVAR(IONETM,1) = 0
      IEVVAR(1) = -999

      NUMINP = 1
      IXLHS = MAXVAR + 1

      NUMEQN = 0
      MAXSTK = 0
      NUMEV = NVAREL
      MAXEV = NUMEV

      WRITE (*, *)
      WRITE (*, 10000) 'Enter the equations'
      WRITE (*, *)
      PRTEQN = .FALSE.

  100 CONTINUE
      IF (.TRUE.) THEN

C      --Write out all lines input so far

         IF (PRTEQN) THEN
            WRITE (*, *)
            CALL SHOW ('TMIN', ' ',
     &         NAMECO, BLKTYP, NAMES,
     &         TIMES, IPTIMS, IDELB, VISELB, SELELB)
            CALL SHOW ('ZOOM', ' ',
     &         NAMECO, BLKTYP, NAMES,
     &         TIMES, IPTIMS, IDELB, VISELB, SELELB)
            CALL SHOW ('BLOCKS', ' ',
     &         NAMECO, BLKTYP, NAMES,
     &         TIMES, IPTIMS, IDELB, VISELB, SELELB)
            CALL SHOW ('SAVE', ' ',
     &         NAMECO, BLKTYP, NAMES,
     &         TIMES, IPTIMS, IDELB, VISELB, SELELB)
            CALL SHOW ('DELETE', ' ',
     &         NAMECO, BLKTYP, NAMES,
     &         TIMES, IPTIMS, IDELB, VISELB, SELELB)

            WRITE (*, *)
            DO 110 I = 1, NUMEQN
               WRITE (*, 10000) 'ALG> ', EQNLIN(I)(:LENSTR(EQNLIN(I)))
  110       CONTINUE

            PRTEQN = .FALSE.
         END IF

C      --Read in equations and commands

         CALL GETINP (0, 0, 'ALG> ', LINE, IOSTAT)
         
C      --Check for comment ''''
         ICMT = INDEX(LINE, '''')
         IF (ICMT .EQ. 1) THEN
C          ... Only comment on line
           GO TO 100
         ELSE IF (ICMT .GT. 1) THEN
           LINE(ICMT-1:) = ' '
         END IF

         EQNLIN(NUMEQN+1) = LINE
  120    CONTINUE
         I = INDEX (EQNLIN(NUMEQN+1), '>')
         IF (I .GT. 0) THEN
            CALL GETINP (0, 0, '   > ', LINE, IOSTAT)
            EQNLIN(NUMEQN+1)(I:) = LINE
            GOTO 120
         END IF
         IF (IOSTAT .LT. 0) THEN
            CALL PRTERR ('CMDSPEC', 'End of command stream before END')
            RETVRB = 'END'
            GOTO 180
         END IF

C      --Parse the equation or command

         ISEQN = (INDEX (EQNLIN(NUMEQN+1), '=') .GT. 0)

         IF (ISEQN) THEN
c           CALL APARSE(input_line, max_num_parsed_flds,
C              num_parsed_flds, field_type, alphnum_flds, numeric_fld)
            CALL APARSE (EQNLIN(NUMEQN+1), MAXENT,
     &         NENT, IENTYP, CENTRY, RENTRY)
         ELSE
            ICONT = 0
            CALL FFISTR (EQNLIN(NUMEQN+1), MAXENT, ICONT,
     &         NENT, IENTYP, CENTRY, IENTRY, RENTRY)
            IENTYP(MIN(NENT,MAXENT)+1) = -999
         END IF
         IF (NENT .EQ. 0) GOTO 100

         IF (.NOT. ISEQN) THEN
C         --Input was a command not an equation
C         --Process command line

            INLINE(1) = ' '
            DO 130 I = 2, 5
               INLINE(I) = CHAR(0)
  130       CONTINUE

C         --Process command
            CALL COMAND (A, INLINE, IENTYP, CENTRY, IENTRY, RENTRY,
     &        NAMECO, BLKTYP, NAMES, TIMES, IPTIMS, IDELB,
     &        VISELB, SELELB, QAREC, INFREC, RETVRB, MERR)

C         --Write command to log file
            IF ((NLOG .GT. 0) .AND. (INLINE(1) .NE. ' ')) THEN
               DO 140 I = 1, 5
                  IF (INLINE(I)(1:1) .EQ. CHAR(0)) GOTO 150
                  WRITE (NLOG, '(A)') INLINE(I)(:LENSTR(INLINE(I)))
  140          CONTINUE
  150          CONTINUE
            END IF

C         --Handle special commands

            PRTEQN = .FALSE.
            IF (RETVRB .EQ. 'PRINT') THEN
               PRTEQN = .TRUE.

            ELSE IF (RETVRB .EQ. 'END') THEN
               GOTO 170

            ELSE IF (RETVRB .EQ. 'QUIT') THEN
               GOTO 170

            ELSE IF (RETVRB .EQ. 'LOG') THEN
               IF (NLOG .LE. 0) THEN
                  CALL PRTERR ('CMDERR', 'Log file cannot be opened')
               ELSE
                  SAVLOG = .TRUE.
                  WRITE (*, 10000) 'Log file will be saved'
               END IF

            ELSE IF (RETVRB .EQ. 'BLOCKS') THEN
C            --unable to move to main program
C            --Expand the ISEVOK array for the new default selection
               NUMEV = NUMEV + 1
               IF (NUMEV .GT. MAXEV) THEN
                  MAXEV = MAXEV + 50
c                  CALL MDGET (NELBLK * MAXEV)
                  CALL MDLONG ('ISEVOK', KIEVOK, NELBLK * MAXEV)
                  CALL MDSTAT (NERR, MEM)
                  IF (NERR .GT. 0) THEN
                     CALL MEMERR
                     MERR = 1
                     RETURN
                  END IF
               END IF

C            --Put default selection in ISEVOK array
               IDEFEV = NUMEV
               K = NELBLK * (IDEFEV-1)
               CALL CPYLOG (NELBLK, SELELB, A(KIEVOK+K))

C            --Make this command a dummy equation
               IF (NUMEQN .GE. MAXEQN) THEN
                  CALL INTSTR (1, 0, MAXEQN, STRA, LSTRA)
                  CALL PRTERR ('CMDSPEC', 'Only ' // STRA(:LSTRA)
     &               // ' equations will be accepted')
                  GOTO 170
               END IF
               NUMEQN = NUMEQN + 1
               NUMENT(NUMEQN) = 0
            END IF

         ELSE

C         --Process equation

            IF (NUMEQN .GE. MAXEQN) THEN
               CALL INTSTR (1, 0, MAXEQN, STRA, LSTRA)
               CALL PRTERR ('CMDSPEC', 'Only ' // STRA(:LSTRA)
     &            // ' equations will be accepted')
               GOTO 170
            END IF
C           Equation is accepted - increment number of equations
            NUMEQN = NUMEQN + 1

C           Initialize number of equation errors
            NEQERR = 0

C         --Assign the fields to the equation entries and check syntax
C           CALL CHKSYN(num_parsed_flds, fld_type, alphanum_fld, real,fld,
C                       GNEnames, max_ent_line, nument_eqn, eqn_entries,
C                       type_each_eqn_entry, inxent_typent, valent_typent,
C                       itment_typent, num_error_eqn)
            CALL CHKSYN (NENT, IENTYP, CENTRY, RENTRY, NAMES, MAXENT,
     &           NUMENT(NUMEQN), NAMENT(1,NUMEQN), TYPENT(1,NUMEQN),
     &           INXENT(1,NUMEQN), VALENT(1,NUMEQN), ITMENT(1,NUMEQN),
     &           NEQERR)

C         --Put the equation into postfix form
            CALL POSTFX (
     &           NUMENT(NUMEQN), NAMENT(1,NUMEQN), TYPENT(1,NUMEQN),
     &           INXENT(1,NUMEQN), VALENT(1,NUMEQN), ITMENT(1,NUMEQN),
     &           NEQERR)

C         --unable to move to main program
C         --Expand the ISEVOK array for temporary and output variables
            IF (NUMEV + 20 .GT. MAXEV) THEN
               MAXEV = MAXEV + 50
c               CALL MDGET (NELBLK * MAXEV)
               CALL MDLONG ('ISEVOK', KIEVOK, NELBLK * MAXEV)
               CALL MDSTAT (NERR, MEM)
               IF (NERR .GT. 0) THEN
                  CALL MEMERR
                  MERR = 1
                  RETURN
               END IF
            END IF

C         --Determine if the variables are in the database and where
C         --and enter the new variables into the /VAR../ arrays
            CALL LOCEQV (NAMECO, NAMES,
     &         NUMENT(NUMEQN), NAMENT(1,NUMEQN), TYPENT(1,NUMEQN),
     &         INXENT(1,NUMEQN), VALENT(1,NUMEQN),
     &         ITMENT(1,NUMEQN), IEVENT(1,NUMEQN), VSZENT(1,NUMEQN),
     &         A(KIEVOK), IDEFEV, NUMEV, NSTK, NEQERR)

C         --Check for syntax errors
            IF (NEQERR .GT. 0) THEN
               CALL PRTERR ('CMDSPEC', 'Equation ignored')
               NUMEQN = NUMEQN - 1
               GOTO 160
            END IF

C         --Handle special case of X = Y
            IF (IDEFEV .LE. 0) THEN
               IF ((NUMENT(NUMEQN) .EQ. 3) .AND.
     &            (NAMENT(1,NUMEQN) .EQ. NAMENT(3,NUMEQN))) THEN
                  NUMENT(NUMEQN) = 0
               END IF
            END IF

C         --Adjust the maximum stack size
            MAXSTK = MAX (MAXSTK, NSTK)

C         --Write equation to log file
            IF (NLOG .GT. 0) THEN
              WRITE (NLOG,'(A)') EQNLIN(NUMEQN)(:LENSTR(EQNLIN(NUMEQN)))
            END IF
         END IF

  160    CONTINUE
         GOTO 100
      END IF

  170 CONTINUE

C   --Scan lines after END
      CALL SCNEOF

  180 CONTINUE

      IF (SAVLOG) THEN
        IF (NLOG .GT. 0) THEN
          CLOSE (NLOG, IOSTAT=IDUM)
        END IF
      ELSE
        IF (NLOG .GT. 0) THEN
          CLOSE (NLOG, STATUS='DELETE', IOSTAT=IDUM)
        END IF
      END IF
          

C   --Adjust the ISEVOK array length
c      CALL MDGET (NELBLK * MAXEV)
      CALL MDLONG ('ISEVOK', KIEVOK, NELBLK * MAXEV)
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) THEN
         CALL MEMERR
         MERR = 1
         RETURN
      END IF



  190 CONTINUE
C     Simulate error so that no output file is generated
      IF (RETVRB .EQ. 'QUIT') MERR =  1
      RETURN

10000  FORMAT (1X, 5A)
      END
