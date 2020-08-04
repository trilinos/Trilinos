C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE CHKSYN (NINENT, IENTYP, CENTRY, RENTRY, NAMES,
     &   MAXENT, NUMENT, NAMENT, TYPENT, INXENT, VALENT, ITMENT,
     &   NERR)
C=======================================================================

C   --*** CHKSYN *** (ALGEBRA) Assign equation entries and checks syntax
C   --   Written by Amy Gilkey - revised 05/23/88
C   --
C   --CHKSYN assigns the equation entries from the input fields.
C   --It assigns NAMENT, TYPENT, INXENT, VALENT, and ITMENT for each entry
C   --(INXENT for a variable is assigned later).  Entries may be added
C   --or deleted from the equation.  The equation syntax is also checked.
C   --
C   --Aliased variables are expanded in this routine.
C   --
C   --Parameters:
C   --   NINENT - IN - the number of input entries
C   --   IENTYP - IN - the free-field reader entry type
C   --            0    = name string
C   --            1,2  = number
C   --            3    = character or **
C   --   CENTRY - IN - character string fields
C   --   RENTRY - IN - numeric fields
C   --   NAMES  - IN - the global, nodal, and element variable names
C   --   MAXENT - IN - the maximum number of entries + 1
C   --   NUMENT - OUT - the number of entries (as in /ENT../)
C   --   NAMENT - OUT - the equation entries (as in /ENT../)
C   --   TYPENT - OUT - the type of each equation entry (as in /ENT../)
C   --   INXENT - OUT - based on TYPENT (as in /ENT../)
C   --   VALENT - OUT - based on TYPENT (as in /ENT../)
C   --   ITMENT - OUT - based on TYPENT (as in /ENT../)
C   --   NERR   - I/O - the number of errors in the equation, may be set
C   --
C   --Common Variables:
C   --   Uses NUMALI, NAMALI, NIXALI, IXALI of /ALIAS../
C   --   Uses FNCNAM of /FNCTB./
C   --   Sets (by DATA) FNCNAM, FNCTYP, FNCSTO of /FNCTB./

      PARAMETER (ICURTM = 1, ILSTTM = 2, IONETM = 3)

      include 'exodusII.inc'
      include 'ag_namlen.blk'

      include 'ag_alias.blk'
      include 'ag_fnctbc.blk'

      INTEGER       IENTYP(*)
      CHARACTER*(*) CENTRY(*)
      REAL          RENTRY(*)
      CHARACTER*(*) NAMES(*)
      CHARACTER*(*) NAMENT(*)
      CHARACTER TYPENT(*)
      INTEGER INXENT(*)
      REAL VALENT(*)
      INTEGER ITMENT(*)

      LOGICAL NOTDON
      CHARACTER*5 STRA

      CHARACTER*(MAXNAM+20) ENT20
C      --ENT20 - long equation entry, before trucation

      CHARACTER*(mxstln) ENT, LSTENT
C      --ENT, LSTENT - single equation entries
      CHARACTER LSTTYP
C      --LSTTYP - single equation entry types

      CHARACTER*10 OPSYM
      SAVE OPSYM
C      --OPSYM - the ordered operators, one to a character, order is used
C      --   to determine the operator INXENT, space is reserved for
C      --   end of line operator (# is 0) and unary minus (~ is 1)

      INTEGER NPARM(MAXFNC)
      SAVE NPARM
C      --NPARM(i) - the required number of parameters for the function i

      CHARACTER*(mxstln) FNCOLD(10), FNCNEW(10)
      INTEGER IXFNC(10)
      SAVE FNCOLD, FNCNEW, IXFNC
C      --FNCOLD - old function names which have been replaced
C      --FNCNEW - the corresponding new function names (in FNCNAM)
C      --IXFNC - the index of the new function names in FNCTB.

      LOGICAL FIRST
      SAVE FIRST
      SAVE IADD
C      --FIRST - true only on the first time through routine
      DATA FIRST / .TRUE. /

      DATA OPSYM / ' +-*/^()=,' /
C      --Note first space to reserve index 1 for unary minus

      call infunc(nparm, fncold, fncnew)

      IF (FIRST) THEN
        IADD = 0
        NUMFNO = 0
 100    CONTINUE
        IF (FNCOLD(NUMFNO+1) .NE. ' ') THEN
          NUMFNO = NUMFNO + 1
          IXFNC(NUMFNO) = LOCSTR (FNCNEW(NUMFNO), NUMFNC, FNCNAM)
          GOTO 100
        END IF

        FIRST = .FALSE.
      END IF

C   --Assign TYPENT Constant, Operator, Variable or Function and
C   --set INXENT, VALENT, and ITMENT

      NUMENT = NINENT
      IF (NUMENT .GE. MAXENT) THEN
         CALL PRTERR ('CMDSPEC', 'Too many equation entries')
         NERR = NERR + 1
         NUMENT = MAXENT-1
      END IF

      DO 110 NENT = 1, NUMENT

         ENT20 = CENTRY(NENT)

         IF ((IENTYP(NENT) .EQ. 1) .OR. (IENTYP(NENT) .EQ. 2)) THEN

C         --Constant - store value

            NAMENT(NENT) = ENT20(1:MAXNAM)
            TYPENT(NENT) = 'C'
            INXENT(NENT) = -999
            VALENT(NENT) = RENTRY(NENT)
            ITMENT(NENT) = -999

         ELSE IF (IENTYP(NENT) .EQ. 3) THEN

C         --Operator - store index

            NAMENT(NENT) = ENT20(1:1)
            TYPENT(NENT) = 'O'
            ISYM = INDEX (OPSYM, ENT20(1:1))
            INXENT(NENT) = ISYM
            VALENT(NENT) = -999
            ITMENT(NENT) = -999
            IF (ISYM .EQ. 0) THEN
               NERR = NERR + 1
               CALL PRTERR ('CMDSPEC',
     &            '"' // ENT20(1:1) // '" is not a valid operator')
            END IF

         ELSE

C         --Variable (may be function) - store name and node/element
C         --specifier and time specifier

            IS = INDEX (ENT20, '$')
            IF (IS .GT. 0) THEN
               CALL CNVNUM (ENT20(IS+1:), NENUM, IERR)
               IF ((IERR .NE. 0) .OR. (NENUM .LE. 0)) THEN
                  NERR = NERR + 1
                  CALL PRTERR ('CMDSPEC', ENT20(:LENSTR(ENT20))
     &               // ' has an invalid node/element specifier')
                  NENUM = 1
               END IF
               ENT20(IS:) = ' '
            ELSE
               NENUM = 0
            END IF
            IS = INDEX (ENT20, ':')
            IF (IS .GT. 0) THEN
               ENT = ENT20(IS+1:)
               IF (ENT .EQ. ' ') THEN
                  ITIME = ILSTTM
               ELSE IF (ENT .EQ. '1') THEN
                  ITIME = IONETM
               ELSE
                  NERR = NERR + 1
                  CALL PRTERR ('CMDSPEC', ENT20(:LENSTR(ENT20))
     &               // ' has an invalid time specifier')
                  ITIME = ICURTM
               END IF
               ENT20(IS:) = ' '
            ELSE
               ITIME = ICURTM
            END IF

            IF (ENT20(MAXNAM+1:) .NE. ' ')
     &         CALL PRTERR ('CMDSPEC',
     &         'WARNING - ' // ENT20(:LENSTR(ENT20))
     &         // ' is truncated to "' // ENT20(1:MAXNAM) // '"')

            NAMENT(NENT) = ENT20(1:MAXNAM)
            TYPENT(NENT) = 'V'
            INXENT(NENT) = -999
            VALENT(NENT) = NENUM
            ITMENT(NENT) = ITIME
         END IF

  110 CONTINUE

C   --Temporarily add a '#' to mark the end of the equation to eliminate
C   --special end checks

      NAMENT(NUMENT+1) = 'line end'
      TYPENT(NUMENT+1) = 'O'

C   --Distinguish between variables and functions (functions names are
C   --followed by "(")

      DO 120 NENT = 1, NUMENT

         IF (TYPENT(NENT) .EQ. 'V') THEN
            IF (NAMENT(NENT+1) .EQ. '(') THEN

               ENT = NAMENT(NENT)

C            --Locate the function name in the function table

               INXF = LOCSTR (ENT, NUMFNC, FNCNAM)
               IF (INXF .LE. 0) THEN
                  IX = LOCSTR (ENT, NUMFNO, FNCOLD)
                  IF (IX .GT. 0) THEN
                     INXF = IXFNC(IX)
                     CALL PRTERR ('CMDSPEC',
     &                  'WARNING - Please use function name "'
     &                  // FNCNAM(INXF)(:LENSTR(FNCNAM(INXF)))
     &                  // '" instead of "' // ENT(:LENSTR(ENT)) // '"')
                     ENT = FNCNAM(INXF)
                  END IF
               END IF

               IF (INXF .GT. 0) THEN

                  IF ((VALENT(NENT) .NE. 0)
     &               .OR. (ITMENT(NENT) .NE. ICURTM)) THEN
                     NERR = NERR + 1
                     CALL PRTERR ('CMDSPEC', 'Function "'
     &                  // ENT(:LENSTR(ENT)) // '" has a specifier')
                  END IF

C               --Function - store name and index and number of parameters

                  NAMENT(NENT) = ENT
                  TYPENT(NENT) = 'F'
                  INXENT(NENT) = INXF
                  VALENT(NENT) = NPARM(INXF)
                  ITMENT(NENT) = -999

               ELSE
                  NERR = NERR + 1
                  CALL PRTERR ('CMDSPEC',
     &               'Function "' // ENT(:LENSTR(ENT)) //
     &               '" is undefined or a variable is followed by "("')
               END IF
            END IF
         END IF
  120 CONTINUE

C   --Find unary '+' or '-' and convert to constant or unary operator

      NENT = 1
  130 CONTINUE
      IF (NENT .LT. NUMENT) THEN
         NENT = NENT + 1

         ENT = NAMENT(NENT)

         IF ((ENT .EQ. '+') .OR. (ENT .EQ. '-')) THEN

            LSTENT = NAMENT(NENT-1)

            IF ((LSTENT .EQ. '=') .OR. (LSTENT .EQ. '(')
     &         .OR. (LSTENT .EQ. ',')) THEN
C            --A unary '+' or '-' is found

               IF ((ENT .EQ. '+') .OR. (TYPENT(NENT+1) .EQ. 'C')) THEN

C               --If it is a '+', delete the sign
C               --If it is a negative constant, delete the sign and alter
C               --the constant

                  DO 140 J = NENT, NUMENT+1-1
                     NAMENT(J) = NAMENT(J+1)
                     TYPENT(J) = TYPENT(J+1)
                     INXENT(J) = INXENT(J+1)
                     VALENT(J) = VALENT(J+1)
                     ITMENT(J) = ITMENT(J+1)
  140             CONTINUE
                  NUMENT = NUMENT - 1

                  IF (ENT .EQ. '-') VALENT(NENT) = - VALENT(NENT)
                  NENT = NENT - 1

               ELSE

C               --It is a signed expression, change to a unary minus

                  NAMENT(NENT) = '~'
                  INXENT(NENT) = 1

               END IF
            END IF
         END IF
         GOTO 130
      END IF

C   --Expand variable aliases

      NENT = 0
  150 CONTINUE
      IF (NENT .LT. NUMENT) THEN
         NENT = NENT + 1

         IF (TYPENT(NENT) .EQ. 'V') THEN

            IALI = LOCSTR (NAMENT(NENT), NUMALI, NAMALI)

            IF (IALI .GE. 1) THEN
               NALI = NIXALI(IALI)

C            --Make room for alias expansion

               IF (NUMENT+IADD .GE. MAXENT) THEN
                  CALL PRTERR ('CMDSPEC', 'Too many equation entries')
                  NERR = NERR + 1
                  NUMENT = MAXENT-1 - IADD
               END IF
               IADD = 2 * (NALI-1)
               DO 160 J = NUMENT+1, NENT, -1
                  NAMENT(J+IADD) = NAMENT(J)
                  TYPENT(J+IADD) = TYPENT(J)
                  INXENT(J+IADD) = INXENT(J)
                  VALENT(J+IADD) = VALENT(J)
                  ITMENT(J+IADD) = ITMENT(J)
  160          CONTINUE
               NUMENT = NUMENT + IADD

C            --Expand alias to "var, var, ..."

               NAMENT(NENT) = NAMES(IXALI(1,IALI))
               IX = NENT
               ISYM = INDEX (OPSYM, ',')

               DO 170 I = 2, NALI
                  NENT = NENT + 1
C               --Mark comma as alias comma
                  NAMENT(NENT) = ',.ALIAS.'
                  TYPENT(NENT) = 'O'
                  INXENT(NENT) = ISYM
                  VALENT(NENT) = -999
                  ITMENT(NENT) = -999
                  NENT = NENT + 1
                  NAMENT(NENT) = NAMES(IXALI(I,IALI))
                  TYPENT(NENT) = TYPENT(IX)
                  INXENT(NENT) = INXENT(IX)
                  VALENT(NENT) = VALENT(IX)
                  ITMENT(NENT) = ITMENT(IX)
  170          CONTINUE
            END IF
         END IF
         GOTO 150
      END IF

C   --Check that the first entry is a variable

      IF (TYPENT(1) .NE. 'V') THEN
         NERR = NERR + 1
         CALL PRTERR ('CMDSPEC',
     &      'The assigned variable is not specified')
      END IF

C   --Check for proper '='

      NES = 0
C      --NES - the number of equal signs

      DO 180 NENT = 1, NUMENT
         IF (TYPENT(NENT) .EQ. 'O') THEN
            IF (NAMENT(NENT) .EQ. '=') NES = NES + 1
         END IF
  180 CONTINUE

      IF (NES .EQ. 1) THEN
         IF (NAMENT(2) .NE. '=') THEN
            NERR = NERR + 1
            CALL PRTERR ('CMDSPEC', 'The "=" is misplaced')
         END IF
      ELSE IF (NES .LT. 1) THEN
         NERR = NERR + 1
         CALL PRTERR ('CMDSPEC', 'The "=" is missing')
      ELSE IF (NES .GT. 1) THEN
         NERR = NERR + 1
         CALL PRTERR ('CMDSPEC', 'Extra "=" found')
      END IF

C   --Check for matched parenthesis, and mark the "level" on parenthesis
C   --and commas

      NPERR = 0
C      --NPERR - set if too many right parenthesis at any time
      NP = 0
C      --NP - the number of parenthesis outstanding (+left, -right)

      DO 190 NENT = 1, NUMENT
         IF (TYPENT(NENT) .EQ. 'O') THEN
            ENT = NAMENT(NENT)(1:1)
            IF (ENT .EQ. '(') THEN
               NP = NP + 1
               ITMENT(NENT) = NP
            ELSE IF (ENT .EQ. ')') THEN
               ITMENT(NENT) = NP
               NP = NP - 1
               IF (NP .LT. 0) NPERR = NPERR + 1
            ELSE IF (ENT .EQ. ',') THEN
               ITMENT(NENT) = NP
C            --Mark comma
               IF (NAMENT(NENT) .EQ. ',') NAMENT(NENT) = ',.MARK.'
            END IF
         END IF
  190 CONTINUE

      IF ((NP .NE. 0) .OR. (NPERR .GT. 0)) THEN
         NERR = NERR + 1
         CALL PRTERR ('CMDSPEC', 'Parenthesis do not balance')
      END IF

C   --Check the number of parameters for all functions

      DO 210 NENT = 1, NUMENT

         IF (TYPENT(NENT) .EQ. 'F') THEN

            IE = NENT + 1
            IF (NAMENT(IE) .EQ. '(') THEN
               NPFND = 1
C               --NPFND - the number of parameters given for the function
               NLEVEL = ITMENT(IE)
C               --NLEVEL - the parenthesis level for the function
               NOTDON = .TRUE.

  200          CONTINUE
               IF ((IE .LT. NUMENT) .AND. NOTDON) THEN
                  IE = IE + 1
                  IF (TYPENT(IE) .EQ. 'O') THEN
                     ENT = NAMENT(IE)(1:1)
                     IF (ENT .EQ. ')') THEN
                        IF (ITMENT(IE) .EQ. NLEVEL) NOTDON = .FALSE.
                     ELSE IF (ENT .EQ. ',') THEN
                        IF (ITMENT(IE) .EQ. NLEVEL) THEN
C                        --Parameter divider, count it and mark it
                           NPFND = NPFND + 1
                           NAMENT(IE) = ','
                        END IF
                     END IF
                  END IF
                  GOTO 200
               END IF

            ELSE
               NPFND = 0
            END IF

C         --Set number of parameters for function if may vary
            IF (VALENT(NENT) .LT. 0) VALENT(NENT) = NPFND

            IF (NPFND .NE. VALENT(NENT)) THEN
               NERR = NERR + 1
               ENT = NAMENT(NENT)
               NPEXP = VALENT(NENT)
               CALL INTSTR (1, 0, NPEXP, STRA, LSTRA)
               CALL PRTERR ('CMDSPEC', 'Expected ' // STRA(:LSTRA)
     &            // ' parameter(s) for function ' // ENT(:LENSTR(ENT)))
               VALENT(NENT) = NPFND
            END IF
         END IF

  210 CONTINUE

C   --Check for non-parameter comma and reset scratch "level" and name
C   --on operators

      DO 230 NENT = 1, NUMENT

         IF (TYPENT(NENT) .EQ. 'O') THEN

            ITMENT(NENT) = -999

            IF (NAMENT(NENT) .EQ. ',.ALIAS.') THEN
               NERR = NERR + 1
               ENT = NAMENT(NENT-1)
               CALL PRTERR ('CMDSPEC',
     &            'Invalid alias expansion starting with "'
     &            // ENT(:LENSTR(ENT)) // '"')

C            --Wipe out commas for this alias to prevent repeated errors
               IE = NENT
  220          CONTINUE
               IF (NAMENT(IE) .EQ. ',.ALIAS.') THEN
                  NAMENT(IE) = ','
                  IE = IE + 2
                  GOTO 220
               END IF

            ELSE IF (NAMENT(NENT) .EQ. ',.MARK.') THEN
               NERR = NERR + 1
               ENT = NAMENT(NENT+1)
               CALL PRTERR ('CMDSPEC', 'Invalid "," before "'
     &            // ENT(:LENSTR(ENT)) // '"')
               NAMENT(NENT) = ','
            END IF
         END IF

  230 CONTINUE

C   --Check for consecutive variables or operators, etc.

      DO 240 NENT = 2, NUMENT

         LSTENT = NAMENT(NENT-1)
         ENT = NAMENT(NENT)
         LSTTYP = TYPENT(NENT-1)
C      --Set up a special type for ")" to prevent special checks
         IF ((LSTTYP .EQ. 'O') .AND. (LSTENT .EQ. ')')) LSTTYP = ')'

         IF (TYPENT(NENT) .EQ. 'O') THEN
            IF (ENT .EQ. '(') THEN
               IF ((LSTTYP .EQ. 'C') .OR. (LSTTYP .EQ. ')')) THEN
                  NERR = NERR + 1
                  CALL PRTERR ('CMDSPEC',
     &               'Consecutive "' // LSTENT(:LENSTR(LSTENT))
     &               // '" and "' // ENT(:LENSTR(ENT)) // '" found')
               END IF
            ELSE IF (ENT .NE. '~') THEN
               IF (LSTTYP .EQ. 'O') THEN
                  NERR = NERR + 1
                  CALL PRTERR ('CMDSPEC',
     &               'Consecutive "' // LSTENT(:LENSTR(LSTENT))
     &               // '" and "' // ENT(:LENSTR(ENT)) // '" found')
               END IF
            END IF

         ELSE
            IF (LSTTYP .NE. 'O') THEN
               NERR = NERR + 1
               CALL PRTERR ('CMDSPEC',
     &            'Consecutive "' // LSTENT(:LENSTR(LSTENT))
     &            // '" and "' // ENT(:LENSTR(ENT)) // '" found')
            END IF
         END IF
  240 CONTINUE

      IF ((TYPENT(NUMENT) .EQ. 'O')
     &   .AND. (NAMENT(NUMENT) .NE. ')')) THEN
         NERR = NERR + 1
         ENT = NAMENT(NUMENT)
         CALL PRTERR ('CMDSPEC',
     &      'Ending "' // ENT(:LENSTR(ENT)) // '" found')
      END IF

      RETURN
      END

      subroutine infunc(nparm, fncold, fncnew)

C     This subroutine was created in order to initialize parameters
C     that are in named common blocks.  The reason these variables
C     are set in the most labor intensive manner is that the
C     current fortran compiler on the SGI does not seem to recognize
C     BLOCK DATA files.

      include 'exodusII.inc'
      include 'ag_fnctbc.blk'

C     PARAMETER (MAXFNC=50)
C     PARAMETER (NUMFNC=39)
C      --MAXFNC - the number of defined functions
C     COMMON /FNCTBC/ FNCNAM(MAXFNC), FNCTYP(MAXFNC)
C     COMMON /FNCTBL/ FNCSTO(MAXFNC)
C     CHARACTER*8 FNCNAM
C     CHARACTER FNCTYP
C     LOGICAL FNCSTO
C      --   Assigned (by DATA) in CHKALG
C      --FNCNAM - the ordered function names, order is used to determine
C      --   the function INXENT
C      --FNCTYP(i) - the type of function i; ' ' for type determined by
C      --   parameters
C      --FNCSTO(i) - true iff function i needs storage (a time fC)
C      --FNCSTO(i) - true iff function i needs storage (a time function)

      integer nparm(numfnc)
      character*(MXSTLN) fncold(10), fncnew(10)

C                            1         2         3
C                   12345678901234567890123456789012
      FNCNAM(1)  = 'AINT                            '
      FNCNAM(2)  = 'ANINT                           '
      FNCNAM(3)  = 'ABS                             '
      FNCNAM(4)  = 'MOD                             '
      FNCNAM(5)  = 'SIGN                            '
      FNCNAM(6)  = 'DIM                             '
      FNCNAM(7)  = 'MAX                             '
      FNCNAM(8)  = 'MIN                             '
      FNCNAM(9)  = 'SQRT                            '
      FNCNAM(10) = 'EXP                             '
      FNCNAM(11) = 'LOG                             '
      FNCNAM(12) = 'LOG10                           '
      FNCNAM(13) = 'SIN                             '
      FNCNAM(14) = 'COS                             '
      FNCNAM(15) = 'TAN                             '
      FNCNAM(16) = 'ASIN                            '
      FNCNAM(17) = 'ACOS                            '
      FNCNAM(18) = 'ATAN                            '
      FNCNAM(19) = 'ATAN2                           '
      FNCNAM(20) = 'SINH                            '
      FNCNAM(21) = 'COSH                            '
      FNCNAM(22) = 'TANH                            '
      FNCNAM(23) = 'TMAG                            '
      FNCNAM(24) = 'PMAX                            '
      FNCNAM(25) = 'PMIN                            '
      FNCNAM(26) = 'PMAX2                           '
      FNCNAM(27) = 'PMIN2                           '
      FNCNAM(28) = 'IFLZ                            '
      FNCNAM(29) = 'IFEZ                            '
      FNCNAM(30) = 'IFGZ                            '
      FNCNAM(31) = 'SUM                             '
      FNCNAM(32) = 'SMAX                            '
      FNCNAM(33) = 'SMIN                            '
      FNCNAM(34) = 'ENVMAX                          '
      FNCNAM(35) = 'ENVMIN                          '
      FNCNAM(36) = 'UHIST                           '
      FNCNAM(37) = 'UGLOB                           '
      FNCNAM(38) = 'UNODE                           '
      FNCNAM(39) = 'UELEM                           '

      NPARM(1)  = 1
      NPARM(2)  = 1
      NPARM(3)  = 1
      NPARM(4)  = 2
      NPARM(5)  = 2
      NPARM(6)  = 2
      NPARM(7)  = -1
      NPARM(8)  = -1
      NPARM(9)  = 1
      NPARM(10) = 1
      NPARM(11) = 1
      NPARM(12) = 1
      NPARM(13) = 1
      NPARM(14) = 1
      NPARM(15) = 1
      NPARM(16) = 1
      NPARM(17) = 1
      NPARM(18) = 1
      NPARM(19) = 2
      NPARM(20) = 1
      NPARM(21) = 1
      NPARM(22) = 1
      NPARM(23) = 6
      NPARM(24) = 6
      NPARM(25) = 6
      NPARM(26) = 3
      NPARM(27) = 3
      NPARM(28) = 3
      NPARM(29) = 3
      NPARM(30) = 3
      NPARM(31) = 1
      NPARM(32) = 1
      NPARM(33) = 1
      NPARM(34) = 1
      NPARM(35) = 1
      NPARM(36) = -1
      NPARM(37) = -1
      NPARM(38) = -1
      NPARM(39) = -1

      do 100 i = 1, NUMFNC
         fnctyp(i) = ' '
  100 continue
      FNCTYP(31) = 'G'
      FNCTYP(32) = 'G'
      FNCTYP(33) = 'G'
      FNCTYP(36) = 'H'
      FNCTYP(37) = 'G'
      FNCTYP(38) = 'N'
      FNCTYP(39) = 'E'

      do 110 i = 1, NUMFNC
         fncsto(i) = .FALSE.
  110 continue
      FNCSTO(32) = .TRUE.
      FNCSTO(33) = .TRUE.
      FNCSTO(34) = .TRUE.
      FNCSTO(35) = .TRUE.

      FNCOLD(1) = 'ATN2                            '
      FNCOLD(2) = 'PMX2                            '
      FNCOLD(3) = 'PMN2                            '
      FNCOLD(4) = '                                '

      FNCNEW(1) = 'ATAN2                           '
      FNCNEW(2) = 'PMAX2                           '
      FNCNEW(3) = 'PMIN2                           '
      FNCNEW(4) = '                                '

      return
      end

