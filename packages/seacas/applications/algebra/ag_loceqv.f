C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE LOCEQV (NAMECO, NAMES,
     &   NUMENT, NAMENT, TYPENT, INXENT, VALENT, ITMENT, IEVENT, VSZENT,
     &   ISEVOK, IDEFEV, NUMEV, MAXSTK, NERR)
C=======================================================================

C   --*** LOCEQV *** (ALGEBRA) Locate equation variables and add to list
C   --   Written by Amy Gilkey - revised 05/25/88
C   --
C   --LOCEQV determines the type, location, and node/element specific (if any)
C   --of an equation's variables.  The variable names, types, and
C   --database indices (for input variables) are added to the /VAR../
C   --arrays.
C   --
C   --The variable name in the NAMVAR array is stripped of any specifier
C   --to eliminate multiple references to a single variable.
C   --The node/element specifier is stored in VALENT and associated with
C   --the equation reference.  Similarly, all global input variables are
C   --assigned one dummy name, and the variable number is stored in VALENT.
C   --
C   --If there have been no errors up to this point, the equation is in
C   --postfix form.  Each entry in the equation is assigned a number
C   --representing the number of items expected after the entry is put
C   --on the stack.  With a nodal or element variable (or an action
C   --involving such a variable), this is simply the number of nodes or
C   --elements.  With a constant, TIME or global variable, this may be
C   --one or the number of nodes or element if the value must be converted
C   --to an array.
C   --
C   --Parameters:
C   --   NAMECO - IN - the coordinate names
C   --   NAMES - IN - the global, nodal, and element variable names
C   --   NUMENT - IN - the number of entries (as in /ENT../)
C   --   NAMENT - IN/OUT - the equation entries (as in /ENT../)
C   --   TYPENT - IN - the type of each equation entry (as in /ENT../)
C   --   INXENT - IN/OUT - based on TYPENT (as in /ENT../)
C   --   VALENT - IN/OUT - based on TYPENT (as in /ENT../)
C   --   ITMENT - IN/OUT - based on TYPENT (as in /ENT../)
C   --   IEVENT - IN/OUT - the ISEVOK index for this entry (as in /ENT../)
C   --   VSZENT - OUT - the type of variable on the stack after processing
C   --      this entry (as in /ENT../)
C   --   ISEVOK - IN/OUT - the element variable truth table;
C   --      includes input, temporary, and output variables
C   --   IDEFEV - IN - the default ISEVOK variable index
C   --   NUMEV - IN/OUT - the maximum ISEVOK variable index
C   --   MAXSTK - OUT - the maximum stack size for this equation
C   --   NERR - IN/OUT - the number of errors in the equation, may be set
C   --
C   --Common Variables:
C   --   Sets NUMINP, IXLHS, NAMVAR, TYPVAR, IDVAR, ISTVAR, IEVVAR of /VAR../
C   --   Uses NDIM, NUMNP, NUMEL, NVARNP, NVAREL, NVARGL of /DBNUMS/
C   --   Uses FNCTYP of /FNCTB./

      include 'exodusII.inc'
      include 'ag_namlen.blk'
      include 'ag_numeqn.blk'
      include 'ag_var.blk'
      include 'ag_dbnums.blk'
      include 'ag_fnctbc.blk'
      PARAMETER (ICURTM = 1, ILSTTM = 2, IONETM = 3)

      CHARACTER*(*) NAMECO(*)
      CHARACTER*(*) NAMES(*)
      CHARACTER*(*) NAMENT(*)
      CHARACTER TYPENT(*)
      INTEGER INXENT(*)
      REAL VALENT(*)
      INTEGER ITMENT(*)
      INTEGER IEVENT(*)
      CHARACTER VSZENT(*)
      LOGICAL ISEVOK(NELBLK,*)

      CHARACTER*(maxnam) NAME
      CHARACTER TYPEQV
      CHARACTER VTYP
      LOGICAL ALLSAM
      CHARACTER*5 STRA, STRB

      CHARACTER TYPSTK(MAXENT)
      INTEGER INXSTK(MAXENT)
C      --TYPSTK - type of variables for stack item
C      --INXSTK - equation entry which pushed this stack item

C   --Save the /VAR../ indices so they can be restored in case of error
      NINP = NUMINP
      ILHS = IXLHS

C   --Set up to check postfix stack if no errors so far
      ITOS = 0
      INERR = NERR
      IF (INERR .NE. 0) TYPSTK(1) = 'G'

      MAXSTK = 0

      VTYP = ' '
      DO 100 NENT = 1, NUMENT
         IEVENT(NENT) = -999
  100 CONTINUE

      DO 180 NENT = 3, NUMENT

         IF (TYPENT(NENT) .EQ. 'V') THEN

C         --Get the name and the node/element specifier, if any

            NAME = NAMENT(NENT)
            NENUM = INT(VALENT(NENT))

C         --Locate where the variable name is in the list of variables

            CALL LOCNAM (NAME, NAMECO, NAMES,
     &         NENUM, TYPEQV, IDEQV, NERR)

C         --Reference variable from entry

            VALENT(NENT) = NENUM
            IF (IDEQV .LT. 0) THEN
               INXENT(NENT) = -1
            ELSE
               INXENT(NENT) = +1
            END IF
C            --Set INXENT -1 for reference to LHS variable, +1 for
C            --input variable; used to link variables
            IF (TYPEQV .EQ. 'E') THEN
               NEWEV = NUMEV + 1
               IEVOLD = IDEFEV
               IF (IDEQV .GE. 0) THEN
                  IEVNEW = IDEQV
               ELSE
                  IEVNEW = IEVVAR(IABS(IDEQV))
               END IF
               CALL MAKEVO (NELBLK, ISEVOK,
     &            IEVOLD, IEVNEW, NEWEV, IEVOLD)
               NEWEV = IEVOLD
               IF (NEWEV .GT. NUMEV) NUMEV = NEWEV
               IEVENT(NENT) = NEWEV
            ELSE
               IEVENT(NENT) = -999
            END IF

C         --Correct variable name to eliminate duplicate references in
C         --the variable array

            IF (((TYPEQV .EQ. 'H') .OR. (TYPEQV .EQ. 'G'))
     &         .AND. (IDEQV .GE. 0)) THEN
C            --Input history and global variables are all in the same entry
               IF (TYPEQV .EQ. 'H') THEN
                  NAME = '.GLOBAL'
               ELSE
                  NAME = '.GLOBAL'
               END IF
               NAMENT(NENT) = NAME
            END IF

C         --Insert name in list of input variables, if not LHS and not
C         --already there; update time reference

            IF (TYPEQV .EQ. 'C') ITMENT(NENT) = ICURTM
            IF (IDEQV .GE. 0) THEN
               IVAR = LOCSTR (NAME, NUMINP, NAMVAR)
               IF (IVAR .LE. 0) THEN
                  NUMINP = NUMINP + 1
                  IVAR = NUMINP
                  IF (NUMINP .LT. IXLHS) THEN
                     NAMVAR(IVAR) = NAME
                     TYPVAR(IVAR) = TYPEQV
                     IDVAR(IVAR) = IDEQV
                     ISTVAR(ICURTM,IVAR) = 0
                     ISTVAR(ILSTTM,IVAR) = 0
                     ISTVAR(IONETM,IVAR) = 0
                     IF (TYPEQV .EQ. 'E') THEN
                        IEVVAR(IVAR) = IDEQV
                     ELSE
                        IEVVAR(IVAR) = -999
                     END IF
                  END IF
               END IF
               ITM = ITMENT(NENT)
               ISTVAR(ITM,IVAR) = 1
            ELSE
               ITM = ITMENT(NENT)
               IF (ITM .NE. ICURTM) THEN
                  NUMLHS = MAXVAR - IXLHS + 1
                  IVAR = LOCSTR (NAME, NUMLHS, NAMVAR(IXLHS))
                  IF (IVAR .GT. 0) THEN
                     IVAR = IXLHS + IVAR - 1
                     ISTVAR(ITM,IVAR) = 1
                  END IF
               END IF
            END IF

C         --Get the variable type to check arrays

            IF ((TYPEQV .EQ. ' ') .OR. (TYPEQV .EQ. 'T') .OR.
     *               (TYPEQV .EQ. 'G')) THEN
               VTYP = 'G'
             ELSE IF (TYPEQV .EQ. 'H') THEN
               VTYP = 'H'
            ELSE
               IF (NENUM .NE. 0) THEN
                  VTYP = 'G'
               ELSE IF (TYPEQV .EQ. 'C') THEN
                  VTYP = 'N'
               ELSE
                  VTYP = TYPEQV
               END IF
            END IF
         END IF

C      --Get the equation variable type, and check for mixed array errors

         IF (INERR .EQ. 0) THEN

C         --Check by manipulating the postfix string

            IF (TYPENT(NENT) .EQ. 'C') THEN

C            --Constant - push global type

               ITOS = ITOS + 1
               IF (ITOS .GT. MAXENT) GOTO 190
               MAXSTK = MAX (MAXSTK, ITOS)

               TYPSTK(ITOS) = 'G'
               INXSTK(ITOS) = NENT
               VSZENT(NENT) = 'G'

            ELSE IF (TYPENT(NENT) .EQ. 'V') THEN

C            --Variable - push variable array type

               ITOS = ITOS + 1
               IF (ITOS .GT. MAXENT) GOTO 190
               MAXSTK = MAX (MAXSTK, ITOS)

               TYPSTK(ITOS) = VTYP
               INXSTK(ITOS) = NENT
               VSZENT(NENT) = VTYP

            ELSE IF (TYPENT(NENT) .EQ. 'O') THEN

C            --Operator - pop operands, push type of any operand array

               IF (INXENT(NENT) .LE. 1) THEN
                  NOP = 1
               ELSE
                  NOP = 2
               END IF
               IOP1 = ITOS - NOP + 1
               ILAST = ITOS
               ITOS = IOP1
               IF (ITOS .LE. 0) GOTO 190

               ALLSAM = .TRUE.
               VTYP = TYPSTK(IOP1)
               DO 110 I = IOP1+1, ILAST
                  IF (VTYP .NE. TYPSTK(I)) ALLSAM = .FALSE.
                  IF (TYPSTK(I) .NE. VTYP) THEN
                     IF (TYPSTK(I) .EQ. 'H') THEN
                        CONTINUE
                     ELSE IF (TYPSTK(I) .EQ. 'G') THEN
                        IF ((VTYP .EQ. 'T') .OR. (VTYP .EQ. 'H'))
     &                     VTYP = TYPSTK(I)
                     ELSE IF ((TYPSTK(I) .EQ. 'N')
     &                  .OR. (TYPSTK(I) .EQ. 'E')) THEN
                        IF ((VTYP .EQ. 'N') .OR. (VTYP .EQ. 'E')) THEN
                           IF (TYPSTK(I) .NE. VTYP) VTYP = '?'
                        ELSE
                           IF (VTYP .NE. '?') VTYP = TYPSTK(I)
                        END IF
                     ELSE IF (TYPSTK(I) .EQ. '?') THEN
                        VTYP = TYPSTK(I)
                     END IF
                  END IF
  110          CONTINUE

               IF (VTYP .EQ. 'E') THEN
                  NEWEV = NUMEV + 1
                  IEVOLD = IDEFEV
                  DO 120 I = IOP1, ILAST
                     IF (TYPSTK(I) .EQ. 'E') THEN
                        IEVNEW = IEVENT(INXSTK(I))
                        CALL MAKEVO (NELBLK, ISEVOK,
     &                     IEVOLD, IEVNEW, NEWEV, IEVOLD)
                     END IF
  120             CONTINUE
                  NEWEV = IEVOLD
                  IF (NEWEV .GT. NUMEV) NUMEV = NEWEV
                  IEVENT(NENT) = NEWEV
               END IF

               IF (.NOT. ALLSAM) THEN
                  DO 130 I = IOP1, ILAST
                     IF (VTYP .NE. TYPSTK(I)) THEN
                        VSZENT(INXSTK(I)) = VTYP
                     END IF
  130             CONTINUE
               END IF

               TYPSTK(ITOS) = VTYP
               INXSTK(ITOS) = NENT
               VSZENT(NENT) = VTYP

            ELSE IF (TYPENT(NENT) .EQ. 'F') THEN

C            --Function - pop parameters and push type of function

               NPARM = INT(VALENT(NENT))
               IPARM = ITOS - NPARM + 1
               ILAST = ITOS
               ITOS = IPARM
               IF (ITOS .LE. 0) GOTO 190

               IF (FNCTYP(INXENT(NENT)) .EQ. ' ') THEN

C               --Get function type from parameter types

                  ALLSAM = .TRUE.
                  VTYP = TYPSTK(IPARM)
                  DO 140 I = IPARM+1, ILAST
                     IF (VTYP .NE. TYPSTK(I)) ALLSAM = .FALSE.
                     IF (TYPSTK(I) .NE. VTYP) THEN
                        IF (TYPSTK(I) .EQ. 'H') THEN
                           CONTINUE
                        ELSE IF (TYPSTK(I) .EQ. 'G') THEN
                           IF (VTYP .EQ. 'H') VTYP = TYPSTK(I)
                        ELSE IF ((TYPSTK(I) .EQ. 'N')
     &                     .OR. (TYPSTK(I) .EQ. 'E')) THEN
                           IF ((VTYP .EQ. 'N')
     &                        .OR. (VTYP .EQ. 'E')) THEN
                              IF (TYPSTK(I) .NE. VTYP) VTYP = '?'
                           ELSE
                              IF (VTYP .NE. '?') VTYP = TYPSTK(I)
                           END IF
                        ELSE IF (TYPSTK(I) .EQ. '?') THEN
                           VTYP = TYPSTK(I)
                        END IF
                     END IF
  140             CONTINUE

                  IF (.NOT. ALLSAM) THEN
                     DO 150 I = IPARM, ILAST
                        IF (VTYP .NE. TYPSTK(I)) THEN
                           VSZENT(INXSTK(I)) = VTYP
                        END IF
  150                CONTINUE
                  END IF

               ELSE

C               --Function has a defined type, parameters may be mixed types

                  VTYP = FNCTYP(INXENT(NENT))
                  DO 160 I = IPARM, ILAST
                     IF (TYPSTK(I) .EQ. '?') VTYP = TYPSTK(I)
                     IF ((TYPSTK(I) .EQ. 'H')
     &                  .OR. (TYPSTK(I) .EQ. 'G')) THEN
                        VSZENT(INXSTK(I)) = '*'
                     END IF
  160             CONTINUE
               END IF

               NEV = 0
               NEWEV = NUMEV + 1
               IEVOLD = IDEFEV
               DO 170 I = IPARM, ILAST
                  IF (TYPSTK(I) .EQ. 'E') THEN
                     NEV = NEV + 1
                     IEVNEW = IEVENT(INXSTK(I))
                     CALL MAKEVO (NELBLK, ISEVOK,
     &                  IEVOLD, IEVNEW, NEWEV, IEVOLD)
                  END IF
  170          CONTINUE
               IF (NEV .GT. 0) THEN
                  NEWEV = IEVOLD
                  IF (NEWEV .GT. NUMEV) NUMEV = NEWEV
                  IEVENT(NENT) = NEWEV
               END IF

               TYPSTK(ITOS) = VTYP
               INXSTK(ITOS) = NENT
               VSZENT(NENT) = VTYP
            END IF

         ELSE

C         --Do simple check that variable types match since postfix string
C         --is incorrect

            IF (TYPENT(NENT) .EQ. 'V') THEN

C            --Check that there are no mixed array types

               IF (TYPSTK(1) .NE. ' ') THEN
                  IF (TYPSTK(1) .NE. VTYP) THEN
                     IF (TYPSTK(1) .EQ. 'H') THEN
                        CONTINUE
                     ELSE IF (TYPSTK(1) .EQ. 'G') THEN
                        IF (VTYP .EQ. 'H') VTYP = TYPSTK(1)
                     ELSE IF ((TYPSTK(1) .EQ. 'N')
     &                  .OR. (TYPSTK(1) .EQ. 'E')) THEN
                        IF ((VTYP .EQ. 'N') .OR. (VTYP .EQ. 'E')) THEN
                           IF (TYPSTK(1) .NE. VTYP) VTYP = '?'
                        ELSE
                           IF (VTYP .NE. '?') VTYP = '?'
                        END IF
                     ELSE IF (TYPSTK(1) .EQ. '?') THEN
                        VTYP = TYPSTK(1)
                     END IF
                  END IF
               END IF

            ELSE IF (TYPENT(NENT) .EQ. 'F') THEN

C            --If defined-type function is present, do not check array types

               IF (FNCTYP(INXENT(NENT)) .NE. ' ') THEN
                  IF (TYPSTK(1) .NE. '?') TYPSTK(1) = ' '
               END IF
            END IF

         END IF
  180 CONTINUE

      IF (TYPENT(1) .EQ. 'V') THEN

         NAME = NAMENT(1)

C      --Specifier is not valid on the assigned variable

         IF (VALENT(1) .GT. 0) THEN
            NERR = NERR + 1
            CALL PRTERR ('CMDSPEC', 'Node/element specifier not allowed'
     &         // ' on assigned variable')
         END IF

         IF (ITMENT(1) .NE. ICURTM) THEN
            NERR = NERR + 1
            CALL PRTERR ('CMDSPEC', 'Time specifier not allowed'
     &         // ' on assigned variable')
         END IF

         INXENT(1) = -1
C         --Set INXENT -1 for reference to LHS variable
         VALENT(1) = 0
         IF (TYPEQV .EQ. 'E') THEN
            IEVENT(1) = NEWEV
         ELSE
            IEVENT(1) = -999
         END IF

      ELSE
         NAME = ' '
      END IF

C   --Assign and check the type for the assigned variable

      IF (TYPSTK(1) .EQ. '?') THEN
         TYPEQV = 'G'
         NERR = NERR + 1
         CALL PRTERR ('CMDSPEC',
     &      'Nodal and element variables are mixed')
      ELSE
         TYPEQV = TYPSTK(1)
      END IF
      VSZENT(1) = TYPEQV

      IF (NAME .EQ. 'TIME') THEN
         IF (TYPEQV .NE. 'G') THEN
            NERR = NERR + 1
            CALL PRTERR ('CMDSPEC', 'TIME must be a global variable')
         END IF
         TYPEQV = 'T'
      END IF

C   --If the assigned variable does not exist, add it to the list,
C   --else check that its type has not changed

      NUMLHS = MAXVAR - IXLHS + 1
      IVAR = LOCSTR (NAME, NUMLHS, NAMVAR(IXLHS))

      IF (IVAR .LE. 0) THEN
         IXLHS = IXLHS - 1
         IF (NUMINP .LT. IXLHS) THEN
            NAMVAR(IXLHS) = NAME
            TYPVAR(IXLHS) = TYPEQV
            IDVAR(IXLHS) = -999
            ISTVAR(ICURTM,IXLHS) = 1
            ISTVAR(ILSTTM,IXLHS) = 0
            ISTVAR(IONETM,IXLHS) = 0
            IF (TYPEQV .EQ. 'E') THEN
               IEVVAR(IXLHS) = NEWEV
            ELSE
               IEVVAR(IXLHS) = -999
            END IF
         END IF
      ELSE
         IVAR = IXLHS + IVAR - 1
         IF (TYPEQV .NE. TYPVAR(IVAR)) THEN
            IF (TYPEQV .NE. ' ') THEN
               NERR = NERR + 1
               CALL PRTERR ('CMDSPEC', NAME(:LENSTR(NAME))
     &            // ' was previously assigned as a different type')
            END IF
         END IF

C      --Change DELETEd or SAVEd variable to assigned variable
         IF ((NERR .EQ. 0) .AND. (ISTVAR(ICURTM,IVAR) .NE. 1))
     &      ISTVAR(ICURTM,IVAR) = 1
      END IF

      GOTO 200

  190 CONTINUE
      NERR = NERR + 1
      CALL PRTERR ('PROGRAM', 'Program stack problem in LOCEQV')

  200 CONTINUE

      IF (NUMINP .GE. IXLHS) THEN
         NERR = NERR + 1
         N = NUMINP + (MAXVAR - IXLHS + 1)
         CALL INTSTR (1, 0, N, STRA, LSTRA)
         CALL INTSTR (1, 0, MAXVAR, STRB, LSTRB)
         CALL PRTERR ('CMDSPEC',
     &      'Too many variable names to store, '
     &      // STRA(:LSTRA) // ' > ' // STRB(:LSTRB))
      END IF

      IF (NERR .GT. 0) THEN
         NUMINP = NINP
         IXLHS = ILHS
      END IF

      RETURN
      END
