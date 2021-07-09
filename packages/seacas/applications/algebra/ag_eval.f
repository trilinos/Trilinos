C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details
C=======================================================================
      SUBROUTINE EVAL (STEP1, WSTEP1, MAXNE, MAXSTK,
     &   IXNODE, IXELB, IXELBO, IXELEM, ISEVOK, NEQN, NUMENT,
     &   NAMENT, TYPENT, INXENT, VALENT, ITMENT, IEVENT, VSZENT,
     &   VARVAL, STACK, *)
C=======================================================================

C   --*** EVAL *** (ALGEBRA) Evaluate the equations
C   --   Written by Amy Gilkey - revised 05/23/88
C   --
C   --EVAL evaluates the equations.  It manipulates the input variables
C   --and assigns the output variables.
C   --
C   --Parameters:
C   --   STEP1  - IN - true iff evaluating the first selected time step
C   --   WSTEP1 - IN - true iff evalutaing the first selected whole time step
C   --   MAXNE  - IN - the VARVAL and STACK dimension
C   --   MAXSTK - IN - the maximum STACK size
C   --   IXNODE - IN - the indices of the zoomed nodes (iff ISZOOM)
C   --   IXELB  - IN - the number of elements in each block
C   --   IXELBO - IN - the number of elements in each output block
C   --   IXELEM - IN - the indices of the zoomed elements (iff ISZOOM)
C   --   ISEVOK - IN - the element variable truth table for the input,
C   --                 temporary, and output variables
C   --   NEQN   - IN - the equation number
C   --   NUMENT - IN - the number of entries (as in /ENT../)
C   --   NAMENT - IN - the equation entries (as in /ENT../)
C   --   TYPENT - IN - the type of each equation entry (as in /ENT../)
C   --   INXENT - IN - based on TYPENT (as in /ENT../)
C   --   VALENT - IN - based on TYPENT (as in /ENT../)
C   --   ITMENT - IN - based on TYPENT (as in /ENT../)
C   --   IEVENT - IN - the ISEVOK index for this entry (as in /ENT../)
C   --   VSZENT - IN - the type of variable on the stack after processing
C   --                 this entry (as in /ENT../)
C   --   VARVAL - IN/OUT - the input variables and returned output variables
C   --   STACK  - SCRATCH - the evaluation stack
C   --   MERR   - OUT - error flag
C   --
C   --Common Variables:
C   --   Uses ISTVAR of /VAR../
C   --   Uses NUMNP, NUMEL of /DBNUMS/
C   --   Uses EQNLIN of /EQNLNS/
C   --   Uses NUMELO, NUMNPO of /DBOUT/

      include 'ag_namlen.blk'
      include 'ag_numeqn.blk'
      include 'ag_var.blk'
      include 'ag_dbnums.blk'
      include 'ag_dbout.blk'
      include 'ag_eqnlns.blk'

      PARAMETER (ICURTM = 1, ILSTTM = 2, IONETM = 3)

      LOGICAL STEP1, WSTEP1
      INTEGER IXNODE(*)
      INTEGER IXELB(0:NELBLK)
      INTEGER IXELBO(0:NELBLK)
      INTEGER IXELEM(*)
      LOGICAL ISEVOK(NELBLK,*)
      CHARACTER*(maxnam) NAMENT(*)
      CHARACTER TYPENT(*)
      INTEGER INXENT(*)
      REAL VALENT(*)
      INTEGER ITMENT(*)
      INTEGER IEVENT(*)
      CHARACTER VSZENT(*)
      REAL VARVAL(MAXNE,*)
      REAL STACK(MAXNE,*)

      CHARACTER IVTYPE
      CHARACTER*5 STRA
      INTEGER MERR

      EXTERNAL MYSUM, MYMAX, MYMIN, MYSIGN, MYAMOD
      EXTERNAL myint, mynint, myabs, mydim, myexp, mylog,
     *  mylog10, mysin, mycos, mytan, myasin, myacos,
     *  myatan, myatan2, mysinh, mycosh, mytanh, mysqrt

      MERR = 0
      IOP2 = 0

      IF (NUMENT .LT. 3) RETURN

C   --Evaluate the equation

C     ITOS - the current top of the stack
      ITOS = 0

      DO 600 NENT = 3, NUMENT

         IF (TYPENT(NENT) .EQ. 'C') THEN

C         --Constant - store in stack(1,itos)

            ITOS = ITOS + 1
            IF (ITOS .GT. MAXSTK) GOTO 620
            IVTYPE = 'G'
            LENARY = 1
            STACK(1,ITOS) = VALENT(NENT)

         ELSE IF (TYPENT(NENT) .EQ. 'V') THEN

C         --Variable - store all its values on the stack

            ITOS = ITOS + 1
            IF (ITOS .GT. MAXSTK) GOTO 620
            ITM = ITMENT(NENT)
            IVAR = INXENT(NENT)
            ISTO = ISTVAR(ITM,IVAR)

            IVTYPE = TYPVAR(IVAR)
            LENARY = ISIZE (IVTYPE)

C         --Store either the array or single value to the stack

            NE = INT(VALENT(NENT))
            IF (NE .LE. 0) THEN
               IF (IVTYPE .NE. 'E') THEN
                  CALL CPYREA (LENARY, VARVAL(1,ISTO), STACK(1,ITOS))
               ELSE
                  IEV = IEVENT(NENT)
                  DO 100 IELB = 1, NELBLK
                     IF (ISEVOK(IELB,IEV)) THEN
                        N = IXELB(IELB) - IXELB(IELB-1)
                        J = IXELB(IELB-1)+1
                        CALL CPYREA (N, VARVAL(J,ISTO), STACK(J,ITOS))
                     END IF
  100             CONTINUE
               END IF
            ELSE
               IF ((IVTYPE .EQ. 'N') .OR. (IVTYPE .EQ. 'E'))
     &            IVTYPE = 'G'
               STACK(1,ITOS) = VARVAL(NE,ISTO)
C              Change length of array to indicate that only
C              1 value is being stored
               LENARY = 1
            END IF

         ELSE IF (TYPENT(NENT) .EQ. 'O') THEN

C         --Operator - pop operands and push result

            IF (INXENT(NENT) .LE. 1) THEN
C            --Handle unary operator

               IOP1 = ITOS
C               --IOP1 - the operand index

            ELSE
C            --Handle binary operator

               IOP1 = ITOS - 1
               IOP2 = ITOS
C               --IOP1, IOP2 - the operand indices

               ITOS = ITOS - 2 + 1
               IF (ITOS .LE. 0) GOTO 620
            END IF

C         --Get the appropriate ISEVOK index
            IF (IVTYPE .EQ. 'E') THEN
               IEV = IEVENT(NENT)
            ELSE
               IEV = 1
            END IF

C         --Get the number of selected values, <0 if all
            IF ((IVTYPE .EQ. 'N') .AND. (NUMNP .NE. NUMNPO)) THEN
               NUMIX = NUMNPO
            ELSE IF ((IVTYPE .EQ. 'E') .AND. (NUMEL .NE. NUMELO))
     &         THEN
               NUMIX = NUMELO
            ELSE
               NUMIX = -999
            END IF

            L = INXENT(NENT)
            GOTO (110, 120, 130, 140, 150, 160), L

            CALL INTSTR (1, 0, L, STRA, LSTRA)
            CALL PRTERR ('PROGRAM', 'Invalid operator "'
     &         // NAMENT(NENT)(1:1) // '" on stack, index '
     &         // STRA(:LSTRA))
            GOTO 630

C         --Unary minus
  110       CONTINUE
            call chkfnc ('~', nament(nent))
            CALL DOOPER ('~',
     &         IVTYPE, LENARY, NELBLK, IXELBO, ISEVOK(1,IEV),
     &         NUMIX, IXNODE, IXELEM,
     &         STACK(1,IOP1), RDUM, STACK(1,ITOS))
            GOTO 170

C         --Addition
  120       CONTINUE
            call chkfnc ('+', nament(nent))
            CALL DOOPER ('+',
     &         IVTYPE, LENARY, NELBLK, IXELBO, ISEVOK(1,IEV),
     &         NUMIX, IXNODE, IXELEM,
     &         STACK(1,IOP1), STACK(1,IOP2), STACK(1,ITOS))
            GOTO 170

C         --Subtraction
  130       CONTINUE
            call chkfnc ('-', nament(nent))
            CALL DOOPER ('-',
     &         IVTYPE, LENARY, NELBLK, IXELBO, ISEVOK(1,IEV),
     &         NUMIX, IXNODE, IXELEM,
     &         STACK(1,IOP1), STACK(1,IOP2), STACK(1,ITOS))
            GOTO 170

C         --Multiplication
  140       CONTINUE
            call chkfnc ('*', nament(nent))
            CALL DOOPER ('*',
     &         IVTYPE, LENARY, NELBLK, IXELBO, ISEVOK(1,IEV),
     &         NUMIX, IXNODE, IXELEM,
     &         STACK(1,IOP1), STACK(1,IOP2), STACK(1,ITOS))
            GOTO 170

C         --Division
  150       CONTINUE
            call chkfnc ('/', nament(nent))
            CALL CHKPAR ('EQ0', 'Divide by 0.0',
     &         IVTYPE, LENARY, NELBLK, IXELBO, ISEVOK(1,IEV),
     &         NUMIX, IXNODE, IXELEM,
     &         STACK(1,IOP2), RDUM, *630)
            CALL DOOPER ('/',
     &         IVTYPE, LENARY, NELBLK, IXELBO, ISEVOK(1,IEV),
     &         NUMIX, IXNODE, IXELEM,
     &         STACK(1,IOP1), STACK(1,IOP2), STACK(1,ITOS))
            GOTO 170

C         --Exponentiation
  160       CONTINUE
            call chkfnc ('^', nament(nent))
            CALL CHKPAR ('2EQ0', '0.0 raised to the 0.0 power',
     &         IVTYPE, LENARY, NELBLK, IXELBO, ISEVOK(1,IEV),
     &         NUMIX, IXNODE, IXELEM,
     &         STACK(1,IOP1), STACK(1,IOP2), *630)
            CALL DOOPER ('^',
     &         IVTYPE, LENARY, NELBLK, IXELBO, ISEVOK(1,IEV),
     &         NUMIX, IXNODE, IXELEM,
     &         STACK(1,IOP1), STACK(1,IOP2), STACK(1,ITOS))
            GOTO 170

  170       CONTINUE

         ELSE IF (TYPENT(NENT) .EQ. 'F') THEN

C         --Function - pop parameters and push result

            NPARM = INT(VALENT(NENT))
            IPARM = ITOS - NPARM + 1
C            --IPARM - the index of the first function parameter
            ISTO = ITMENT(NENT)
C            --ISTO - the storage index for time function results

C         --The equality ITOS = IPARM is assumed in many functions
            ITOS = IPARM
            IF (ITOS .LE. 0) GOTO 620

C         --Get the appropriate ISEVOK index
            IF (IVTYPE .EQ. 'E') THEN
               IEV = IEVENT(NENT)
            ELSE
               IEV = 1
            END IF

C         --Get the number of selected values, <0 if all
            IF ((IVTYPE .EQ. 'N') .AND. (NUMNP .NE. NUMNPO)) THEN
               NUMIX = NUMNPO
            ELSE IF ((IVTYPE .EQ. 'E') .AND. (NUMEL .NE. NUMELO))
     &         THEN
               NUMIX = NUMELO
            ELSE
               NUMIX = -999
            END IF

            L = INXENT(NENT)
            GOTO (180, 190, 200, 210, 220, 230, 240, 260,
     &         280, 290, 300, 310,
     &         320, 330, 340, 350, 360, 370, 380, 390, 400, 410,
     &         420, 430, 440, 450, 460,
     &         470, 480, 490,
     &         500, 510, 520, 530, 540,
     &         550, 560, 570, 580), L

            CALL INTSTR (1, 0, L, STRA, LSTRA)
            CALL PRTERR ('PROGRAM', 'Function "' // NAMENT(NENT)
     &         // '" is undefined, index ' // STRA(:LSTRA))
            GOTO 630

C         --AINT
  180       CONTINUE
            call chkfnc ('AINT', nament(nent))

            CALL DOFNC1 (MYINT,
     &         IVTYPE, LENARY, NELBLK, IXELBO, ISEVOK(1,IEV),
     &         NUMIX, IXNODE, IXELEM,
     &         STACK(1,IPARM), STACK(1,ITOS))
            GOTO 590

C         --ANINT
  190       CONTINUE
            call chkfnc ('ANINT', nament(nent))

            CALL DOFNC1 (MYNINT,
     &         IVTYPE, LENARY, NELBLK, IXELBO, ISEVOK(1,IEV),
     &         NUMIX, IXNODE, IXELEM,
     &         STACK(1,IPARM), STACK(1,ITOS))
            GOTO 590

C         --ABS
  200       CONTINUE
            call chkfnc ('ABS', nament(nent))

            CALL DOFNC1 (MYABS,
     &         IVTYPE, LENARY, NELBLK, IXELBO, ISEVOK(1,IEV),
     &         NUMIX, IXNODE, IXELEM,
     &         STACK(1,IPARM), STACK(1,ITOS))
            GOTO 590

C         --MOD
  210       CONTINUE
            call chkfnc ('MOD', nament(nent))
            CALL CHKPAR ('EQ0', 'MOD divide by 0.0',
     &         IVTYPE, LENARY, NELBLK, IXELBO, ISEVOK(1,IEV),
     &         NUMIX, IXNODE, IXELEM,
     &         STACK(1,IPARM+1), RDUM, *630)
            CALL DOFNC2 (MYAMOD,
     &         IVTYPE, LENARY, NELBLK, IXELBO, ISEVOK(1,IEV),
     &         NUMIX, IXNODE, IXELEM,
     &         STACK(1,IPARM+0), STACK(1,IPARM+1), STACK(1,ITOS))
            GOTO 590

C         --SIGN
  220       CONTINUE
            call chkfnc ('SIGN', nament(nent))
            CALL DOFNC2 (MYSIGN,
     &         IVTYPE, LENARY, NELBLK, IXELBO, ISEVOK(1,IEV),
     &         NUMIX, IXNODE, IXELEM,
     &         STACK(1,IPARM+0), STACK(1,IPARM+1), STACK(1,ITOS))
            GOTO 590

C         --DIM
  230       CONTINUE
            call chkfnc ('DIM', nament(nent))

            CALL DOFNC2 (MYDIM,
     &         IVTYPE, LENARY, NELBLK, IXELBO, ISEVOK(1,IEV),
     &         NUMIX, IXNODE, IXELEM,
     &         STACK(1,IPARM+0), STACK(1,IPARM+1), STACK(1,ITOS))
            GOTO 590

C         --MAX
  240       CONTINUE
            call chkfnc ('MAX', nament(nent))
            DO 250 IX = 1, NPARM-1
               CALL DOFNC2 (MYMAX,
     &            IVTYPE, LENARY, NELBLK, IXELBO, ISEVOK(1,IEV),
     &            NUMIX, IXNODE, IXELEM,
     &            STACK(1,IPARM+0), STACK(1,IPARM+IX), STACK(1,ITOS))
  250       CONTINUE
            GOTO 590

C         --MIN
  260       CONTINUE
            call chkfnc ('MIN', nament(nent))
            DO 270 IX = 1, NPARM-1
               CALL DOFNC2 (MYMIN,
     &            IVTYPE, LENARY, NELBLK, IXELBO, ISEVOK(1,IEV),
     &            NUMIX, IXNODE, IXELEM,
     &            STACK(1,IPARM+0), STACK(1,IPARM+IX), STACK(1,ITOS))
  270       CONTINUE
            GOTO 590

C         --SQRT
  280       CONTINUE
            call chkfnc ('SQRT', nament(nent))
            CALL CHKPAR ('LT0', 'SQRT of a negative number',
     &         IVTYPE, LENARY, NELBLK, IXELBO, ISEVOK(1,IEV),
     &         NUMIX, IXNODE, IXELEM,
     &         STACK(1,IPARM), RDUM, *630)

            CALL DOFNC1 (MYSQRT,
     &         IVTYPE, LENARY, NELBLK, IXELBO, ISEVOK(1,IEV),
     &         NUMIX, IXNODE, IXELEM,
     &         STACK(1,IPARM), STACK(1,ITOS))
            GOTO 590

C         --EXP
  290       CONTINUE
            call chkfnc ('EXP', nament(nent))

            CALL DOFNC1 (MYEXP,
     &         IVTYPE, LENARY, NELBLK, IXELBO, ISEVOK(1,IEV),
     &         NUMIX, IXNODE, IXELEM,
     &         STACK(1,IPARM), STACK(1,ITOS))
            GOTO 590

C         --LOG
  300       CONTINUE
            call chkfnc ('LOG', nament(nent))
            CALL CHKPAR ('LE0', 'LOG of a non-positive number',
     &         IVTYPE, LENARY, NELBLK, IXELBO, ISEVOK(1,IEV),
     &         NUMIX, IXNODE, IXELEM,
     &         STACK(1,IPARM), RDUM, *630)

            CALL DOFNC1 (MYLOG,
     &         IVTYPE, LENARY, NELBLK, IXELBO, ISEVOK(1,IEV),
     &         NUMIX, IXNODE, IXELEM,
     &         STACK(1,IPARM), STACK(1,ITOS))
            GOTO 590

C         --LOG10
  310       CONTINUE
            call chkfnc ('LOG10', nament(nent))
            CALL CHKPAR ('LE0', 'LOG10 of a non-positive number',
     &         IVTYPE, LENARY, NELBLK, IXELBO, ISEVOK(1,IEV),
     &         NUMIX, IXNODE, IXELEM,
     &         STACK(1,IPARM), RDUM, *630)

            CALL DOFNC1 (MYLOG10,
     &         IVTYPE, LENARY, NELBLK, IXELBO, ISEVOK(1,IEV),
     &         NUMIX, IXNODE, IXELEM,
     &         STACK(1,IPARM), STACK(1,ITOS))
            GOTO 590

C         --SIN
  320       CONTINUE
            call chkfnc ('SIN', nament(nent))

            CALL DOFNC1 (MYSIN,
     &         IVTYPE, LENARY, NELBLK, IXELBO, ISEVOK(1,IEV),
     &         NUMIX, IXNODE, IXELEM,
     &         STACK(1,IPARM), STACK(1,ITOS))
            GOTO 590

C         --COS
  330       CONTINUE
            call chkfnc ('COS', nament(nent))

            CALL DOFNC1 (MYCOS,
     &         IVTYPE, LENARY, NELBLK, IXELBO, ISEVOK(1,IEV),
     &         NUMIX, IXNODE, IXELEM,
     &         STACK(1,IPARM), STACK(1,ITOS))
            GOTO 590

C         --TAN
  340       CONTINUE
            call chkfnc ('TAN', nament(nent))

            CALL DOFNC1 (MYTAN,
     &         IVTYPE, LENARY, NELBLK, IXELBO, ISEVOK(1,IEV),
     &         NUMIX, IXNODE, IXELEM,
     &         STACK(1,IPARM), STACK(1,ITOS))
            GOTO 590

C         --ASIN
  350       CONTINUE
            call chkfnc ('ASIN', nament(nent))
            CALL CHKPAR ('GTABS1',
     *        'ASIN of a number greater than abs(1.0)',
     &         IVTYPE, LENARY, NELBLK, IXELBO, ISEVOK(1,IEV),
     &         NUMIX, IXNODE, IXELEM,
     &         STACK(1,IPARM), RDUM, *630)

            CALL DOFNC1 (MYASIN,
     &         IVTYPE, LENARY, NELBLK, IXELBO, ISEVOK(1,IEV),
     &         NUMIX, IXNODE, IXELEM,
     &         STACK(1,IPARM), STACK(1,ITOS))
            GOTO 590

C         --ACOS
  360       CONTINUE
            call chkfnc ('ACOS', nament(nent))
            CALL CHKPAR ('GTABS1',
     *        'ACOS of a number greater than abs(1.0)',
     &         IVTYPE, LENARY, NELBLK, IXELBO, ISEVOK(1,IEV),
     &         NUMIX, IXNODE, IXELEM,
     &         STACK(1,IPARM), RDUM, *630)

            CALL DOFNC1 (MYACOS,
     &         IVTYPE, LENARY, NELBLK, IXELBO, ISEVOK(1,IEV),
     &         NUMIX, IXNODE, IXELEM,
     &         STACK(1,IPARM), STACK(1,ITOS))
            GOTO 590

C         --ATAN
  370       CONTINUE
            call chkfnc ('ATAN', nament(nent))

            CALL DOFNC1 (MYATAN,
     &         IVTYPE, LENARY, NELBLK, IXELBO, ISEVOK(1,IEV),
     &         NUMIX, IXNODE, IXELEM,
     &         STACK(1,IPARM), STACK(1,ITOS))
            GOTO 590

C         --ATAN2
  380       CONTINUE
            call chkfnc ('ATAN2', nament(nent))
            CALL CHKPAR ('2EQ0', 'ATAN2 parameters are both 0.0',
     &         IVTYPE, LENARY, NELBLK, IXELBO, ISEVOK(1,IEV),
     &         NUMIX, IXNODE, IXELEM,
     &         STACK(1,IPARM+0), STACK(1,IPARM+1), *630)

            CALL DOFNC2 (MYATAN2,
     &         IVTYPE, LENARY, NELBLK, IXELBO, ISEVOK(1,IEV),
     &         NUMIX, IXNODE, IXELEM,
     &         STACK(1,IPARM+0), STACK(1,IPARM+1), STACK(1,ITOS))
            GOTO 590

C         --SINH
  390       CONTINUE
            call chkfnc ('SINH', nament(nent))

            CALL DOFNC1 (MYSINH,
     &         IVTYPE, LENARY, NELBLK, IXELBO, ISEVOK(1,IEV),
     &         NUMIX, IXNODE, IXELEM,
     &         STACK(1,IPARM), STACK(1,ITOS))
            GOTO 590

C         --COSH
  400       CONTINUE
            call chkfnc ('COSH', nament(nent))

            CALL DOFNC1 (MYCOSH,
     &         IVTYPE, LENARY, NELBLK, IXELBO, ISEVOK(1,IEV),
     &         NUMIX, IXNODE, IXELEM,
     &         STACK(1,IPARM), STACK(1,ITOS))
            GOTO 590

C         --TANH
  410       CONTINUE
            call chkfnc ('TANH', nament(nent))

            CALL DOFNC1 (MYTANH,
     &         IVTYPE, LENARY, NELBLK, IXELBO, ISEVOK(1,IEV),
     &         NUMIX, IXNODE, IXELEM,
     &         STACK(1,IPARM), STACK(1,ITOS))
            GOTO 590

C         --TMAG
  420       CONTINUE
            call chkfnc ('TMAG', nament(nent))
            CALL TMAG (IVTYPE, LENARY, NELBLK, IXELBO, ISEVOK(1,IEV),
     &         STACK(1,IPARM), STACK(1,IPARM+1), STACK(1,IPARM+2),
     &         STACK(1,IPARM+3), STACK(1,IPARM+4), STACK(1,IPARM+5),
     &         NUMIX, IXNODE, IXELEM, STACK(1,ITOS))
            GOTO 590

C         --PMAX
  430       CONTINUE
            call chkfnc ('PMAX', nament(nent))
            CALL PXN (.FALSE.,
     &         IVTYPE, LENARY, NELBLK, IXELBO, ISEVOK(1,IEV),
     &         NUMIX, IXNODE, IXELEM,
     &         STACK(1,IPARM+0), STACK(1,IPARM+1), STACK(1,IPARM+2),
     &         STACK(1,IPARM+3), STACK(1,IPARM+4), STACK(1,IPARM+5),
     &         STACK(1,ITOS), *630)
            GOTO 590

C         --PMIN
  440       CONTINUE
            call chkfnc ('PMIN', nament(nent))
            CALL PXN (.TRUE.,
     &         IVTYPE, LENARY, NELBLK, IXELBO, ISEVOK(1,IEV),
     &         NUMIX, IXNODE, IXELEM,
     &         STACK(1,IPARM+0), STACK(1,IPARM+1), STACK(1,IPARM+2),
     &         STACK(1,IPARM+3), STACK(1,IPARM+4), STACK(1,IPARM+5),
     &         STACK(1,ITOS), *630)
            GOTO 590

C         --PMAX2
  450       CONTINUE
            call chkfnc ('PMAX2', nament(nent))
            CALL PXN2 (.FALSE.,
     &         IVTYPE, LENARY, NELBLK, IXELBO, ISEVOK(1,IEV),
     &         NUMIX, IXNODE, IXELEM,
     &         STACK(1,IPARM+0), STACK(1,IPARM+1), STACK(1,IPARM+2),
     &         STACK(1,ITOS))
            GOTO 590

C         --PMIN2
  460       CONTINUE
            call chkfnc ('PMIN2', nament(nent))
            CALL PXN2 (.TRUE.,
     &         IVTYPE, LENARY, NELBLK, IXELBO, ISEVOK(1,IEV),
     &         NUMIX, IXNODE, IXELEM,
     &         STACK(1,IPARM+0), STACK(1,IPARM+1), STACK(1,IPARM+2),
     &         STACK(1,ITOS))
            GOTO 590

C         --IFLZ
  470       CONTINUE
            call chkfnc ('IFLZ', nament(nent))
            CALL DOIF ('IFLZ',
     &         IVTYPE, LENARY, NELBLK, IXELBO, ISEVOK(1,IEV),
     &         NUMIX, IXNODE, IXELEM,
     &         STACK(1,IPARM+0), STACK(1,IPARM+1), STACK(1,IPARM+2),
     &         STACK(1,ITOS))
            GOTO 590

C         --IFEZ
  480       CONTINUE
            call chkfnc ('IFEZ', nament(nent))
            CALL DOIF ('IFEZ',
     &         IVTYPE, LENARY, NELBLK, IXELBO, ISEVOK(1,IEV),
     &         NUMIX, IXNODE, IXELEM,
     &         STACK(1,IPARM+0), STACK(1,IPARM+1), STACK(1,IPARM+2),
     &         STACK(1,ITOS))
            GOTO 590

C         --IFGZ
  490       CONTINUE
            call chkfnc ('IFGZ', nament(nent))
            CALL DOIF ('IFGZ',
     &         IVTYPE, LENARY, NELBLK, IXELBO, ISEVOK(1,IEV),
     &         NUMIX, IXNODE, IXELEM,
     &         STACK(1,IPARM+0), STACK(1,IPARM+1), STACK(1,IPARM+2),
     &         STACK(1,ITOS))
            GOTO 590

C         --SUM
  500       CONTINUE
            CALL DOFNCG (MYSUM,
     &         IVTYPE, LENARY, NELBLK, IXELBO, ISEVOK(1,IEV),
     &         NUMIX, IXNODE, IXELEM,
     &         STACK(1,IPARM), STACK(1,ITOS))
            IVTYPE = 'G'
            LENARY = 1
            GOTO 590

C         --SMAX
  510       CONTINUE
            call chkfnc ('SMAX', nament(nent))
            CALL DOFNCG (MYMAX,
     &         IVTYPE, LENARY, NELBLK, IXELBO, ISEVOK(1,IEV),
     &         NUMIX, IXNODE, IXELEM,
     &         STACK(1,IPARM), STACK(1,ITOS))
            IVTYPE = 'G'
            LENARY = 1
            GOTO 590

C         --SMIN
  520       CONTINUE
            call chkfnc ('SMIN', nament(nent))
            CALL DOFNCG (MYMIN,
     &         IVTYPE, LENARY, NELBLK, IXELBO, ISEVOK(1,IEV),
     &         NUMIX, IXNODE, IXELEM,
     &         STACK(1,IPARM), STACK(1,ITOS))
            IVTYPE = 'G'
            LENARY = 1
            GOTO 590

C         --ENVMAX
  530       CONTINUE
            call chkfnc ('ENVMAX', nament(nent))
            ISTO = ITMENT(NENT)
            IF ((.NOT. STEP1) .OR. (.NOT. WSTEP1)) THEN
               CALL DOFNC2 (MYMAX,
     &            IVTYPE, LENARY, NELBLK, IXELBO, ISEVOK(1,IEV),
     &            NUMIX, IXNODE, IXELEM,
     &            STACK(1,IPARM), VARVAL(1,ISTO), STACK(1,ITOS))
            END IF
            CALL CPYREA (LENARY, STACK(1,ITOS), VARVAL(1,ISTO))
            GOTO 590

C         --ENVMIN
  540       CONTINUE
            call chkfnc ('ENVMIN', nament(nent))
            ISTO = ITMENT(NENT)
            IF ((.NOT. STEP1) .OR. (.NOT. WSTEP1)) THEN
               CALL DOFNC2 (MYMIN,
     &            IVTYPE, LENARY, NELBLK, IXELBO, ISEVOK(1,IEV),
     &            NUMIX, IXNODE, IXELEM,
     &            STACK(1,IPARM), VARVAL(1,ISTO), STACK(1,ITOS))
            END IF
            CALL CPYREA (LENARY, STACK(1,ITOS), VARVAL(1,ISTO))
            GOTO 590

C         --UHIST
  550       CONTINUE
            call chkfnc ('UHIST', nament(nent))
            CALL UHIST (NPARM, MAXNE, STACK(1,IPARM),
     &         NUMNP, NUMEL, *630)
            IVTYPE = 'H'
            LENARY = ISIZE (IVTYPE)
            GOTO 590

C         --UGLOB
  560       CONTINUE
            call chkfnc ('UGLOB', nament(nent))
            CALL UGLOB (NPARM, MAXNE, STACK(1,IPARM),
     &         NUMNP, NUMEL, *630)
            IVTYPE = 'G'
            LENARY = ISIZE (IVTYPE)
            GOTO 590

C         --UNODE
  570       CONTINUE
            call chkfnc ('UNODE', nament(nent))
            CALL UNODE (NPARM, MAXNE, STACK(1,IPARM),
     &         NUMNP, NUMEL, *630)
            IVTYPE = 'N'
            LENARY = ISIZE (IVTYPE)
            GOTO 590

C         --UELEM
  580       CONTINUE
            call chkfnc ('UELEM', nament(nent))
            CALL UELEM (NPARM, MAXNE, STACK(1,IPARM),
     &         NUMNP, NUMEL, *630)
            IVTYPE = 'E'
            LENARY = ISIZE (IVTYPE)
            GOTO 590

  590       CONTINUE
         END IF

         IF (IVTYPE .NE. VSZENT(NENT)) THEN
            IVTYPE = VSZENT(NENT)
            IOLDSZ = LENARY
            LENARY = ISIZE (IVTYPE)
            IF ((IOLDSZ .EQ. 1) .AND. (LENARY .GT. 1)) THEN
               CALL INIREA (LENARY, STACK(1,ITOS), STACK(1,ITOS))
            END IF
         END IF
  600 CONTINUE

C   --Store the equation result into VARVAL

      IANS = INXENT(1)
      ISTO = ISTVAR(ICURTM,IANS)
      NE = INT(VALENT(1))
      IF (NE .LE. 0) THEN
         CALL CPYREA (LENARY, STACK(1,ITOS), VARVAL(1,ISTO))
      ELSE
         VARVAL(NE,ISTO) = STACK(1,ITOS)
      END IF

      RETURN

  620 CONTINUE
      CALL INTSTR (1, 0, NENT, STRA, LSTRA)
      CALL PRTERR ('PROGRAM', 'Program stack problem, type '
     &   // TYPENT(NENT) // ', entry ' // STRA(:LSTRA))

  630 CONTINUE
      WRITE (*, 10000) EQNLIN(NEQN)(:LENSTR(EQNLIN(NEQN)))
10000  FORMAT (4X, A)

      RETURN 1
      END

      real function myint(parm)
      myint = int(parm)
      return
      end

      real function mynint(parm)
      mynint = nint(parm)
      return
      end

      real function myabs(parm)
      myabs = abs(parm)
      return
      end

      real function mydim(parm1, parm2)
      mydim = dim(parm1, parm2)
      return
      end

      real function myexp(parm)
      myexp = exp(parm)
      return
      end

      real function mylog(parm)
      mylog = log(parm)
      return
      end

      real function mylog10(parm)
      mylog10 = log10(parm)
      return
      end

      real function mysin(parm)
      mysin = sin(parm)
      return
      end

      real function mycos(parm)
      mycos = cos(parm)
      return
      end

      real function mytan(parm)
      mytan = tan(parm)
      return
      end

      real function myasin(parm)
      myasin = asin(parm)
      return
      end

      real function myacos(parm)
      myacos = acos(parm)
      return
      end

      real function myatan(parm)
      myatan = atan(parm)
      return
      end

      real function myatan2(parm1, parm2)
      myatan2 = atan2(parm1, parm2)
      return
      end

      real function mysinh(parm)
      mysinh = sinh(parm)
      return
      end

      real function mycosh(parm)
      mycosh = cosh(parm)
      return
      end

      real function mytanh(parm)
      mytanh = tanh(parm)
      return
      end

      real function mysqrt(parm)
      mysqrt = sqrt(parm)
      return
      end

