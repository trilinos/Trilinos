C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE POSTFX (NUMENT, NAMENT, TYPENT, INXENT, VALENT, ITMENT,
     &   NERR)
C=======================================================================

C   --*** POSTFX *** (ALGEBRA) Convert the equation to postfix form
C   --   Written by Amy Gilkey - revised 08/18/87
C   --
C   --POSTFX converts the equation to postfix form.  Parenthesis
C   --are eliminated.  Nothing is done if NERR is non-zero.
C   --
C   --Conversion is done using an operator stack.  Operators are
C   --moved to the operator stack and returned to the equation in
C   --postfix order.  Non-operators are simply moved up in the equation.
C   --Any operator which is moved to the stack has precedence over those
C   --already in the stack.  An operator of lower precedence will move
C   --the preceding stack operator to the equation.  A special case is
C   --made for a function.  It is pushed onto the operator stack, and
C   --popped back into the equation at the end of the parameters.
C   --
C   --Parameters:
C   --   NUMENT - IN/OUT - the number of entries (as in /ENT../)
C   --   NAMENT - IN/OUT - the equation entries (as in /ENT../)
C   --   TYPENT - IN/OUT - the type of each equation entry (as in /ENT../)
C   --   INXENT - IN/OUT - based on TYPENT (as in /ENT../)
C   --   VALENT - IN/OUT - based on TYPENT (as in /ENT../)
C   --   ITMENT - IN/OUT - based on TYPENT (as in /ENT../)
C   --   NERR - IN/OUT - the number of errors in the equation, may be set

      include 'ag_namlen.blk'
      include 'ag_numeqn.blk'
      PARAMETER (ICURTM = 1, ILSTTM = 2, IONETM = 3)

      CHARACTER*(*) NAMENT(*)
      CHARACTER TYPENT(*)
      INTEGER INXENT(*)
      REAL VALENT(*)
      INTEGER ITMENT(*)

      CHARACTER*(maxnam) STKNAM(MAXENT)
      CHARACTER   STKTYP(MAXENT)
      INTEGER     STKINX(MAXENT)
      REAL        STKVAL(MAXENT)
      INTEGER     STKITM(MAXENT)
C      --STK.. - the operand stack; holds the corresponding /ENT../ values
      INTEGER     STKEQU(MAXENT)
C      --STKEQU - the IOPTBL index for this operand

      INTEGER IEQUIV(0:10)
      SAVE IEQUIV
C      --IEQUIV - the IOPTBL index for the given operator

      INTEGER IOPTBL(-1:5,-1:8)
      SAVE IOPTBL
C      --IOPTBL - chart of stack actions to be performed
C      --IOPTBL(i,j) - i = top stack operator, j = current operator
C      --   = 1   Pop operator from stack into the equation
C      --   = 2   Push operator onto the stack
C      --   = 3   Pop matching '(' off the stack
C      --   = 4   Delete from equation
C      --   = 0   Error

C                   #  ~  +  -  *  /  ^  (  )  =  ,
      DATA IEQUIV / 0, 1, 2, 2, 3, 3, 4, 5, 6, 7, 8 /
C      --Index -1 is reserved for the function name

C                                  F  #  ~  +  *  ^  (  )  =  ,
      DATA (IOPTBL(-1,J),J=-1,8) / 0, 1, 0, 1, 1, 1, 2, 1, 0, 1 /
      DATA (IOPTBL(0,J), J=-1,8) / 2, 4, 2, 2, 2, 2, 2, 0, 0, 0 /
      DATA (IOPTBL(1,J), J=-1,8) / 2, 1, 0, 1, 1, 1, 2, 1, 0, 1 /
      DATA (IOPTBL(2,J), J=-1,8) / 2, 1, 2, 1, 2, 2, 2, 1, 0, 1 /
      DATA (IOPTBL(3,J), J=-1,8) / 2, 1, 2, 1, 1, 2, 2, 1, 0, 1 /
      DATA (IOPTBL(4,J), J=-1,8) / 2, 1, 2, 1, 1, 2, 2, 1, 0, 1 /
      DATA (IOPTBL(5,J), J=-1,8) / 2, 0, 2, 2, 2, 2, 2, 3, 0, 4 /

      IF (NERR .NE. 0) GOTO 120

C   --Initialize the operator stack by putting '#' on it

      STKNAM(1) = '#'
      STKTYP(1) = 'O'
      STKINX(1) = 0
      STKVAL(1) = -999
      STKITM(1) = ICURTM
      STKEQU(1) = 0
      IOPTOS = 1
C      --IOPTOS - the top of the operator stack

C   --Put '#' at end of the equation

      NENT = NUMENT+1
      NAMENT(NENT) = STKNAM(1)
      TYPENT(NENT) = STKTYP(1)
      INXENT(NENT) = STKINX(1)
      VALENT(NENT) = STKVAL(1)
      ITMENT(NENT) = STKITM(1)

      IEQTOS = 2
C      --IEQTOS - the top of new postfix equation stack

      DO 110 NENT = 3, NUMENT+1

         IF ((TYPENT(NENT) .NE. 'O')
     &      .AND. (TYPENT(NENT) .NE. 'F')) THEN

C         --Move entry up in the equation

            IEQTOS = IEQTOS + 1
            IF (IEQTOS .NE. NENT) THEN
               NAMENT(IEQTOS) = NAMENT(NENT)
               TYPENT(IEQTOS) = TYPENT(NENT)
               INXENT(IEQTOS) = INXENT(NENT)
               VALENT(IEQTOS) = VALENT(NENT)
               ITMENT(IEQTOS) = ITMENT(NENT)
            END IF

         ELSE

  100       CONTINUE
            IF (TYPENT(NENT) .EQ. 'O') THEN
               IEQU = IEQUIV(INXENT(NENT))
            ELSE
               IEQU = -1
            END IF
            J = STKEQU(IOPTOS)
            IACT = IOPTBL(J,IEQU)

            IF (IACT .EQ. 1) THEN

C            --Pop operator from stack onto the equation, and redo

               IEQTOS = IEQTOS + 1
               NAMENT(IEQTOS) = STKNAM(IOPTOS)
               TYPENT(IEQTOS) = STKTYP(IOPTOS)
               INXENT(IEQTOS) = STKINX(IOPTOS)
               VALENT(IEQTOS) = STKVAL(IOPTOS)
               ITMENT(IEQTOS) = STKITM(IOPTOS)
               IOPTOS = IOPTOS - 1

               GOTO 100

            ELSE IF (IACT .EQ. 2) THEN

C            --Push operator onto the stack

               IOPTOS = IOPTOS + 1
               STKNAM(IOPTOS) = NAMENT(NENT)
               STKTYP(IOPTOS) = TYPENT(NENT)
               STKINX(IOPTOS) = INXENT(NENT)
               STKVAL(IOPTOS) = VALENT(NENT)
               STKITM(IOPTOS) = ITMENT(NENT)
               STKEQU(IOPTOS) = IEQU

            ELSE IF (IACT .EQ. 3) THEN

C            --Pop matching '(' off the stack

               IOPTOS = IOPTOS - 1

            ELSE IF (IACT .EQ. 4) THEN

C            --Delete from equation
               CONTINUE

            ELSE

C            --Error

               NERR = NERR + 1
               CALL PRTERR ('PROGRAM', 'Postfix string problem')
               GOTO 120

            END IF
         END IF

  110 CONTINUE

      NUMENT = IEQTOS

  120 CONTINUE
      RETURN
      END
