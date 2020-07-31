C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE LNKVAR (*)
C=======================================================================

C   --*** LNKVAR *** (ALGEBRA) Link equation references with variables
C   --   Written by Amy Gilkey - revised 11/30/87
C   --
C   --LNKVAR links the variables in the equation with the /VAR../ arrays.
C   --Each equation variable (from /ENT../) is found in the /VAR../
C   --arrays, and the /VAR../ index is stored in the INXENT array.
C   --Note that since the /VAR../ arrays were constructed from the
C   --equation variables, each variable will be found.
C   --
C   --Parameters:
C   --   * - return statement if an error is found; message is printed
C   --
C   --Common Variables:
C   --   Sets INXENT, VALENT of /ENT../
C   --   Uses NAMENT, TYPENT of /ENT../
C   --   Uses NAMVAR, TYPVAR, IDVAR, ISTVAR, NUMINP, IXLHS of /VAR../

      include 'ag_namlen.blk'
      include 'exodusII.inc'
      include 'ag_numeqn.blk'
      include 'ag_ent.blk'
      include 'ag_var.blk'

      PARAMETER (ICURTM = 1, ILSTTM = 2, IONETM = 3)

      NUMLHS = MAXVAR - IXLHS + 1

      DO 110 NEQN = 1, NUMEQN
         DO 100 NENT = 1, NUMENT(NEQN)
            IF (TYPENT(NENT,NEQN) .EQ. 'V') THEN

C            --If searching for a LHS variable (INXENT < 0),
C            --   then search the assigned variables, NAMVAR(IXLHS..MAXVAR)
C            --   else search the input variables, NAMVAR(1..NUMINP)

               IF (INXENT(NENT,NEQN) .LT. 0) THEN
                  IV = LOCSTR (NAMENT(NENT,NEQN),
     &               NUMLHS, NAMVAR(IXLHS))
                  IF (IV .GT. 0) IV = IXLHS + IV - 1
               ELSE
                  IV = LOCSTR (NAMENT(NENT,NEQN), NUMINP, NAMVAR(1))
               END IF

               IF (IV .LE. 0) THEN
                  CALL PRTERR ('PROGRAM', 'Variable "'
     &               // NAMENT(NENT,NEQN) // '" cannot be found')
                  GOTO 120
               END IF

               INXENT(NENT,NEQN) = IV

               IF ((TYPVAR(IV) .EQ. 'T')
     &            .OR. (TYPVAR(IV) .EQ. 'G')) THEN
c               IF ((TYPVAR(IV) .EQ. 'T')
c     &            .OR. (TYPVAR(IV) .EQ. 'H')
c     &            .OR. (TYPVAR(IV) .EQ. 'G')) THEN
                  IF (NAMVAR(IV) .EQ. '.GLOBAL') THEN
c                  IF ((NAMVAR(IV) .EQ. '.HISTORY')
c     &               .OR. (NAMVAR(IV) .EQ. '.GLOBAL')) THEN
                     VALENT(NENT,NEQN) =
     &                  VALENT(NENT,NEQN) + IDVAR(IV) - 1
                  ELSE
                     VALENT(NENT,NEQN) = IDVAR(IV)
                  END IF
               END IF

            END IF
  100    CONTINUE
  110 CONTINUE
      RETURN

  120 CONTINUE
      RETURN 1
      END
