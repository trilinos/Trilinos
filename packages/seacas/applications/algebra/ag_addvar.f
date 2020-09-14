C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details
C=======================================================================
      SUBROUTINE ADDVAR (NAMEQV, NVAR, TYPEQV, IDEQV, NINP, ILHS)
C=======================================================================
C   --*** ADDVAR *** (ALGEBRA) Add set of variables to list
C   --   Written by Amy Gilkey - revised 12/14/87
C   --
C   --ADDVAR adds all non-repeated input (INP) variables from a set of
C   --variables to the /VAR../ arrays as both input variables and assigned
C   --(LHS) variables.  History variables are all assigned the input
C   --variable name ".HISTORY".  Global variables are all assigned the input
C   --variable name ".GLOBAL".
C   --
C   --Parameters:
C   --   NAMEQV - IN - the names of the variables in the set
C   --   NVAR - IN - the number of variables in the set
C   --   TYPEQV - IN - the type of all variables in the set
C   --   IDEQV - IN - the index of the first variable in the set;
C   --      following are consecutive
C   --   NINP - IN - the number of input variables in /VAR../ that need
C   --      to be checked
C   --   ILHS - IN - the starting index of the assigned variables in /VAR../
C   --      that need to be checked
C   --
C   --Common Variables:
C   --   Sets NUMINP, IXLHS, NAMVAR, TYPVAR, IDVAR, ISTVAR, IEVVAR of /VAR../

      PARAMETER (ICURTM = 1, ILSTTM = 2, IONETM = 3)
      include 'ag_namlen.blk'
      include 'ag_var.blk'

      CHARACTER*(*) NAMEQV(*)
      CHARACTER TYPEQV
      CHARACTER*(maxnam) NAME

      NLHS = MAXVAR - ILHS + 1

      DO 100 IVAR = 1, NVAR

         NAME = NAMEQV(IVAR)

C      --Skip if already an assigned variable

         IF (LOCSTR (NAME, NLHS, NAMVAR(ILHS)) .LE. 0) THEN

C         --If not already an input variable, insert it in the input variables

            IF ((TYPEQV .NE. 'G') .OR. (IVAR .EQ. 1)) THEN

C            --Only insert first global variable with name .GLOBAL
               IF (TYPEQV .EQ. 'G') NAME = '.GLOBAL'

               IF (LOCSTR (NAME, NINP, NAMVAR) .LE. 0) THEN
                  NUMINP = NUMINP + 1
                  IF (NUMINP .LT. IXLHS) THEN
                     NAMVAR(NUMINP) = NAME
                     TYPVAR(NUMINP) = TYPEQV
                     IDVAR(NUMINP) = IDEQV + IVAR - 1
                     ISTVAR(ICURTM,NUMINP) = 1
                     ISTVAR(ILSTTM,NUMINP) = 0
                     ISTVAR(IONETM,NUMINP) = 0
                     IF (TYPEQV .EQ. 'E') THEN
                        IEVVAR(NUMINP) = IDEQV + IVAR - 1
                     ELSE
                        IEVVAR(NUMINP) = -999
                     END IF
                  END IF
               END IF
            END IF

C         --Insert into the LHS array as saved variable

            IXLHS = IXLHS - 1
            IF (NUMINP .LT. IXLHS) THEN
               NAMVAR(IXLHS) = NAMEQV(IVAR)
C               --Do not use NAME here since the global name is wrong
               TYPVAR(IXLHS) = TYPEQV
               IDVAR(IXLHS) = -999
               ISTVAR(ICURTM,IXLHS) = -2
               ISTVAR(ILSTTM,IXLHS) = 0
               ISTVAR(IONETM,IXLHS) = 0
               IF (TYPEQV .EQ. 'E') THEN
                  IEVVAR(IXLHS) = IDEQV + IVAR - 1
               ELSE
                  IEVVAR(IXLHS) = -999
               END IF
            END IF
         END IF
  100 CONTINUE

      RETURN
      END
