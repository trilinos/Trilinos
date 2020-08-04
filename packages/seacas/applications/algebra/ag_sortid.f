C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE SORTID (STYP, ITOTAL, IBEGIN, IEND)
C=======================================================================

C   --*** SORTID *** (ALGEBRA) Gather equation variables of one type
C   --   Written by Amy Gilkey - revised 08/19/87
C   --
C   --SORTID gathers the variables of the specified type into one area
C   --of the /VAR../ arrays and orders these variables by IDVAR.
C   --The variables are sorted at the start of the specified range.
C   --
C   --Parameters:
C   --   STYP   - IN - the variable type to be sorted
C   --   ITOTAL - IN - the ending /VAR../ index of the variables to be ordered
C   --   IBEGIN - IN - the starting /VAR../ index of the variables of type STYP
C   --   IEND   - OUT - the ending /VAR../ index of the variables of type STYP
C   --
C   --Common Variables:
C   --   Sets NAMVAR, TYPVAR, IDVAR, ISTVAR, IEVVAR of /VAR../

      include 'ag_namlen.blk'
      include 'exodusII.inc'
      include 'ag_var.blk'

      CHARACTER STYP

      CHARACTER*(maxnam) NAMTMP
      CHARACTER TYPTMP
      INTEGER ISTTMP(3)

      IEND = IBEGIN - 1

      DO 110 NVAR = IBEGIN, ITOTAL

         NID = 0
         DO 100 I = NVAR, ITOTAL
            IF (TYPVAR(I) .EQ. STYP) THEN
               IF (NID .EQ. 0) NID = I
               IF (IDVAR(I) .LT. IDVAR(NID)) NID = I
            END IF
  100    CONTINUE
         IF (NID .EQ. 0) GOTO 120

         IEND = IEND + 1

         IF (NID .NE. NVAR) THEN
            NAMTMP = NAMVAR(NVAR)
            TYPTMP = TYPVAR(NVAR)
            IDTMP = IDVAR(NVAR)
            ISTTMP(1) = ISTVAR(1,NVAR)
            ISTTMP(2) = ISTVAR(2,NVAR)
            ISTTMP(3) = ISTVAR(3,NVAR)
            IEVTMP = IEVVAR(NVAR)

            NAMVAR(NVAR) = NAMVAR(NID)
            TYPVAR(NVAR) = TYPVAR(NID)
            IDVAR(NVAR) = IDVAR(NID)
            ISTVAR(1,NVAR) = ISTVAR(1,NID)
            ISTVAR(2,NVAR) = ISTVAR(2,NID)
            ISTVAR(3,NVAR) = ISTVAR(3,NID)
            IEVVAR(NVAR) = IEVVAR(NID)

            NAMVAR(NID) = NAMTMP
            TYPVAR(NID) = TYPTMP
            IDVAR(NID) = IDTMP
            ISTVAR(1,NID) = ISTTMP(1)
            ISTVAR(2,NID) = ISTTMP(2)
            ISTVAR(3,NID) = ISTTMP(3)
            IEVVAR(NID) = IEVTMP
         END IF
  110 CONTINUE

  120 CONTINUE
      RETURN
      END
