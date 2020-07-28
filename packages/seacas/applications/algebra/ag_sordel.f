C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE SORDEL (ITOTAL, IBEGIN, IEND)
C=======================================================================

C   --*** SORDEL *** (ALGEBRA) Sort deleted variables to end of entries
C   --   Written by Amy Gilkey - revised 12/11/87
C   --
C   --SORDEL sorts the deleted variables so that they appear at the start
C   --of all non-deleted entries.
C   --
C   --Parameters:
C   --   ITOTAL - IN - the ending /VAR../ index of the variables to sort
C   --   IBEGIN - IN - the starting /VAR../ index of the deleted variables
C   --   IEND - OUT - the ending /VAR../ index of the deleted variables
C   --
C   --Common Variables:
C   --   Sets NAMVAR, TYPVAR, IDVAR, ISTVAR, IEVVAR of /VAR../

      PARAMETER (ICURTM = 1, ILSTTM = 2, IONETM = 3)
      include 'ag_namlen.blk'
      include 'ag_var.blk'

      CHARACTER*(maxnam) NAMTMP
      CHARACTER TYPTMP
      INTEGER ISTTMP(3)

      IEND = IBEGIN - 1

      DO 100 NVAR = IBEGIN, ITOTAL

         IF (ISTVAR(ICURTM,NVAR) .EQ. -1) THEN

            IEND = IEND + 1

            IF (IEND .NE. NVAR) THEN
               NAMTMP = NAMVAR(NVAR)
               TYPTMP = TYPVAR(NVAR)
               IDTMP = IDVAR(NVAR)
               ISTTMP(1) = ISTVAR(1,NVAR)
               ISTTMP(2) = ISTVAR(2,NVAR)
               ISTTMP(3) = ISTVAR(3,NVAR)
               IEVTMP = IEVVAR(NVAR)

               NAMVAR(NVAR) = NAMVAR(IEND)
               TYPVAR(NVAR) = TYPVAR(IEND)
               IDVAR(NVAR) = IDVAR(IEND)
               ISTVAR(1,NVAR) = ISTVAR(1,IEND)
               ISTVAR(2,NVAR) = ISTVAR(2,IEND)
               ISTVAR(3,NVAR) = ISTVAR(3,IEND)
               IEVVAR(NVAR) = IEVVAR(IEND)

               NAMVAR(IEND) = NAMTMP
               TYPVAR(IEND) = TYPTMP
               IDVAR(IEND) = IDTMP
               ISTVAR(1,IEND) = ISTTMP(1)
               ISTVAR(2,IEND) = ISTTMP(2)
               ISTVAR(3,IEND) = ISTTMP(3)
               IEVVAR(IEND) = IEVTMP
            END IF
         END IF
  100 CONTINUE

      RETURN
      END
