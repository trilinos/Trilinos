C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE SCACAL (NAME, IVAR, USESEL, IELBST,
     &   ISTMN, ICALC)
C=======================================================================

C   --*** SCACAL *** (BLOT) Return variable scaling if already calculated
C   --   Written by Amy Gilkey - revised 11/03/87
C   --
C   --SCACAL determines if the variables has already been scaled.  If not,
C   --a scaling message is printed.  For element variables, the message
C   --is printed if all element blocks are not selected.
C   --
C   --Parameters:
C   --   NAME - IN - the variable name
C   --   IVAR - IN - the variable index
C   --   USESEL - IN - use the element blocks selected array iff true,
C   --      else all selected
C   --   IELBST - IN - the element block status (>0 if selected)
C   --   ISTMN - IN - the "state" of the minimums and maximums:
C   --      <0 if not calculated
C   --       0 if no elements in element block
C   --      +n if calculated
C   --   ICALC - OUT - the scaling flag
C   --      0 = all element blocks selected (if element) and min/max calculated
C   --      1 = not all element blocks selected (if element), but min/max
C   --          calculated for each selected element block
C   --      2 = min/max not calculated for some element block
C   --
C   --Common Variables:
C   --   Uses NELBLK, NVARNP, NVAREL of /DBNUMS/

      include 'dbnums.blk'

      CHARACTER*(*) NAME
      LOGICAL USESEL
      INTEGER IELBST(*)
      INTEGER ISTMN(0:*)

      CHARACTER TYP

C   --Get variable type

      CALL DBVTYP_BL (IVAR, TYP, IDUM)

      IF ((TYP .EQ. 'H') .OR. (TYP .EQ. 'G')
     &   .OR. (TYP .EQ. 'N')) THEN

C      --History, global or nodal variable, check min/max already computed

         IF (ISTMN(0) .GE. 0) THEN
            ICALC = 0
         ELSE
            ICALC = 2
         END IF

      ELSE IF (TYP .EQ. 'E') THEN

C      --Element variable, check if all element blocks selected and if
C      --all-blocks min/max already computed

         IF (USESEL) THEN
            CALL CNTELB (IELBST, NELBLK, NUMON, NUMSEL)
         ELSE
            NUMSEL = NELBLK
         END IF

         IF ((NUMSEL .GE. NELBLK) .AND. (ISTMN(0) .GE. 0)) THEN
            ICALC = 0

         ELSE

C         --If not all element blocks selected, determine if each min/max for
C         --variable already computed

            ICALC = 1
            DO 100 IELB = 1, NELBLK
               IF (ISTMN(IELB) .LT. 0) ICALC = 2
  100       CONTINUE
         END IF

      END IF

      IF (ICALC .GT. 0) THEN
         WRITE (*, 10000) 'Scaling variable ', NAME(:LENSTR(NAME))
      END IF

      RETURN
10000  FORMAT (1X, 5A)
      END
