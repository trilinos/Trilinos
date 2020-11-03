C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE MAKEVO (NELBLK, ISEVOK, IEV1, IEV2, IEVOUT, IEVSEL)
C=======================================================================

C   --*** MAKEVO *** (ALGEBRA) Make up ISEVOK entry
C   --   Written by Amy Gilkey - revised 05/12/88
C   --
C   --MAKEVO makes up an ISEVOK entry based on two input entries.  The output
C   --entry is true iff the two input entries are true.  A flag indicates
C   --if the output entry is equal to one of the input entries, but the
C   --output entry is still set.
C   --
C   --Parameters:
C   --   NELBLK - IN - the number of element blocks
C   --   ISEVOK - IN/OUT - the variable truth table
C   --   IEV1, IEV2 - IN - the ISEVOK indices of the two input entries;
C   --      either index may be 0, but not both
C   --   IEVOUT - OUT - the ISEVOK index of the output entry
C   --   IEVSEL - OUT - the ISEVOK index of the input entry that is equal
C   --      to the output entry or the index of the output entry if
C   --      neither entry is equal

      LOGICAL ISEVOK(NELBLK,*)

      LOGICAL OK1, OK2

      OK1 = (IEV1 .GE. 1)
      OK2 = (IEV2 .GE. 1)

      IF ((IEV1 .GE. 1) .AND. (IEV2 .GE. 1)) THEN
         DO 100 IELB = 1, NELBLK
            IF (ISEVOK(IELB,IEV1) .EQV. ISEVOK(IELB,IEV2)) THEN
               ISEVOK(IELB,IEVOUT) = ISEVOK(IELB,IEV1)
            ELSE
               IF (ISEVOK(IELB,IEV1)) THEN
                  OK1 = .FALSE.
               ELSE
                  OK2 = .FALSE.
               END IF
               ISEVOK(IELB,IEVOUT) = .FALSE.
            END IF
  100    CONTINUE
      ELSE IF (IEV1 .GE. 1) THEN
         DO 110 IELB = 1, NELBLK
            ISEVOK(IELB,IEVOUT) = ISEVOK(IELB,IEV1)
  110    CONTINUE
      ELSE IF (IEV2 .GE. 1) THEN
         DO 120 IELB = 1, NELBLK
            ISEVOK(IELB,IEVOUT) = ISEVOK(IELB,IEV2)
  120    CONTINUE
      ELSE
         DO 130 IELB = 1, NELBLK
            ISEVOK(IELB,IEVOUT) = .TRUE.
  130    CONTINUE
      END IF

      IF (OK1) THEN
         IEVSEL = IEV1
      ELSE IF (OK2) THEN
         IEVSEL = IEV2
      ELSE
         IEVSEL = IEVOUT
      END IF

      RETURN
      END
