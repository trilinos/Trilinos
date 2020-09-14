C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE USIDS (IFLD, INTYP, CFIELD, IFIELD,
     &   LOLD1, IDOLD1, LOLD2, IDOLD2, LENNEW, IDNEW, *)
C=======================================================================

C   --*** USIDS *** (GEN3D) Read list of IDs
C   --   Written by Amy Gilkey - revised 05/21/86
C   --
C   --USIDS processes a list of IDs.  The ID list is checked for
C   --repetitions with the existing list and with itself.  Repetitions
C   --are flagged with an error message and ignored.
C   --
C   --Parameters:
C   --   IFLD - IN/OUT - the free-field index
C   --   INTYP - IN - the free-field type
C   --   CFIELD - IN - the free-field characters
C   --   IFIELD - IN - the free-field integers
C   --   LOLD1, LOLD2 - IN - the length of the existing IDs
C   --   IDOLD1, IDOLD2 - IN - the existing IDs
C   --   LENNEW - OUT - the length of the returned IDs
C   --   IDNEW - OUT - the returned IDs
C   --   * - return statement iff serious error

      PARAMETER (MAXSET=10)

      INTEGER INTYP(*)
      CHARACTER*8 CFIELD(*)
      INTEGER IFIELD(*)
      INTEGER IDOLD1(*), IDOLD2(*)
      INTEGER IDNEW(*)

      CHARACTER*5 STRA
      LOGICAL DUPID

      CALL INIINT (MAXSET, 0, IDNEW)
      LENNEW = 0

   10 CONTINUE
      IF (INTYP(IFLD) .GE. -1) THEN

         CALL FFINTG (IFLD, INTYP, IFIELD,
     &      'set id', 0, ID, *50)

         IF (LENNEW .GE. MAXSET) THEN
            CALL INTSTR (1, 0, MAXSET, STRA, LSTRA)
            CALL PRTERR ('CMDERR',
     &         'Number of IDs must be less than ' // STRA(:LSTRA))
            GOTO 60
         END IF

         DUPID = .FALSE.
         DO 20 I = 1, LOLD1
            IF (ID .EQ. IDOLD1(I)) DUPID = .TRUE.
   20    CONTINUE
         DO 30 I = 1, LOLD2
            IF (ID .EQ. IDOLD2(I)) DUPID = .TRUE.
   30    CONTINUE
         DO 40 I = 1, LENNEW
            IF (ID .EQ. IDNEW(I)) DUPID = .TRUE.
   40    CONTINUE
         IF (DUPID) THEN
            CALL INTSTR (1, 0, ID, STRA, LSTRA)
            CALL PRTERR ('CMDWARN',
     &         'Duplicate ID ' // STRA(:LSTRA) // ' ignored')
            GOTO 50
         END IF

         LENNEW = LENNEW + 1
         IDNEW(LENNEW) = ID
   50    CONTINUE
         GOTO 10
      END IF

   60 CONTINUE
      RETURN
      END
