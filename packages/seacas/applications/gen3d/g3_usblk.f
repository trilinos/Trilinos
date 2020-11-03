C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE USBLK (IFLD, INTYP, CFIELD, IFIELD,
     &   NEWTYP, NELBLK, IDELB, BLKTYP, *)
C=======================================================================

C   --*** USBLK *** (GEN3D) Read list of element block IDs
C   --   Written by Amy Gilkey - revised 05/21/86
C   --
C   --USBLK processes a list of element block IDs.  If an element block
C   --is in the list, the block type is changed to the new block type.
C   --
C   --Parameters:
C   --   IFLD - IN/OUT - the free-field index
C   --   INTYP - IN - the free-field type
C   --   CFIELD - IN - the free-field characters
C   --   IFIELD - IN - the free-field integers
C   --   NEWTYP - IN - the element block type to be set
C   --   NELBLK - IN - the number of element blocks
C   --   IDELB - IN - the ids for each block
C   --   BLKTYP - IN/OUT - the element block type
C   --   * - return statement iff serious error

      PARAMETER (MAXSET=10)

      INTEGER INTYP(*)
      CHARACTER*8 CFIELD(*)
      INTEGER IFIELD(*)
      CHARACTER NEWTYP
      INTEGER IDELB(NELBLK)
      CHARACTER BLKTYP(NELBLK)

      LOGICAL FFEXST
      CHARACTER*5 ISTR

      IF (.NOT. FFEXST (IFLD, INTYP)) THEN
         CALL INISTR (NELBLK, NEWTYP, BLKTYP)
      END IF

   10 CONTINUE
      IF (FFEXST (IFLD, INTYP)) THEN
         CALL FFINTG (IFLD, INTYP, IFIELD,
     &      'block id', 0, ID, *20)
         IELB = LOCINT (ID, NELBLK, IDELB)
         IF (IELB .LE. 0) THEN
            CALL INTSTR (1, 0, ID, ISTR, LSTR)
            CALL PRTERR ('CMDERR',
     &         'Invalid block id ' // ISTR(:LSTR) // ', ignored')
            GOTO 20
         END IF
         BLKTYP(IELB) = NEWTYP
   20    CONTINUE
         GOTO 10
      END IF

      RETURN
      END
