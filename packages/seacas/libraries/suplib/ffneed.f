C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE FFNEED (IFLD, INTYP, FTYPE, NFLD, EXPECT, *)
C=======================================================================
C$Id: ffneed.f,v 1.3 2009/03/25 12:46:02 gdsjaar Exp $
C$Log: ffneed.f,v $
CRevision 1.3  2009/03/25 12:46:02  gdsjaar
CAdd copyright and license notice to all files.
C
CRevision 1.2  1993/08/19 18:40:54  gdsjaar
CFixed incorrect itype for integers
C
c Revision 1.1.1.1  1990/08/14  16:14:31  gdsjaar
c Testing
c
c Revision 1.1  90/08/14  16:14:29  gdsjaar
c Initial revision
c
c Revision 1.1  90/08/09  13:39:25  gdsjaar
c Initial revision
c

C   --*** FFNEED *** (FFLIB) Check free-field fields for type
C   --   Written by Amy Gilkey - revised 10/21/86
C   --
C   --FFNEED checks that the next free-format fields exist and are of the
C   --appropriate type.
C   --
C   --Parameters:
C   --   IFLD - IN - the index of the current field number, NOT incremented
C   --   INTYP - IN - the input types from the free-field reader
C   --   FTYPE - IN - the expected field type:
C   --      C for character, R for real, I for integer, other for character
C   --   NFLD - IN - the number of expected fields
C   --   EXPECT - IN - the value to expect string, for error
C   --   * - return statement if the fields do not exist or are not of the
C   --      expected type; message is printed

      INTEGER IFLD
      INTEGER INTYP(*)
      CHARACTER*(*) FTYPE
      INTEGER NFLD
      CHARACTER*(*) EXPECT

      CHARACTER*80 ERRMSG

      IF (FTYPE(1:1) .EQ. 'R') THEN
         ITYPE = 1
      ELSE IF (FTYPE(1:1) .EQ. 'I') THEN
         ITYPE = 2
      ELSE
         ITYPE = 0
      END IF

      DO 100 ICHK = IFLD, IFLD + NFLD - 1
         IF (INTYP(ICHK) .NE. ITYPE) THEN
            IF ((ITYPE .NE. 1) .OR. (INTYP(ICHK) .NE. 2)) THEN
               ERRMSG = 'Expected ' // EXPECT
               CALL PRTERR ('CMDERR', ERRMSG(:LENSTR(ERRMSG)))
               GOTO 110
            END IF
         END IF
  100 CONTINUE

      RETURN

  110 CONTINUE
      RETURN 1
      END
