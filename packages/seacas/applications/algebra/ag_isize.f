C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

C=======================================================================
      INTEGER FUNCTION ISIZE (VTYP)
C=======================================================================

C   --*** ISIZE *** (ALGEBRA) Return the size of a stack variable
C   --   Written by Amy Gilkey - revised 06/22/87
C   --
C   --ISIZE returns the size of a stack variable type.
C   --
C   --Parameters:
C   --   VTYP - IN - the stack variable type (as in VSZENT of /ENT/)
C   --
C   --Common Variables:
C   --   Uses NUMNP, NUMEL of /DBNUMS/

      include 'ag_dbnums.blk'

      CHARACTER VTYP

      IF ((VTYP .EQ. 'T') .OR.
     &   (VTYP .EQ. 'H') .OR. (VTYP .EQ. 'G')) THEN
         ISIZE = 1
      ELSE IF (VTYP .EQ. 'N') THEN
         ISIZE = NUMNP
      ELSE IF (VTYP .EQ. 'E') THEN
         ISIZE = NUMEL
      ELSE
         ISIZE = MAX (NUMNP, NUMEL, 1)
      END IF

      RETURN
      END
