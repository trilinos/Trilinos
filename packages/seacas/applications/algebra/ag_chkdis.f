C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details
C=======================================================================
      SUBROUTINE CHKDIS (NDIM, NAMECO, NVARNP, NAMENV, LN1, LN2)
C=======================================================================

C   --*** CHKDIS *** (ALGEBRA) Check displacement variables
C   --   Written by Amy Gilkey - revised 03/02/88
C   --
C   --CHKDIS finds the displacement variables.  The first two/three nodal
C   --variables are displacement variables if and only if they begin with
C   --'D' and end with the last character of the corresponding coordinate
C   --name.
C   --
C   --Parameters:
C   --   NDIM   - IN - the number of coordinates
C   --   NAMECO - IN - the coordinate names
C   --   NVARNP - IN - the number of nodal variables
C   --   NAMENV - IN - the nodal variable names

      include 'exodusII.inc'

      CHARACTER*(LN1) NAMECO(*)
      CHARACTER*(LN2) NAMENV(*)

      LOGICAL DEFOK

      IF (NVARNP .GE. NDIM) THEN
         DEFOK = .TRUE.
         LN = MAX (LENSTR (NAMENV(1)), 2)
         DO 100 I = 1, NDIM
            LC = LENSTR (NAMECO(I))
            IF ((NAMENV(I)(1:1) .NE. 'D')
     &         .OR. (NAMENV(I)(1:LN-1) .NE. NAMENV(1)(1:LN-1))
     &         .OR. (NAMENV(I)(LN:LN) .NE. NAMECO(I)(LC:LC)))
     &         DEFOK = .FALSE.
  100    CONTINUE

      ELSE
         DEFOK = .FALSE.
      END IF

      IF (.NOT. DEFOK) THEN
         CALL PRTERR ('WARNING', 'Output database will not have'
     &      // ' a valid set of displacement functions.')
      END IF

      RETURN
      END
