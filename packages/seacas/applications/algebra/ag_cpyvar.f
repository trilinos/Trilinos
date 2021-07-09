C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE CPYVAR (TYP, MAXNE, FROM, TO)
C=======================================================================

C   --*** CPYVAR *** (ALGEBRA) Copy a database variable entry
C   --   Written by Amy Gilkey - revised 03/15/88
C   --
C   --CPYVAR copies a variable entry into another entry.
C   --
C   --Parameters:
C   --   TYP - IN - the variable type, used to determine entry size
C   --   MAXNE - IN - the maximum entry size
C   --   FROM - IN - the variables to be copied
C   --   TO - OUT - the copied variables
C   --
C   --Common Variables:
C   --   Uses NVARNP, NVAREL of /DBNUMS/

      include 'ag_dbnums.blk'

      CHARACTER TYP
      REAL FROM(*), TO(*)

      IF ((TYP .EQ. 'T') .OR. (TYP .EQ. 'H') .OR. (TYP .EQ. 'G')) THEN
         N = MAXNE
      ELSE IF ((TYP .EQ. 'C') .OR. (TYP .EQ. 'N')) THEN
         N = NUMNP
      ELSE IF (TYP .EQ. 'E') THEN
         N = NUMEL
      END IF

      CALL CPYREA (N, FROM, TO)

      RETURN
      END
