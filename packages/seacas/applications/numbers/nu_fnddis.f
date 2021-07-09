C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE FNDDIS (NAMECO, NAMENV,
     &   DEFOK, NDIM, NVARNP, IXDEF, IYDEF, IZDEF)
C=======================================================================

C   --*** FNDDIS *** (MESH) Find displacement variables
C   --   Written by Amy Gilkey - revised 01/25/88
C   --
C   --FNDDIS finds the displacement variables.  The first two/three nodal
C   --variables are displacement variables if and only if they begin with
C   --'D' and end with the last character of the corresponding coordinate
C   --name.
C   --
C   --Parameters:
C   --   NAMECO - IN - the coordinate names
C   --   NAMENV - IN - the nodal variable names
C   --   DEFOK - OUT - true iff the displacement variables were found
C   --   NDIM - IN   - number of dimensions
C   --   NVARNP - IN - number of nodal variables
C   --   IXDEF, IYDEF, IZDEF - OUT - the indices of the displacement variables
C   --

      LOGICAL IS3DIM

      include 'exodusII.inc'
      CHARACTER*(MXSTLN) NAMECO(NDIM)
      CHARACTER*(MXSTLN) NAMENV(NVARNP)
      LOGICAL DEFOK

C ... Size of this string must be large enough to hold 3 names plus other info
C     See print near bottom.
      CHARACTER*400 STRING

      IS3DIM = (NDIM .EQ. 3)

      IF (NVARNP .GE. NDIM) THEN

C      --Locate displacement variables

         DEFOK = .TRUE.
         NDEF0 = 0
         LN = MAX (LENSTR (NAMENV(1)), 2)
         DO 100 I = 1, NDIM
            LC = LENSTR (NAMECO(I))
            IF ((NAMENV(I)(1:1) .NE. 'D')
     &         .OR. (NAMENV(I)(1:LN-1) .NE. NAMENV(1)(1:LN-1))
     &         .OR. (NAMENV(I)(LN:LN) .NE. NAMECO(I)(LC:LC)))
     &         DEFOK = .FALSE.
  100    CONTINUE

         IF (DEFOK) THEN
           ixdef = ndef0 + 1
           iydef = ndef0 + 2
           izdef = ndef0 + 3
         END IF

      ELSE
C      --If not at least NDIM vars, no displacement variables
         DEFOK = .FALSE.
      END IF

      IF (.NOT. DEFOK) THEN
         WRITE (*, 10000)
10000     FORMAT (/,
     &      ' A valid set of displacement functions cannot be found.', /
     &      ' Your mesh plots will not be deformed.')

         IXDEF = 0
         IYDEF = 0
         IZDEF = 0

      ELSE

         WRITE (*, *)
         IF (.NOT. IS3DIM) THEN
            WRITE (STRING, '(10A)') 'Displacement variables are ',
     &         NAMENV(IXDEF), ' ', NAMENV(IYDEF)
         ELSE
            WRITE (STRING, '(10A)') 'Displacement variables are ',
     &         NAMENV(IXDEF), ' ', NAMENV(IYDEF), ' ', NAMENV(IZDEF)
         END IF
         CALL SQZSTR (STRING, LSTR)
         WRITE (*, 10010) STRING(:LSTR)
10010     FORMAT (1X, 5A)

      END IF

      RETURN
      END
