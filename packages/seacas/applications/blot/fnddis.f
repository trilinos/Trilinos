C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE FNDDIS (NAMECO, NAMENV, DEFOK, IXDEF,
     &                   IYDEF, IZDEF, DEFFAC, NAMLEN)
C=======================================================================

C   --*** FNDDIS *** (MESH) Find displacement variables
C   --   Written by Amy Gilkey - revised 01/25/88
C   --
C   --FNDDIS finds the displacement variables.  The first two/three nodal
C   --variables are displacement variables if and only if they begin with
C   --'D' and end with the last character of the corresponding coordinate
C   --name.
C     NAMING CONVENTION FOR DISPLACEMENT VARIABLES
C     FIRST OPTION
C     1) The nodal variables must be the first two or three variables
C        in the nodal variable names list.
C     2) Nodal variables must begin with a 'D'
C     3) The last character of the nodal variable must be equal to
C        the last character in the coordinate name.
C     4) The names of the displacement variables must be the same
C        except for the last character
C     SECOND OPTION
C     1) The nodal variables must appear together in the nodal
C        variable names list.
C     2) The nodal variables must begin with "DISP"
C     3) The last character of the nodal variable must be equal to
C        the last character in the coordinate name.
C     4) The names of the displacement variables must be the same
C        except for the last character
C     Sample of valid displacement variable names
C     Coordinate Name: X, Y, Z
C     Nodal Variable Names: DISPLX, DISPLY, DISPLZ
C     Note: Nodal variable names should be written with capital letters
C   --
C   --Parameters:
C   --   NAMECO - IN  - the coordinate names
C   --   NAMENV - IN  - the nodal variable names
C   --   DEFOK  - OUT - true iff the displacement variables were found
C   --   IXDEF  - OUT - index of x displacement variable
C   --   IYDEF  - OUT - index of y displacement variable
C   --   IZDEF  - OUT - index of z displacement variable
C   --   DEFFAC - OUT - the default displacement magnification
C   --
C   --Common Variables:
C   --   Uses NDIM, NVARNP, NSTEPW of /DBNUMS/
C   --   Uses IS3DIM of /D3NUMS/

      include 'params.blk'
      include 'dbnums.blk'
      include 'd3nums.blk'

      CHARACTER*(MXSTLN) NAMECO(*)
      CHARACTER*(NAMLEN) NAMENV(*)
      LOGICAL DEFOK

      CHARACTER CDUM
      CHARACTER*1024 STRING

      IF ((NSTEPW .GT. 0) .AND. (NVARNP .GE. NDIM)) THEN

C      --Locate displacement variables

         DEFOK = .TRUE.
         NDEF0 = 0
         LN = MAX (LENSTR (NAMENV(1)), 2)
         DO 100 I = 1, NDIM
            LC = LENSTR (NAMECO(I))
C           Check if first character is a 'D'
C           Check if all displacement names are the same other
C           than last character
C           Check if the last character in the displacement nodal variable
C           name is equal to the last character in the coordinate name.
            IF ((NAMENV(I)(1:1) .NE. 'D')
     &         .OR. (NAMENV(I)(1:LN-1) .NE. NAMENV(1)(1:LN-1))
     &         .OR. (NAMENV(I)(LN:LN) .NE. NAMECO(I)(LC:LC)))
     &         DEFOK = .FALSE.
  100    CONTINUE

         IF (.NOT. DEFOK) THEN
           DO 300 J = 2, NVARNP
             IF (LENSTR(NAMENV(J)) .GE. 5) THEN
               IF ((NAMENV(J)(1:4) .EQ. "DISP")
     &           .AND. (J+NDIM-1 .LE. NVARNP)) THEN
                 DEFOK = .TRUE.
                 NDEF0 = J-1
                 LN = MAX (LENSTR (NAMENV(J)), 2)
                 DO 200 I = 1, NDIM
                   LC = LENSTR (NAMECO(I))
                   IF ((NAMENV(I+J-1)(1:4) .NE. "DISP")
     &               .OR. (NAMENV(I+J-1)(1:LN-1) .NE.
     &                     NAMENV(J)(1:LN-1))
     &               .OR. (NAMENV(I+J-1)(LN:LN) .NE.
     &                     NAMECO(I)(LC:LC))) THEN
                     DEFOK = .FALSE.
                   ENDIF
 200             CONTINUE
                 IF (DEFOK) THEN
                    EXIT
                 ENDIF
               ENDIF
             ENDIF
 300       CONTINUE
         ENDIF

         IF (DEFOK) THEN
C           find the index for the x displacement variable
            CALL DBVIX_BL ('N', NDEF0+1, IXDEF)
C           find the index for the y displacement variable
            CALL DBVIX_BL ('N', NDEF0+2, IYDEF)
            IF (IS3DIM) THEN
C              find the index for the z displacement variable
               CALL DBVIX_BL ('N', NDEF0+3, IZDEF)
            ELSE
               IZDEF = 0
            END IF
         END IF

      ELSE
C      --If no time steps, no displacement variables
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

         DEFFAC = 0.0

      ELSE

         WRITE (*, *)
         CALL DBVTYP_BL (IXDEF, CDUM, NXDEF)
         CALL DBVTYP_BL (IYDEF, CDUM, NYDEF)
         IF (.NOT. IS3DIM) THEN
            WRITE (STRING, '(10A)') 'Displacement variables are ',
     &         NAMENV(NXDEF), ' ', NAMENV(NYDEF)
         ELSE
            CALL DBVTYP_BL (IZDEF, CDUM, NZDEF)
            WRITE (STRING, '(10A)') 'Displacement variables are ',
     &         NAMENV(NXDEF), ' ', NAMENV(NYDEF), ' ', NAMENV(NZDEF)
         END IF
         CALL SQZSTR (STRING, LSTR)
         WRITE (*, 10010) STRING(:LSTR)
10010     FORMAT (1X, 5A)

         DEFFAC = -999.0
      END IF

      RETURN
      END
