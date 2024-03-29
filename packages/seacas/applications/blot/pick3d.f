C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE PICK3D (PROMPT, CURSOR,
     &   NUMNPF, XN, YN, ZN, HIDENP, MIDDEF, IFLD, INTYP, RFIELD,
     &   XRES, YRES, ZRES, *)
C=======================================================================

C   --*** PICK3D *** (MESH) Pick an object point
C   --   Written by Amy Gilkey - revised 02/05/88
C   --
C   --PICK3D picks a user-selected object point.  For cursor input, PICK3D
C   --puts a cursor on the screen, waits until the user positions the cursor,
C   --then translates the position into a 3D point by finding the nearest
C   --node and returning its unrotated, deformed coordinates.  For character
C   --input, PICK3D gets an X,Y,Z set from the free-format fields (defaults
C   --are not allowed).
C   --
C   --Parameters:
C   --   PROMPT - IN - the point requested
C   --   CURSOR - IN - true iff cursor input
C   --   NUMNPF - IN - the number of nodes (cursor input only)
C   --   XN, YN, ZN - IN - the displayed mesh nodal coordinates
C   --      (cursor input only)
C   --   HIDENP - IN - node i is defined iff HIDENP(i) is true
C   --      (cursor input only)
C   --   MIDDEF - IN - true if the middle of the displayed mesh is to be the
C   --      starting position, else use the last selected position
C   --      (cursor input only)
C   --   IFLD - IN/OUT - the free-format field index (character input only)
C   --   INTYP - IN - the free-format field types (character input only)
C   --   RFIELD - IN - the free-format real fields (character input only)
C   --   NPRES - OUT - the returned node
C   --   XRES, YRES, ZRES - OUT - the returned point
C   --   * - return statement if error picking point
C   --
C   --Common Variables:
C   --   Sets DXLAST, DYLAST of /PICK/
C   --   Uses /PICK/

      PARAMETER (KLFT=1, KRGT=2, KBOT=3, KTOP=4, KNEA=5, KFAR=6)

      include 'pick.blk'

      CHARACTER*(*) PROMPT
      LOGICAL CURSOR
      REAL XN(*), YN(*), ZN(*)
      LOGICAL HIDENP(*)
      LOGICAL MIDDEF
      INTEGER INTYP(*)
      REAL RFIELD(*)

      CHARACTER*80 XPRMPT
      CHARACTER CH

      IF (CURSOR) THEN
         IF (.NOT. INITP) THEN
            CALL PRTERR ('CMDERR', 'No mesh is displayed')
            GOTO 110
         END IF

C      --Put up the cursor at the default location and wait for the user
C      --to select a point

         IF (PROMPT .NE. ' ') THEN
            XPRMPT = 'Pick node nearest ' // PROMPT // '...'
         ELSE
            XPRMPT = ' '
         END IF

         IF (MIDDEF) THEN
            DXRES = DXMID
            DYRES = DYMID
         ELSE
            DXRES = DXLAST
            DYRES = DYLAST
         END IF
         CALL GRIKEY (XPRMPT, DXRES, DYRES, CH, *110)

C      --Save the selected position

         DXLAST = DXRES
         DYLAST = DYRES

C      --Translate the point from device coordinates into object coordinates

         XRES = PKMESH(KLFT) + (DXRES - DMESH(KLFT)) * DVRAT
         YRES = PKMESH(KBOT) + (DYRES - DMESH(KBOT)) * DVRAT

C      --Get the coordinates of the node nearest the point

         DISMIN = 10E+30
         DO 100 INP = 1, NUMNPF
            IF (HIDENP(INP)) GOTO 100
            X = XRES - XN(INP)
            Y = YRES - YN(INP)
            DIS = X*X + Y*Y
            IF (DIS .LT. DISMIN) THEN
               DISMIN = DIS
               NPRES = INP
            END IF
  100    CONTINUE

         CALL UNROT (1, 1, PKRMAT, PKRCEN,
     &      XN(NPRES), YN(NPRES), ZN(NPRES), XRES, YRES, ZRES)

      ELSE

C      --Get the X,Y,Z point from the free-format fields, no defaults

         XPRMPT = PROMPT // ' X, Y, Z'
         CALL FFNEED (IFLD, INTYP, 'R', 3,
     &      XPRMPT(:LENSTR(XPRMPT)), *110)
         CALL FFREAL (IFLD, INTYP, RFIELD,
     &      'X', 0.0, XRES, *110)
         CALL FFREAL (IFLD, INTYP, RFIELD,
     &      'Y', 0.0, YRES, *110)
         CALL FFREAL (IFLD, INTYP, RFIELD,
     &      'Z', 0.0, ZRES, *110)
      END IF

      RETURN

  110 CONTINUE
      RETURN 1
      END
