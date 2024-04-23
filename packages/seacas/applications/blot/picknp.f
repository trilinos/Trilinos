C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE PICKNP (PROMPT, CURSOR, WINDOW,
     &   NDIM, NUMNPF, XN, YN, ZN, HIDENP, MIDDEF,
     &   IFLD, INTYP, RFIELD, NPRES, *)
C=======================================================================

C   --*** PICKNP *** (MESH) Pick the nearest node
C   --   Written by Amy Gilkey - revised 03/07/88
C   --
C   --PICKNP picks a user-selected node.  For cursor input, PICKNP puts
C   --a cursor on the screen, waits until the user positions the cursor,
C   --then translates the position into a 3D point by finding the nearest
C   --node.  For character input, PICKNP gets an X,Y,Z set from the
C   --free-format fields (defaults are not allowed).
C   --
C   --Parameters:
C   --   PROMPT - IN - the point requested
C   --   CURSOR - IN - true iff cursor input
C   --   WINDOW - IN - true iff the coordinates are window versus object
C   --      coordinates
C   --   NDIM - IN - the number of dimensions
C   --   NUMNPF - IN - the number of nodes
C   --   XN, YN, ZN - IN - the displayed mesh nodal coordinates
C   --   HIDENP - IN - node i is defined iff HIDENP(i) is true
C   --      (only if WINDOW)
C   --   MIDDEF - IN - true if the middle of the displayed mesh is to be the
C   --      starting position, else use the last selected position
C   --      (cursor input only)
C   --   IFLD - IN/OUT - the free-format field index (character input only)
C   --   INTYP - IN - the free-format field types (character input only)
C   --   RFIELD - IN - the free-format real fields (character input only)
C   --   NPRES - OUT - the returned node
C   --   * - return statement if error picking point
C   --
C   --Common Variables:
C   --   Sets DXLAST, DYLAST of /PICK/
C   --   Uses /PICK/

      PARAMETER (KLFT=1, KRGT=2, KBOT=3, KTOP=4, KNEA=5, KFAR=6)

      include 'pick.blk'

      CHARACTER*(*) PROMPT
      LOGICAL CURSOR
      LOGICAL WINDOW
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

      ELSE

C      --Get the X,Y,Z point from the free-format fields, no defaults

         IF (WINDOW .OR. (NDIM .NE. 3)) THEN
            XPRMPT = PROMPT // ' X, Y'
         ELSE
            XPRMPT = PROMPT // ' X, Y, Z'
         END IF
         NEED = NDIM
         IF (WINDOW) NEED = 2
         CALL FFNEED (IFLD, INTYP, 'R', NEED,
     &      XPRMPT(:LENSTR(XPRMPT)), *110)
         CALL FFREAL (IFLD, INTYP, RFIELD,
     &      'X', 0.0, XRES, *110)
         CALL FFREAL (IFLD, INTYP, RFIELD,
     &      'Y', 0.0, YRES, *110)
         IF (.NOT. (WINDOW .OR. (NDIM .NE. 3))) THEN
            CALL FFREAL (IFLD, INTYP, RFIELD,
     &         'Z', 0.0, ZRES, *110)
         END IF
      END IF

C   --Get the coordinates of the node nearest the point

      DISMIN = 10E+30
      DO 100 INP = 1, NUMNPF
         IF (WINDOW .AND. NDIM .EQ. 3) THEN
            IF (HIDENP(INP)) GOTO 100
         END IF
         X = XRES - XN(INP)
         Y = YRES - YN(INP)
         IF (WINDOW .OR. (NDIM .NE. 3)) THEN
            DIS = X*X + Y*Y
         ELSE
            Z = ZRES - ZN(INP)
            DIS = X*X + Y*Y + Z*Z
         END IF
         IF (DIS .LT. DISMIN) THEN
            DISMIN = DIS
            NPRES = INP
         END IF
  100 CONTINUE

      RETURN

  110 CONTINUE
      RETURN 1
      END
