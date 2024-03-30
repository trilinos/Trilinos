C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE PICKN3 (PROMPT, NUMNPF, XN, YN, ZN, HIDENP,
     &                   MIDDEF, NODRES , *)
C=======================================================================

C   --*** PICKND *** (MESH) Pick a node and return node id
C   --   Written by Ray J. Meyers - 20 June, 1990
C   --
C   --PICK3D picks a user-selected node.  PICKND
C   --puts a cursor on the screen, waits until the user positions the cursor,
C   --then translates finds the nearest
C   --node and returns its node id.
C   --
C   --Parameters:
C   --   PROMPT - IN - the point requested
C   --   NUMNPF - IN - the number of nodes desired
C   --   XN, YN, ZN - IN - the displayed mesh nodal coordinates
C   --   HIDENP - IN - node i is defined iff HIDENP(i) is true
C   --   MIDDEF - IN - true if the middle of the displayed mesh is to be the
C   --      starting position, else use the last selected position
C   --   NODRES - the returned node number
C   --   * - return statement if error picking point
C   --
C   --Common Variables:
C   --   Sets DXLAST, DYLAST of /PICK/
C   --   Uses /PICK/

      PARAMETER (KLFT=1, KRGT=2, KBOT=3, KTOP=4, KNEA=5, KFAR=6)

      include 'pick.blk'

      CHARACTER*(*) PROMPT
      REAL XN(*), YN(*), ZN(*)
      LOGICAL HIDENP(*)
      LOGICAL MIDDEF

      CHARACTER*80 XPRMPT
      CHARACTER CH

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
            NODRES = INP
         END IF
  100  CONTINUE

      RETURN

  110 CONTINUE
      RETURN 1
      END
