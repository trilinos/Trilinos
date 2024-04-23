C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE PICKN2 (PROMPT, NUMNPF, XN, YN, MIDDEF, NDRES, *)
C=======================================================================

C   --*** PICKN2 *** (MESH) Pick a node from window
C   --  written by Ray J. Meyers    21 June, 1990
C   --
C   --PICKN2 puts a cursor on the screen, finds the associated world
C   --       coordinates and then locates and returns the node closest
C   --       to those coordinates.
C   --
C   --Parameters:
C   --   PROMPT - IN - the point requested
C   --   NUMNPF - IN - the number of nodes
C   --   XN - IN - the x coordinates of the nodes
C   --   yN - IN - the y coordinates of the nodes
C   --   MIDDEF - IN - true if the middle of the displayed mesh is to be the
C   --      starting position, else use the last selected position
C   --   NDRES - OUT - the node id selected
C   --   * - return statement if error picking point
C   --
C   --Common Variables:
C   --   Sets DXLAST, DYLAST of /PICK/
C   --   Uses /PICK/

      DIMENSION XN(*), YN(*)
      PARAMETER (KLFT=1, KRGT=2, KBOT=3, KTOP=4, KNEA=5, KFAR=6)

      include 'pick.blk'

      CHARACTER*(*) PROMPT
      LOGICAL MIDDEF

      CHARACTER*80 XPRMPT
      CHARACTER CH

      IF (.NOT. INITP) THEN
         CALL PRTERR ('CMDERR', 'No mesh is displayed')
         GOTO 100
      END IF

C      --Put up the cursor at the default location and wait for the user
C      --to select a point

      IF (PROMPT .NE. ' ') THEN
         XPRMPT = 'Pick ' // PROMPT // '...'
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
      CALL GRIKEY (XPRMPT, DXRES, DYRES, CH, *100)

C      --Save the selected position

      DXLAST = DXRES
      DYLAST = DYRES

C      --Translate the point from device coordinates into mesh coordinates

      XRES = PKMESH(KLFT) + (DXRES - DMESH(KLFT)) * DVRAT
      YRES = PKMESH(KBOT) + (DYRES - DMESH(KBOT)) * DVRAT

C  -- FIND THE NEAREST NODE

      DISMIN = 10E+30
      DO 10 I = 1, NUMNPF
         X = XRES - XN(I)
         Y = YRES - YN(I)
         DIS = X*X + Y*Y
         IF(DIS .LT. DISMIN) THEN
             DISMIN = DIS
             NDRES = I
         END IF
10    CONTINUE

      RETURN

  100 CONTINUE
      RETURN 1
      END
