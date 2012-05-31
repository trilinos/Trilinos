C Copyright(C) 2009 Sandia Corporation. Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
C certain rights in this software.
C         
C Redistribution and use in source and binary forms, with or without
C modification, are permitted provided that the following conditions are
C met:
C 
C     * Redistributions of source code must retain the above copyright
C       notice, this list of conditions and the following disclaimer.
C 
C     * Redistributions in binary form must reproduce the above
C       copyright notice, this list of conditions and the following
C       disclaimer in the documentation and/or other materials provided
C       with the distribution.
C     * Neither the name of Sandia Corporation nor the names of its
C       contributors may be used to endorse or promote products derived
C       from this software without specific prior written permission.
C 
C THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
C "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
C LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
C A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
C OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
C SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
C LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
C DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
C THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
C (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
C OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

C $Log: picknp.f,v $
C Revision 1.3  2009/03/25 12:36:46  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.2  2000/05/31 22:22:02  gdsjaar
C Fixed handling of the what and what3 commands.  It was accessing the
C hidenp array when it shouldn't have and passed in numnp instead of
C numel in one place.
C
C Revision 1.1  1994/04/07 20:07:02  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.3  1993/08/19  19:10:23  gdsjaar
c Nodes are not hidden in 2d, test was wrong.
c
c Revision 1.2  1990/12/14  08:54:46  gdsjaar
c Added RCS Id and Log to all files
c
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

      COMMON /PICK/   INITP, PKDEF,
     &   PKMESH(KTOP), PKRMAT(3,3), PKRCEN(3),
     &   DMESH(KTOP), DVRAT, DXMID, DYMID, DXLAST, DYLAST
      LOGICAL INITP, PKDEF

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
