C Copyright(C) 2009-2017 National Technology & Engineering Solutions of
C Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
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
C     * Neither the name of NTESS nor the names of its
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

C $Log: pick2d.f,v $
C Revision 1.2  2009/03/25 12:36:46  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.1  1994/04/07 20:06:50  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:54:36  gdsjaar
c Added RCS Id and Log to all files
c
C=======================================================================
      SUBROUTINE PICK2D (PROMPT, CURSOR, MIDDEF, IFLD, INTYP, RFIELD,
     &   XRES, YRES, *)
C=======================================================================

C   --*** PICK2D *** (MESH) Pick a window point
C   --   Written by Amy Gilkey - revised 02/05/88
C   --
C   --PICK2D picks a user-selected window point.  For cursor input, PICK2D
C   --puts a cursor on the screen, waits until the user positions the cursor,
C   --then translates the position into a 2D point.  For character input,
C   --PICK2D gets an X,Y pair from the free-format fields (defaults are
C   --not allowed).
C   --
C   --Parameters:
C   --   PROMPT - IN - the point requested
C   --   CURSOR - IN - true iff cursor input
C   --   MIDDEF - IN - true if the middle of the displayed mesh is to be the
C   --      starting position, else use the last selected position
C   --      (cursor input only)
C   --   IFLD - IN/OUT - the free-format field index (character input only)
C   --   INTYP - IN - the free-format field types (character input only)
C   --   RFIELD - IN - the free-format real fields (character input only)
C   --   XRES, YRES - OUT - the returned point
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
      LOGICAL MIDDEF
      INTEGER INTYP(*)
      REAL RFIELD(*)

      CHARACTER*80 XPRMPT
      CHARACTER CH

      IF (CURSOR) THEN
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

      ELSE

C      --Get the X,Y pair from the free-format fields, no defaults

         XPRMPT = PROMPT // ' X, Y'
         CALL FFNEED (IFLD, INTYP, 'R', 2,
     &      XPRMPT(:LENSTR(XPRMPT)), *100)
         CALL FFREAL (IFLD, INTYP, RFIELD,
     &      'X', 0.0, XRES, *100)
         CALL FFREAL (IFLD, INTYP, RFIELD,
     &      'Y', 0.0, YRES, *100)
      END IF

      RETURN

  100 CONTINUE
      RETURN 1
      END
