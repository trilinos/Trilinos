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

C $Log: outlin.f,v $
C Revision 1.3  2009/03/25 12:36:46  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.2  2009/01/22 21:34:22  gdsjaar
C There were several inline dbnums common blocks. Replaced with the
C include so they all have the same size with the added variable types.
C
C Added minor support for nodeset and sideset variables.
C
C It can print the count and the names, but that is all
C at this time.
C
C Revision 1.1  1994/04/07 20:06:34  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.3  1993/09/24  17:32:40  gdsjaar
c Added an outline off/on command to toggle drawing of the view window outline
c
c Revision 1.2  1990/12/14  08:54:27  gdsjaar
c Added RCS Id and Log to all files
c
C=======================================================================
      SUBROUTINE OUTLIN (BLKCOL, *)
C=======================================================================

C   --*** OUTLIN *** (MESH) Outline viewport windows
C   --   Written by John Glick - 11/11/88
C   --
C   --OUTLIN outlines the viewport windows (symmetric views as one
C   --window, copies as separate windows)
C   --
C   --Parameters:
C   --   BLKCOL - IN/OUT - the user selected colors of the element blocks.
C   --                    BLKCOL(0) = 1 if the user defined material
C   --                                colors should be used in mesh plots.
C   --                              = -1 if program selected colors should
C   --                                be used.
C   --                    BLKCOL(i) = the user selected color of element
C   --                               block i:
C   --                                  -2 - no color selected by user.
C   --                                  -1 - black
C   --                                   0 - white
C   --                                   1 - red
C   --                                   2 - green
C   --                                   3 - yellow
C   --                                   4 - blue
C   --                                   5 - cyan
C   --                                   6 - magenta
C   --   * - return statement if the cancel function is active
C   --
C   --Common Variables:
C   --   Uses MSHDEF of /MSHOPC/
C   --   Uses XISSYM, YISSYM of /VIEWS/
C   --   Uses DVIEW of /LAYOUT/

      PARAMETER (KLFT=1, KRGT=2, KBOT=3, KTOP=4)

      include 'dbnums.blk'
      COMMON /MSHOPC/ MSHDEF(4), MSHNUM(4)
      CHARACTER*8 MSHDEF, MSHNUM
      COMMON /VIEWS/  MULTIM,
     &   XISSYM, YISSYM, XAXSYM, YAXSYM, LFTSYM, BOTSYM
      LOGICAL MULTIM, XISSYM, YISSYM, LFTSYM, BOTSYM
      COMMON /LAYOUD/ DVIEW(KTOP,4), WVIEW(KTOP,4)

      COMMON /LEGOPT/ DOQA(2), DOLEG(2), DOAXIS(2), DOBOX
      LOGICAL DOQA, DOLEG, DOAXIS, DOBOX

      INTEGER BLKCOL(0:NELBLK)

      LOGICAL GRABRT
      INTEGER NDEFVW, IXVW

      LOGICAL XISCOP, YISCOP
      REAL DAVIEW(KTOP)

      IF (.NOT. DOBOX) RETURN

      XISCOP = (.NOT. XISSYM) .AND. (MSHDEF(1) .NE. 'NONE')
      YISCOP = (.NOT. YISSYM) .AND. (MSHDEF(4) .NE. 'NONE')

      CALL UGRCOL (0, BLKCOL)

      IF (XISSYM) THEN
         DAVIEW(KLFT) = DVIEW(KLFT,2) - (DVIEW(KRGT,2) - DVIEW(KLFT,2))
      ELSE
         DAVIEW(KLFT) = DVIEW(KLFT,2)
      END IF
      DAVIEW(KRGT) = DVIEW(KRGT,2)
      IF (YISSYM) THEN
         DAVIEW(KBOT) = DVIEW(KBOT,2) - (DVIEW(KTOP,2) - DVIEW(KBOT,2))
      ELSE
         DAVIEW(KBOT) = DVIEW(KBOT,2)
      END IF
      DAVIEW(KTOP) = DVIEW(KTOP,2)

      IF (XISSYM .AND. YISSYM) THEN

         IF (GRABRT ()) RETURN 1
         CALL GRBOX ('L',
     &      DAVIEW(KLFT), DAVIEW(KRGT), DAVIEW(KBOT), DAVIEW(KTOP))

      ELSE IF (XISSYM) THEN

         IF (GRABRT ()) RETURN 1
         CALL GRBOX ('L',
     &      DAVIEW(KLFT), DAVIEW(KRGT), DVIEW(KBOT,2), DVIEW(KTOP,2))
         IF (YISCOP) THEN
            IF (GRABRT ()) RETURN 1
            CALL GRBOX ('L',
     &         DAVIEW(KLFT), DAVIEW(KRGT), DVIEW(KBOT,4), DVIEW(KTOP,4))
         END IF

      ELSE IF (YISSYM) THEN

         IF (XISCOP) THEN
            IF (GRABRT ()) RETURN 1
            CALL GRBOX ('L',
     &         DVIEW(KLFT,1), DVIEW(KRGT,1), DAVIEW(KBOT), DAVIEW(KTOP))
         END IF
         IF (GRABRT ()) RETURN 1
         CALL GRBOX ('L',
     &      DVIEW(KLFT,2), DVIEW(KRGT,2), DAVIEW(KBOT), DAVIEW(KTOP))

      ELSE

         DO 100 IVW = 1, NDEFVW (.TRUE.)
            IVIEW = IXVW (.TRUE., IVW)
            IF (GRABRT ()) RETURN 1
            CALL GRBOX ('L', DVIEW(KLFT,IVIEW), DVIEW(KRGT,IVIEW),
     &         DVIEW(KBOT,IVIEW), DVIEW(KTOP,IVIEW))
  100    CONTINUE
      END IF

C   --Flush buffer, so label is complete at this point
      CALL PLTFLU

      RETURN
      END
