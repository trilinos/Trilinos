C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

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
