C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE MRKFAC (NUMFAC, F, HIDEF, XF, YF, ZF,
     &   NMIN, NMAX, FMIN, FMAX, BLKCOL, *)
C=======================================================================

C   --*** MRKFAC *** (DETOUR) Mark minimum/maximum face values
C   --   Written by Amy Gilkey - revised 03/09/88
C   --
C   --MRKFAC marks the faces which have the minimum or maximum values
C   --with symbols.
C   --
C   --Parameters:
C   --   NUMFAC - IN - the number of surface faces
C   --   F - IN - the face variable values
C   --   HIDEF(i) - IN - true iff face i is hidden (3D only)
C   --   XF, YF, ZF - IN - the face center coordinates
C   --   NMIN, NMAX - IN - the number of variables values matching the
C   --      minimum and the maximum
C   --   FMIN, FMAX - IN - the minimum and maximum function values
C   --   BLKCOL - IN - the user selected colors of the element blocks.
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
C   --   Uses IS3DIM, NNPSUR of /D3NUMS/
C   --   Uses MAXMIN, MAXMAX of /CNTR/

      include 'dbnums.blk'
      include 'd3nums.blk'
      include 'cntr.blk'

      REAL F(*)
      LOGICAL HIDEF(*)
      REAL XF(*), YF(*), ZF(*)
      INTEGER BLKCOL(0:NELBLK)

      LOGICAL GRABRT
      LOGICAL MRKMIN, MRKMAX

      CALL UGRCOL (0, BLKCOL)

      MRKMIN = (MAXMIN .GE. NMIN)
      MRKMAX = (MAXMAX .GE. NMAX)

      IF ((FMIN .NE. FMAX) .AND. (MRKMIN .OR. MRKMAX)) THEN

C      --Mark minimum and maximum values with centered symbols

         DO 100 IFAC = 1, NUMFAC
            IF (IS3DIM) THEN
               IF (HIDEF(IFAC)) GOTO 100
            END IF

            IF (FMIN .EQ. F(IFAC)) THEN
               IF (MRKMIN) THEN
                  IF (GRABRT ()) RETURN 1
#if NeedsDoubleEscape
                  CALL MPD2SY (1, XF(IFAC), YF(IFAC), '\\CCI')
                  CALL MPD2SY (1, XF(IFAC), YF(IFAC), '\\CCS')
#else
                  CALL MPD2SY (1, XF(IFAC), YF(IFAC), '\CCI')
                  CALL MPD2SY (1, XF(IFAC), YF(IFAC), '\CCS')
#endif
               END IF
            ELSE IF (FMAX .EQ. F(IFAC)) THEN
               IF (MRKMAX) THEN
                  IF (GRABRT ()) RETURN 1
#if NeedsDoubleEscape
                  CALL MPD2SY (1, XF(IFAC), YF(IFAC), '\\CX')
                  CALL MPD2SY (1, XF(IFAC), YF(IFAC), '\\CCS')
#else
                  CALL MPD2SY (1, XF(IFAC), YF(IFAC), '\CX')
                  CALL MPD2SY (1, XF(IFAC), YF(IFAC), '\CCS')
#endif
               END IF
            END IF
  100    CONTINUE

         CALL PLTFLU
      END IF

      RETURN
      END
