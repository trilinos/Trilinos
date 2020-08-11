C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE SPPLT1 (A, NEUTRL, IPTIMS, TIMES, IDTIME,
     &   N, OVER, NCRV, NUMCRV,
     &   NENUM, DIST, PLTVAL, TXLAB, TYLAB, NAMES,
     &   NNE, ISEGEL, NPDON, NPTOT, LIDSP, BLKCOL,
     *  MAPEL, MAPND, *)
C=======================================================================

C   --*** SPPLT1 *** (SPLOT) Plot the curve for a time and a variable
C   --   Modified by John Glick - 11/4/88
C   --   Written by Amy Gilkey - revised 01/22/88
C   --
C   --SPPLT1 plots the curve requested for a single time step and for
C   --a single variable, it assumes that the X axis (and possibly the Y
C   --axis) has been scaled.
C   --
C   --Parameters:
C   --   A - IN - the dynamic memory base array
C   --   NEUTRL - IN  - the type of neutral file to write.
C   --   IPTIMS - IN - the selected time steps
C   --   TIMES - IN - the database times
C   --   IDTIME - IN - the id in the IPTIMS array identifying the
C   --                 timeplane of the time being plotted.
C   --   N - IN - the curve variable number
C   --   OVER - IN - true iff the curve is overlaid
C   --   NCRV - IN - the curve number
C   --   NUMCRV - IN - true iff the curve limits needs to be numbered
C   --   NENUM - IN - the selected node/element numbers
C   --   DIST - IN - the plot distances (compressed)
C   --   PLTVAL - IN - the plot data for this curve (compressed)
C   --   TXLAB, TYLAB - IN - the X and Y axis labels (scratch for neutral file)
C   --   NAMES - IN - the variable names
C   --   NNE - IN - the number of selected nodes/elements with defined values
C   --   ISEGEL - IN - the NENUM indices of the defined elements
C   --      (only if NNE < NNENUM)
C   --   NPDON - IN - the current plot number (only if OVER)
C   --   NPTOT - IN - the total number of plots (only if OVER)
C   --   LIDSP(0:*)  - IN - the indices of the selected variables
C   --          whose values will be displayed on the plot legend.
C   --          LIDSP(0) = the number of variables in the list.
C   --          LIDSP(i) identifies the ith variable in the list.
C   --          If LIDSP(i) > 0, LIDSP(i) is the id of a history variable.
C   --          If LIDSP(i) < 0, -LIDSP(i) is the id of a global variable.
C   --          If LIDSP(i) = 0, TIME is to be displayed on the plot legend
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
C   --   * - return statement if cancel or quit
C   --
C   --Common Variables:
C   --   Uses NNENUM of /SELNE/
C   --   Uses NSPVAR of /SPVARS/
C   --   Uses DOGRID, LINTYP, ISYTYP, LABSID of /XYOPT/
C   --   Uses IYSCAL of /XYLIM/
C   --   Sets YMIN, YMAX of /XYLIM/

      include 'params.blk'
      include 'neutral.blk'
      include 'dbnums.blk'
      include 'xyopt.blk'
      include 'selne.blk'
      include 'spvars.blk'
      include 'xylim.blk'

      DIMENSION A(*)
      INTEGER IPTIMS(*)
      REAL TIMES(*)
      INTEGER IDTIME
      LOGICAL OVER
      LOGICAL NUMCRV
      INTEGER NENUM(NNENUM)
      REAL DIST(NNENUM)
      REAL PLTVAL(NNENUM)
      CHARACTER*(*) TXLAB, TYLAB
      CHARACTER*(*) NAMES(*)
      INTEGER ISEGEL(*)
      INTEGER LIDSP(0:*)
      INTEGER BLKCOL(0:NELBLK)
      INTEGER MAPEL(*), MAPND(*)

      LOGICAL GRABRT
      CHARACTER*80 PLTITL

C   --Calculate min/max if needed

      IF (IYSCAL .EQ. 'CURVE') THEN
         CALL MINMAX (NNE, PLTVAL, YMIN, YMAX)
         IF (NEUTRL .EQ. 0)
     *     CALL EXPMAX (' ', YMIN, YMAX)
      END IF
      IF (OVER) THEN
         IF (IYSCAL .EQ. 'CURVE')
     &      CALL XYAXIS (NCRV, DOGRID, TXLAB, TYLAB, BLKCOL, *120)
      END IF

      IF (NEUTRL .EQ. CSV) THEN
         CALL prterr('CMDREQ','CSV neutral format not yet '//
     $        'supported for SPLOT.')
         return 1
      END IF

      IF (NEUTRL .EQ. 0) THEN

C      --Label plot if needed

  100    CONTINUE
         IF (.NOT. OVER) THEN
            CALL SPLAB (A, 1, IPTIMS(IDTIME), TIMES,
     &         NENUM, N, 1, NUMCRV, NAMES,
     &         TXLAB, TYLAB, LIDSP, BLKCOL, MAPEL, MAPND, *120)
            CALL XYAXIS (0, DOGRID, TXLAB, TYLAB, BLKCOL, *120)
         END IF

         IF (GRABRT()) RETURN 1

         IF (OVER) THEN
            CALL GRCOLR (NCRV)
            CALL GRSYMB (LINTYP, ISYTYP, NCRV)
         ELSE
            CALL GRCOLR (1)
            CALL GRSYMB (LINTYP, ISYTYP, 1)
         END IF

         IF (GRABRT()) RETURN 1

C      --Plot variable against distances

         IEND = 0
  110    CONTINUE
         IF (IEND .LT. NNE) THEN
            IF (NNE .LT. NNENUM) THEN
               CALL SPSEGM (NNE, ISEGEL, ISTART, IEND)
            ELSE
               ISTART = 1
               IEND = NNENUM
            END IF
            NLEFT = IEND - ISTART + 1
            CALL PLTCUR (DIST(ISTART), PLTVAL(ISTART), NLEFT)
            GOTO 110
         END IF

         IF (NUMCRV) THEN
            IF (GRABRT()) RETURN 1
            CALL GRNCRV (LABSID, NCRV, NNE, DIST, PLTVAL,
     &         (LINTYP .EQ. 0))
         END IF

C      --Finish plot

         IF (OVER) CALL PLTFLU
         IF (.NOT. OVER) THEN
C         --Set color in case text is requested
            CALL UGRCOL (0, BLKCOL)
            CALL GRPEND (.TRUE., .TRUE., NPDON, NPTOT, .FALSE.,
     *        *100, *120)
         END IF

      ELSE

C      --Get plot labels

         CALL SPLABN (N, TIMES(IPTIMS(IDTIME)), NENUM, NAMES,
     &      PLTITL, TXLAB, TYLAB, MAPEL, MAPND)

C      --Plot variable against distances

         IF (NEUTRL .EQ. XMGR) THEN
           CALL WRTNEU (NNE, DIST, PLTVAL, PLTITL, TXLAB, TYLAB)
         ELSE IF (NEUTRL .EQ. GRAF) THEN
           CALL GRFNEU (NNE, DIST, PLTVAL, PLTITL, TXLAB, TYLAB)
         ELSE IF (NEUTRL .EQ. CSV .or. NEUTRL .EQ. XML) THEN
           CALL prterr('CMDREQ','CSV neutral format not yet '//
     $           'supported for SPLOT.')
         END IF
      END IF

      RETURN

  120 CONTINUE
      RETURN 1
      END
