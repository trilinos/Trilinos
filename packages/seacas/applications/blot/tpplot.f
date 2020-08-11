C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE TPPLOT (NEUTRL, MAXPTS, NPTS, TIMLIM, PLTVAL,
     &  NAMES, BLKCOL, MAPEL, MAPND)
C=======================================================================

C   --*** TPPLOT *** (TPLOT) Plot the curves
C   --   Written by Amy Gilkey - revised 01/22/88
C   --
C   --TPPLOT does all the plotting for a set, including labeling.
C   --It also calculates the scaling information for each curve.
C   --
C   --Parameters:
C   --   NEUTRL - IN  - the type of neutral file to write.
C   --   MAXPTS - IN - the maximum number of points on a curve
C   --   NPTS - IN - the number of points on each curve
C   --   TIMLIM - IN - the starting and ending time
C   --   PLTVAL - IN - the plot data;
C   --      PLTVAL(x,NTPVAR+1) holds the times if TIMPLT
C   --      PLTVAL(x,NTPVAR+2) holds the compressed times if TIMPLT and needed
C   --   NAMES - IN - the variable names
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
C   --
C   --Common Variables:
C   --   Uses NTPCRV, NTPVAR, TIMPLT of /TPVARS/
C   --   Uses DOGRID, LINTYP, ISYTYP, OVERLY of /XYOPT/
C   --   Sets XMIN, XMAX, YMIN, YMAX of /XYLIM/

      PARAMETER (NUMSYM = 6, NUMLIN = 6)

      include 'params.blk'
      include 'neutral.blk'
      include 'dbnums.blk'
      include 'tpvars.blk'
      include 'xyopt.blk'
      include 'xylim.blk'

      INTEGER NPTS(NTPVAR)
      REAL TIMLIM(2)
      REAL PLTVAL(MAXPTS,NTPVAR+2)
      CHARACTER*(*) NAMES(*)
      INTEGER BLKCOL(0:NELBLK)
      INTEGER MAPEL(*), MAPND(*)

      LOGICAL GRABRT, gobck
      LOGICAL SCACRV
      LOGICAL NUMCRV
      LOGICAL DOLEG
      CHARACTER*1024 PLTITL
      CHARACTER*1024 TXLAB, TYLAB

      LOGICAL SVOVER
      CHARACTER*8 SVLSID

C   --Save user-set parameters and set program parameters (to eliminate
C   --checks for senseless conditions)

      SVOVER = OVERLY
      SVLSID = LABSID

      IF ((NEUTRL .NE. 0) .OR. (NTPCRV .LE. 1)) OVERLY = .FALSE.

      IF (IXSCAL .NE. 'SET') THEN
        IXSCAL = IAXSCA
        IF (TIMPLT) THEN
          IXSCAL = 'ALL'
        ELSE IF (NEUTRL .NE. 0) THEN
          IXSCAL = 'CURVE'
        ELSE IF (NTPCRV .EQ. 1) THEN
          IXSCAL = 'ALL'
        ELSE IF (OVERLY) THEN
          IF (IXSCAL .EQ. 'PLOT') IXSCAL = 'ALL'
        ELSE
          IF (IXSCAL .EQ. 'CURVE') IXSCAL = 'PLOT'
        END IF
      END IF
      IF (IYSCAL .NE. 'SET') THEN
        IYSCAL = IAXSCA
        IF (NEUTRL .NE. 0) THEN
          IYSCAL = 'CURVE'
        ELSE IF (NTPCRV .EQ. 1) THEN
          IYSCAL = 'ALL'
        ELSE IF (OVERLY) THEN
          IF (IYSCAL .EQ. 'PLOT') IYSCAL = 'ALL'
        ELSE
          IF (IYSCAL .EQ. 'CURVE') IYSCAL = 'PLOT'
        END IF
      END IF

      IF ((.NOT. OVERLY)
     &  .OR. ((ISYTYP .LT. 0) .AND. (NUMSYM .GE. NTPCRV))
     &  .OR. ((LINTYP .LT. 0) .AND. (NUMLIN .GE. NTPCRV)))
     &  LABSID = 'NONE'

      NUMCRV = (LABSID .NE. 'NONE')

C   --Calculate axis limits if same scale

      IF (IXSCAL .EQ. 'ALL') THEN
        IF (TIMPLT) THEN
          CALL MINMAX (MAXPTS, PLTVAL(1,NTPVAR+1), XMIN, XMAX)
        ELSE
          CALL CRVLIM ('X', TIMPLT, MAXPTS, NPTS, 1, NTPVAR, PLTVAL)
        END IF
        IF (NEUTRL .EQ. 0)
     *    CALL EXPMAX (LABSID, XMIN, XMAX)
      END IF

      IF (IYSCAL .EQ. 'ALL') THEN
        CALL CRVLIM ('Y', TIMPLT, MAXPTS, NPTS, 1, NTPVAR, PLTVAL)
        IF (NEUTRL .EQ. 0)
     *    CALL EXPMAX (' ', YMIN, YMAX)
      END IF

      SCACRV = (IXSCAL .EQ. 'CURVE') .OR. (IYSCAL .EQ. 'CURVE')

C   --Label plot if overlaid

 100  CONTINUE
      IF (OVERLY) THEN
        CALL TPLAB (1, NTPCRV, NUMCRV, TIMLIM, NAMES,
     &    TXLAB, TYLAB, BLKCOL, MAPEL, MAPND, *130)
        IF (.NOT. SCACRV)
     &    CALL XYAXIS (0, DOGRID, TXLAB, TYLAB, BLKCOL, *130)
      END IF

      gobck = .false.
      if (neutrl .ne. csv .and. neutrl .ne. raw) then
        N = 1
        np = 1
 120    continue

          IF (TIMPLT) THEN
            IF (NPTS(N) .EQ. MAXPTS) THEN
              NX = NTPVAR+1
            ELSE
              NX = NTPVAR+2
            END IF
            NY = N
          ELSE
            NX = N
            NY = N+1
          END IF

C      --Calculate min/max if needed

          IF ((IXSCAL .EQ. 'PLOT') .OR. (IXSCAL .EQ. 'CURVE')) THEN
            CALL CRVLIM ('X', TIMPLT, MAXPTS, NPTS, N, NY, PLTVAL)
            IF (NEUTRL .EQ. 0)
     *        CALL EXPMAX (LABSID, XMIN, XMAX)
          END IF
          IF ((IYSCAL .EQ. 'PLOT') .OR. (IYSCAL .EQ. 'CURVE')) THEN
            CALL CRVLIM ('Y', TIMPLT, MAXPTS, NPTS, N, NY, PLTVAL)
            IF (NEUTRL .EQ. 0)
     *        CALL EXPMAX (' ', YMIN, YMAX)
          END IF
          IF (OVERLY) THEN
            IF (SCACRV)
     &        CALL XYAXIS (NP, DOGRID, TXLAB, TYLAB, BLKCOL, *130)
          END IF

          IF (NEUTRL .EQ. 0) THEN

C         --Label plot if needed

 110        CONTINUE
            IF (.NOT. OVERLY) THEN
              CALL TPLAB (N, 1, NUMCRV, TIMLIM, NAMES,
     &          TXLAB, TYLAB, BLKCOL, MAPEL, MAPND, *130)
              CALL XYAXIS (0, DOGRID, TXLAB, TYLAB, BLKCOL, *130)
            END IF

            IF (GRABRT()) GOTO 130
            IF (OVERLY) THEN
              CALL GRCOLR (NP)
              CALL GRSYMB (LINTYP, ISYTYP, NP)
            ELSE
              CALL GRCOLR (1)
              CALL GRSYMB (LINTYP, ISYTYP, 1)
            END IF

C         --Plot variable against time or variable against variable

            IF (GRABRT()) GOTO 130
            CALL PLTCUR (PLTVAL(1,NX), PLTVAL(1,NY), NPTS(N))

            IF (NUMCRV) THEN
              IF (GRABRT()) GOTO 130
              CALL GRNCRV (LABSID, NP, NPTS(N),

     &          PLTVAL(1,NX), PLTVAL(1,NY), (LINTYP .EQ. 0))
            END IF

C         --Finish plot

            IF (OVERLY) THEN
              CALL PLTFLU
              gobck = .false.
            END IF
            IF (.NOT. OVERLY) THEN
C            --Set color in case text is requested
              CALL UGRCOL (0, BLKCOL)
              gobck = .true.
              CALL GRPEND (.TRUE., .TRUE., NP, NTPCRV, GOBCK,
     $             *110, *130)
            END IF

          ELSE

            gobck = .false.

C         --Get plot labels

            CALL TPLABN (N, TIMLIM, NAMES, PLTITL, TXLAB, TYLAB,
     *        MAPEL, MAPND)

C         --Plot variable against time or variable against variable

            IF (NEUTRL .EQ. XMGR) THEN
              CALL WRTNEU (NPTS(N), PLTVAL(1,NX), PLTVAL(1,NY),
     &          PLTITL, TXLAB, TYLAB)
            ELSE IF (NEUTRL .EQ. GRAF) THEN
              CALL GRFNEU (NPTS(N), PLTVAL(1,NX), PLTVAL(1,NY),
     &          PLTITL, TXLAB, TYLAB)
            END IF
          END IF

          if (gobck) then
             n = ny -1
             np = np - 1
             if (n .lt. 1) n = 1
             if (np .lt. 1) np = 1
          else
             N = NY+1
             np = np + 1
          end if
          if (np .le. ntpcrv) go to 120

      else
        if (neutrl .eq. csv) then
          doleg = .true.
        else
          doleg = .false.
        end if
C ... CSV format neutral file selected...
        CALL TPLABN (1, TIMLIM, NAMES, PLTITL, TXLAB, TYLAB,
     *    MAPEL, MAPND)
        call wrtcsv(ntpcrv, maxpts, npts, pltval, txlab, names, doleg,
     *    MAPEL, MAPND)
      end if
C   --Finish overlaid plot

      IF (OVERLY) THEN
C      --Set color in case text is requested
        CALL UGRCOL (0, BLKCOL)
        CALL GRPEND (.TRUE., .TRUE., 0, 0, .FALSE., *100, *130)
      END IF

 130  CONTINUE

C   --Restore user-set parameters
      OVERLY = SVOVER
      LABSID = SVLSID
      IF (IXSCAL .NE. 'SET') IXSCAL = IAXSCA
      IF (IYSCAL .NE. 'SET') IYSCAL = IAXSCA

      RETURN
      END
