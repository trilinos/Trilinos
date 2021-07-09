C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE SPPLOT (A, NEUTRL, NPTIMS, IPTIMS, TIMES,
     *  NENUM, DIST, PLTVAL, NAMES, IXSEGV, ISEGEL, SQDIST, LIDSP,
     &   BLKCOL, MAPEL, MAPND)
C=======================================================================

C   --*** SPPLOT *** (SPLOT) Plot curves for all times
C   --   Modified by John Glick - 11/9/88
C   --   Written by Amy Gilkey - revised 01/22/88
C   --
C   --SPPLOT generates all the plots requested for the plot set (for
C   --all plot times).  It scales the X axis and the Y axis (as requested).
C   --It then processes each curve by reading the data and calling SPPLT1.
C   --
C   --Parameters:
C   --   A - IN - the dynamic memory base array
C   --   NEUTRL - IN  - the type of neutral file to write.
C   --   NPTIMS - IN - the number of selected times
C   --   IPTIMS - IN - the selected time steps
C   --   TIMES - IN - the database times
C   --   NENUM - IN - the selected node/element numbers
C   --   DIST - IN - the plot distances
C   --   PLTVAL - IN - the plot data
C   --   NAMES - IN - the variable names
C   --   IXSEGV - IN - the index into ISEGEL for "segmented" variables;
C   --      0 if all selected elements are defined
C   --   ISEGEL - IN - the NENUM indices of the defined elements;
C   --      ISEGEL(0,i) = the number of defined elements
C   --   SQDIST - SCRATCH - size = NNENUM
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
C   --
C   --Common Variables:
C   --   Uses NNENUM of /SELNE/
C   --   Uses NSPVAR of /SPVARS/
C   --   Uses DOGRID, ISYTYP, OVERLY, OVERTM of /XYOPT/
C   --   Uses IXSCAL, IYSCAL of /XYLIM/
C   --   Sets XMIN, XMAX, YMIN, YMAX of /XYLIM/

      PARAMETER (NUMSYM = 6, NUMLIN = 6)

      include 'params.blk'
      include 'dbnums.blk'
      include 'selne.blk'
      include 'spvars.blk'
      include 'xyopt.blk'
      include 'xylim.blk'

      DIMENSION A(*)
      INTEGER IPTIMS(NPTIMS)
      REAL TIMES(*)
      INTEGER NENUM(NNENUM)
      REAL DIST(NNENUM)
      REAL PLTVAL(NNENUM,NSPVAR,*)
      CHARACTER*(*) NAMES(*)
      INTEGER IXSEGV(*)
      INTEGER ISEGEL(0:NNENUM,*)
      REAL SQDIST(NNENUM)
      INTEGER LIDSP(0:*)
      INTEGER BLKCOL(0:NELBLK)
      INTEGER MAPEL(*), MAPND(*)

      LOGICAL NUMCRV
      CHARACTER*1024 TXLAB, TYLAB

      LOGICAL SVOVER, SVOVTM
      CHARACTER*8 SVLSID

C   --Save user-set parameters and set program parameters (to eliminate
C   --checks for senseless conditions)

      SVOVTM = OVERTM
      SVOVER = OVERLY
      SVLSID = LABSID

      IF ((NEUTRL .NE. 0) .OR. (NSPVAR .LE. 1)) OVERLY = .FALSE.
      IF ((NEUTRL .NE. 0) .OR. (NPTIMS .LE. 1)) OVERTM = .FALSE.

      IF (IXSCAL .NE. 'SET') THEN
         IXSCAL = 'ALL'
      END IF
      IF (IYSCAL .NE. 'SET') THEN
         IYSCAL = IAXSCA
         IF (NEUTRL .GT. 0) THEN
            IYSCAL = 'CURVE'
         ELSE IF ((NSPVAR .LE. 1) .AND. (NPTIMS .LE. 1)) THEN
            IYSCAL = 'ALL'
         ELSE IF ((OVERLY .OR. (NSPVAR .LE. 1))
     &      .AND. (OVERTM .OR. (NPTIMS .LE. 1))) THEN
            IF (IYSCAL .EQ. 'PLOT') IYSCAL = 'ALL'
         ELSE IF (.NOT. (OVERLY .OR. OVERTM)) THEN
            IF (IYSCAL .EQ. 'CURVE') IYSCAL = 'PLOT'
         END IF
      END IF

      NC = 1
      IF (OVERLY) NC = NPTIMS
      IF (OVERTM) NC = NSPVAR
      IF (((.NOT. OVERLY) .AND. (.NOT. OVERTM))
     &   .OR. ((ISYTYP .LT. 0) .AND. (NUMSYM .GE. NC))
     &   .OR. ((LINTYP .LT. 0) .AND. (NUMLIN .GE. NC)))
     &   LABSID = 'NONE'

      NUMCRV = (LABSID .NE. 'NONE')

C   --Compact variables which have undefined element values (distances
C   --must be compacted for each variable at a later time)

      DO 100 N = 1, NSPVAR
         IF (IXSEGV(N) .GT. 0) THEN
            IX = IXSEGV(N)
            CALL SQZIXV (ISEGEL(0,IX), ISEGEL(1,IX),
     &         PLTVAL(1,N,NPT), PLTVAL(1,N,NPT))
         END IF
  100 CONTINUE

C   --Scale distance parameters (only done once)

      IF (IXSCAL .NE. 'SET') THEN
         CALL MINMAX (NNENUM, DIST, XMIN, XMAX)
         IF (NEUTRL .EQ. 0)
     *     CALL EXPMAX (LABSID, XMIN, XMAX)
      END IF

C   --Scale Y axis if same scale requested

      IF (IYSCAL .EQ. 'ALL') THEN
         DO 120 NPT = 1, NPTIMS
            DO 110 N = 1, NSPVAR
               IF (IXSEGV(N) .LE. 0) THEN
                  NNE = NNENUM
               ELSE
                  NNE = ISEGEL(0,IXSEGV(N))
               END IF
               CALL MINMAX (NNE, PLTVAL(1,N,NPT), VMIN, VMAX)
               IF ((NPT .EQ. 1) .AND. (N .EQ. 1)) THEN
                  YMIN = VMIN
                  YMAX = VMAX
               ELSE
                  YMIN = MIN (YMIN, VMIN)
                  YMAX = MAX (YMAX, VMAX)
               END IF
  110       CONTINUE
  120    CONTINUE
         call adjrng(ymin, ymax)
         CALL EXPMAX (' ', YMIN, YMAX)
      END IF

C   --Set the total number of plots

      IF (OVERTM) THEN
         NPTOT = NSPVAR
      ELSE IF (OVERLY) THEN
         NPTOT = NPTIMS
      ELSE
         NPTOT = NSPVAR * NPTIMS
      END IF

      IF (OVERTM) THEN

         DO 160 N = 1, NSPVAR
            IF (IXSEGV(N) .LE. 0) THEN
               NNE = NNENUM
            ELSE
               NNE = ISEGEL(0,IXSEGV(N))
            END IF

C         --Calculate axis limits if to be done once for all times

            IF (IYSCAL .EQ. 'PLOT') THEN
               DO 130 NPT = 1, NPTIMS
                  CALL MINMAX (NNE, PLTVAL(1,N,NPT), VMIN, VMAX)
                  IF (NPT .EQ. 1) THEN
                     YMIN = VMIN
                     YMAX = VMAX
                  ELSE
                     YMIN = MIN (YMIN, VMIN)
                     YMAX = MAX (YMAX, VMAX)
                  END IF
  130          CONTINUE
               call adjrng(ymin, ymax)
               CALL EXPMAX (' ', YMIN, YMAX)
            END IF

C         --Label plot if overlaid by time

  140       CONTINUE
            CALL SPLAB (A, NPTIMS, IPTIMS, TIMES, NENUM,
     &           N, 1, NUMCRV, NAMES, TXLAB, TYLAB,
     &           LIDSP, BLKCOL,  MAPEL, MAPND, *210)
            IF (IYSCAL .NE. 'CURVE')
     &           CALL XYAXIS (0, DOGRID, TXLAB, TYLAB, BLKCOL, *210)

C         --Compact distances to match compacted variables

            IF (IXSEGV(N) .GT. 0) THEN
               CALL SQZIXV (NNE, ISEGEL(1,IXSEGV(N)), DIST, SQDIST)
            END IF

C         --Plot all time curves for this variable

            DO 150 NPT = 1, NPTIMS
               NPDON = N
               IF (IXSEGV(N) .LE. 0) THEN
                  CALL SPPLT1 (A, NEUTRL, IPTIMS, TIMES, NPT,
     &               N, OVERTM, NPT, NUMCRV,
     &               NENUM, DIST, PLTVAL(1,N,NPT),
     &               TXLAB, TYLAB, NAMES,
     &               NNE, ISEGEL(1,IXSEGV(N)), NPDON,
     &               NPTOT, LIDSP, BLKCOL, MAPEL, MAPND, *210)
               ELSE
                  CALL SPPLT1 (A, NEUTRL, IPTIMS, TIMES, NPT,
     &               N, OVERTM, NPT, NUMCRV,
     &               NENUM, SQDIST, PLTVAL(1,N,NPT),
     &               TXLAB, TYLAB, NAMES,
     &               NNE, ISEGEL(1,IXSEGV(N)), NPDON,
     &               NPTOT, LIDSP, BLKCOL, MAPEL, MAPND, *210)
               END IF
  150       CONTINUE

C         --Finish plot

C     --Set color in case text is requested
            CALL UGRCOL (0, BLKCOL)
            CALL GRPEND (.TRUE., .TRUE., N, NSPVAR, .FALSE., *140, *210)

  160    CONTINUE

      ELSE

         DO 200 NPT = 1, NPTIMS

C         --Calculate axis limits if to be done once for all curves

            IF (IYSCAL .EQ. 'PLOT') THEN
               DO 170 N = 1, NSPVAR
                  IF (IXSEGV(N) .LE. 0) THEN
                     NNE = NNENUM
                  ELSE
                     NNE = ISEGEL(0,IXSEGV(N))
                  END IF
                  CALL MINMAX (NNE, PLTVAL(1,N,NPT), VMIN, VMAX)
                  IF (N .EQ. 1) THEN
                     YMIN = VMIN
                     YMAX = VMAX
                  ELSE
                     YMIN = MIN (YMIN, VMIN)
                     YMAX = MAX (YMAX, VMAX)
                  END IF
  170          CONTINUE
               call adjrng(ymin, ymax)
               CALL EXPMAX (' ', YMIN, YMAX)
            END IF

C         --Label plot if overlaid curves

  180       CONTINUE
            IF (OVERLY) THEN
               CALL SPLAB (A, 1, IPTIMS(NPT), TIMES,
     &            NENUM, 1, NSPVAR,
     &            NUMCRV, NAMES, TXLAB, TYLAB, LIDSP,
     &            BLKCOL,  MAPEL, MAPND, *210)
               IF (IYSCAL .NE. 'CURVE')
     &            CALL XYAXIS (0, DOGRID, TXLAB, TYLAB, BLKCOL, *210)
            END IF

C         --Plot all variable curves for this time

            DO 190 N = 1, NSPVAR
               IF (IXSEGV(N) .LE. 0) THEN
                  NNE = NNENUM
               ELSE
                  NNE = ISEGEL(0,IXSEGV(N))
               END IF

C            --Compact distances to match compacted variables

               IF (IXSEGV(N) .GT. 0) THEN
                  CALL SQZIXV (NNE, ISEGEL(1,IXSEGV(N)), DIST, SQDIST)
               END IF

               IF (OVERLY) THEN
                  NPDON = NPT
               ELSE
                  NPDON = (NPT-1)*NSPVAR + N
               END IF
               IF (IXSEGV(N) .LE. 0) THEN
                  CALL SPPLT1 (A, NEUTRL, IPTIMS, TIMES, NPT,
     &               N, OVERLY, N, NUMCRV,
     &               NENUM, DIST, PLTVAL(1,N,NPT),
     &               TXLAB, TYLAB, NAMES,
     &               NNE, ISEGEL(1,IXSEGV(N)), NPDON,
     &               NPTOT, LIDSP, BLKCOL, MAPEL, MAPND, *210)
               ELSE
                  CALL SPPLT1 (A, NEUTRL, IPTIMS, TIMES, NPT,
     &               N, OVERLY, N, NUMCRV,
     &               NENUM, SQDIST, PLTVAL(1,N,NPT),
     &               TXLAB, TYLAB, NAMES,
     &               NNE, ISEGEL(1,IXSEGV(N)), NPDON,
     &               NPTOT, LIDSP, BLKCOL, MAPEL, MAPND, *210)
               END IF
  190       CONTINUE

C         --Finish plot

            IF (OVERLY) THEN
C            --Set color in case text is requested
               CALL UGRCOL (0, BLKCOL)
               CALL GRPEND (.TRUE., .TRUE., NPT, NPTIMS, .FALSE.,
     $              *180, *210)
            END IF

  200    CONTINUE
      END IF

  210 CONTINUE

C   --Restore user-set parameters
      OVERTM = SVOVTM
      OVERLY = SVOVER
      LABSID = SVLSID
      IF (IXSCAL .NE. 'SET') IXSCAL = IAXSCA
      IF (IYSCAL .NE. 'SET') IYSCAL = IAXSCA

      RETURN
      END

      subroutine adjrng(ymin, ymax)
C     .. Check that YMIN and YMAX are "different enough"
      if ((ymax - ymin) / ymax .lt. 1.0e-4) then
         ymax = ymax * 1.001
         ymin = ymin / 1.001
      end if
      return
      end

