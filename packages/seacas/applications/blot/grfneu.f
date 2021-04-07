C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE GRFNEU (NPTS, XPTS, YPTS, PLTITL, TXLAB, TYLAB)
C=======================================================================

C   --*** WRTNEU *** (XYPLOT) Write curve to neutral file
C   --   Written by Amy Gilkey - revised 04/21/88
C   --
C   --WRTNEU writes the data for a curve to a neutral file.  The first
C   --time the routine is called, the neutral file is opened.
C   --
C   --The format of the neutral file is described in "GRAFAID Code User
C   --Manual" under Neutral File Format.
C   --
C   --Parameters:
C   --   NPTS - IN - the number of points on the curve
C   --   XPTS, YPTS - IN - the points on the curve
C   --   PLTITL - IN - the plot title describing the curve
C   --      (e.g. "TIME vs SIGXX at ELEMENT 30")
C   --   TXLAB, TYLAB - IN - the X and Y axis labels, either the
C   --      user-input labels or the plot variable descriptions
C   --
C   --Common Variables:
C   --   Uses TITLE, CREATE, MODIFY, DRAW of /DBTITL/
C   --   Uses DOQA, CAPTN of /LEGOPT/
C   --   Uses XMIN, XMAX, YMIN, YMAX of /LIMITS/
C   --   Uses NEU, NUMCRV, INCCRV, CRVNAM of /NEUTR./
C   --   Uses and sets NEUOPN of /NEUTR./

      include 'dbname.blk'
      include 'dbtitl.blk'
      include 'legopt.blk'
      include 'xylim.blk'
      include 'neutr.blk'

      REAL XPTS(NPTS), YPTS(NPTS)
      CHARACTER*(*) PLTITL
      CHARACTER*(*) TXLAB, TYLAB

      CHARACTER*2048 FILNAM, ERRMSG
      CHARACTER*4 XTYP
      CHARACTER*4 AXTYP
      CHARACTER AUX
      CHARACTER*80 CURVE

      DATA AUX /'F'/
      DATA AXTYP /'NOLO'/

      IF (.NOT. GRFOPN) THEN

C      --Open the neutral file and write the title line

        filnam = basenam(:lenstr(basenam)) // '.grf'

C      --Open the neutral file and write the title line
        open (unit=neugrf, file=filnam(:lenstr(filnam)),
     *    form='formatted', status='unknown', iostat=ierr)
        IF (IERR .NE. 0) THEN
          ERRMSG = 'Neutral file "'//FILNAM(:LENSTR(FILNAM))//
     *      '" could not be opened.'
          CALL PRTERR ('CMDERR', ERRMSG(:LENSTR(ERRMSG)))
          GOTO 170
        END IF
        GRFOPN = .TRUE.

        WRITE (NEUGRF, 10000) DRAW(1)(:LENSTR(DRAW(1))),
     &    CREATE(1)(:LENSTR(CREATE(1))), CREATE(3), CREATE(4),
     &    MODIFY(1)(:LENSTR(MODIFY(1))), MODIFY(3), MODIFY(4)
10000   FORMAT (A, ': CREATED BY ', A, ' ', A, ' ', A,
     &    ', MODIFIED BY ', A, ' ', A, ' ', A)
      END IF

C   --Figure out if the X points are monotonic or not

      XTYP = 'MONO'
      DO 100 I = 2, NPTS
         IF (XPTS(I-1) .GE. XPTS(I)) THEN
            XTYP = 'NONM'
            GOTO 110
         END IF
  100 CONTINUE
  110 CONTINUE

C   --Get the curve name

      WRITE (CURVE, '(A32, I7)', IOSTAT=IDUM) CRVNAM, NUMCRV
      CALL PCKSTR (1, CURVE)
      NUMCRV = NUMCRV + INCCRV

      WRITE (*, 10090) 'Writing ', PLTITL(:LENSTR(PLTITL))

C   --Write the begin curve record with the curve name

      WRITE (NEUGRF, 10010) 'BEGIN CURVE', CURVE
10010  FORMAT (A, ',', A15)

C   --Write the title lines

      DO 120 IEND = 3, 1, -1
         IF (CAPTN(IEND,2) .NE. ' ') GOTO 130
  120 CONTINUE
  130 CONTINUE
      IF (DOQA(2)) THEN
         WRITE (NEUGRF, 10020) 1+MAX(IEND,1), TITLE
         IF (IEND .GT. 0) THEN
            DO 140 I = 1, IEND
               WRITE (NEUGRF, 10030) CAPTN(I,2)
  140       CONTINUE
         ELSE
            WRITE (NEUGRF, 10030) PLTITL
         END IF
      ELSE IF (IEND .GT. 0) THEN
         WRITE (NEUGRF, 10020) IEND, CAPTN(1,2)
         DO 150 I = 2, IEND
            WRITE (NEUGRF, 10030) CAPTN(I,2)
  150    CONTINUE
      ELSE
         WRITE (NEUGRF, 10020) 1, PLTITL
      END IF
10020  FORMAT (I1, ',', A80)
10030  FORMAT (A80)

C   --Write the X and Y labels

      WRITE (NEUGRF, 10040) TXLAB
      WRITE (NEUGRF, 10040) TYLAB
10040  FORMAT (A40)

C   --Write the min/max, the number of points, and the auxiliary data flag

      WRITE (NEUGRF, 10050) XMIN, XMAX, YMIN, YMAX, NPTS, AUX
10050  FORMAT (4 (1PE15.7E3, ','), I5, ',', A1)
      WRITE (NEUGRF, 10060) AXTYP, XTYP, ' '
10060  FORMAT (A4, ',', A4, ',', A4)

C   --Write the data points

      DO 160 I = 1, NPTS
         WRITE (NEUGRF, 10070) XPTS(I), YPTS(I)
10070     FORMAT (2 (1PE15.7E3, :, ','))
  160 CONTINUE

C   --Write the end curve record with the curve name

      WRITE (NEUGRF, 10080) 'END CURVE', CURVE
10080  FORMAT (A, ',', A)

  170 CONTINUE
      RETURN
10090  FORMAT (1X, 5A)
      END
