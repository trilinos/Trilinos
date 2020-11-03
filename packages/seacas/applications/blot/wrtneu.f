C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE WRTNEU (NPTS, XPTS, YPTS, PLTITL, TXLAB, TYLAB)
C=======================================================================

C   --*** WRTNEU *** (XYPLOT) Write curve to neutral file
C   --   Written by Amy Gilkey - revised 04/21/88
C   --
C   --WRTNEU writes the data for a curve to a neutral file which is
C   --readable by the xmgr program.  The first
C   --time the routine is called, the neutral file is opened.
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

      CHARACTER*2048 filnam, errmsg
      REAL XPTS(NPTS), YPTS(NPTS)
      CHARACTER*(*) PLTITL
      CHARACTER*(*) TXLAB, TYLAB

      CHARACTER*15 CURVE
      integer numgrf
      save numgrf

      IF (.NOT. NEUOPN) THEN
        filnam = basenam(:lenstr(basenam)) // '.xmgr'

C      --Open the neutral file and write the title line
        write (*,*) "Neutral File: ", filnam
        open (unit=neu, file=filnam(:lenstr(filnam)), form='formatted',
     *    status='unknown', iostat=ierr)
        IF (IERR .NE. 0) THEN
          ERRMSG = 'Neutral file "'//FILNAM(:LENSTR(FILNAM))//
     *      '" could not be opened.'
          CALL PRTERR ('CMDERR', ERRMSG(:LENSTR(ERRMSG)))
          GOTO 170
        END IF
        NEUOPN = .TRUE.

        WRITE (NEU, 10000) 'Written by:  ',
     *    DRAW(1), DRAW(3)(:LENSTR(DRAW(3))), DRAW(4)(:LENSTR(DRAW(4)))
        WRITE (NEU, 10000) 'Created by:  ',
     *    CREATE(1), CREATE(3)(:LENSTR(CREATE(3))),
     *    CREATE(4)(:LENSTR(CREATE(4)))
        WRITE (NEU, 10000) 'Modified by: ',
     &    MODIFY(1), MODIFY(3)(:LENSTR(MODIFY(3))),
     *    MODIFY(4)(:LENSTR(MODIFY(4)))
10000   FORMAT ('# ',A, A, A, ' ', A, ' ', A)

C--Write header information.
        numgrf = 0
        numcrv = 31
C... Set world min/max values to appropriate values.
        WXMIN =  1.0E30
        WXMAX = -1.0E30
        WYMIN =  1.0E30
        WYMAX = -1.0E30
      END IF
C   --If numcrv >= 30, increment grf number and print header.
      if (numcrv .ge. 30 .or.
     *  (inccrv .lt. 0 .and. numcrv .gt. abs(inccrv))) then
        numcrv = 1
        write (curve, '(A1, I7)') 'g', numgrf
        CALL PCKSTR (1, CURVE)
        WRITE (NEU, 10100) curve
10100   FORMAT ('@ with ',A,/
     *    '@ legend on'/
     *    '@ legend loctype view'/
     *    '@ legend char size 0.75')
        numgrf = numgrf + 1
      end if
C   --Get the curve name
      WRITE (CURVE, '(A1, I7)') 's', NUMCRV-1
      CALL PCKSTR (1, CURVE)

      WRITE (*, 10090) PLTITL(:LENSTR(PLTITL)), numgrf-1, numcrv-1
10090 FORMAT (' Writing "',A,'" to Graph ',i2,', Set ',i2)

C   --Write the begin curve record with the curve name
      WRITE (NEU, 10020) 'title', TITLE(:lenstr(title))
      WRITE (NEU, 10020) 'subtitle', PLTITL(:lenstr(pltitl))

      ncol = mod(numcrv-1,15) + 1
      WRITE (NEU, 10010) CURVE, 'color', ncol
10010 FORMAT ('@ ',A, A, I7)

C   --Write the title lines

10020 FORMAT ('@ ',A,' "', A,'"')

C   --Write the X and Y labels

      WRITE (NEU, 10040) '@ xaxis label', TXLAB(:lenstr(txlab))
      WRITE (NEU, 10040) '@ yaxis label', TYLAB(:lenstr(tylab))
10040 FORMAT (A,' "',A,'"')
      write (neu, 10050) numcrv-1, TYLAB(:lenstr(tylab))
10050 FORMAT ('@ legend string ',i5,' "',A,'"')
      write (neu, 10060) curve, TYLAB(:lenstr(tylab))
10060 FORMAT ('@ ',A,' comment "',A,'"')
C   --Write the min/max, the number of points, and the auxiliary data flag

C   --Write the data points

      DO 160 I = 1, NPTS
        WXMAX = MAX(WXMAX, XPTS(I))
        WXMIN = MIN(WXMIN, XPTS(I))
        WYMAX = MAX(WYMAX, YPTS(I))
        WYMIN = MIN(WYMIN, YPTS(I))
        WRITE (NEU, 10070) XPTS(I), YPTS(I)
10070   FORMAT (1PE15.7E3,4x,1pe15.7E3)
 160  CONTINUE
      WRITE (NEU, '(A)') '&'
 170  CONTINUE
      NUMCRV = NUMCRV + 1
      RETURN
      END

