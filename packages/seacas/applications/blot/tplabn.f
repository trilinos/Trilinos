C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE TPLABN (IPVAR, TIMLIM, NAMES, PLTITL, TXLAB, TYLAB,
     *  MAPEL, MAPND)
C=======================================================================

C   --*** TPLABN *** (TPLOT) Get neutral file plot labels
C   --   Written by Amy Gilkey - revised 03/23/87
C   --
C   --TPLABN makes up the plot titles and labels for the neutral file.
C   --
C   --Parameters:
C   --   IPVAR - IN - the /TPVARS/ index of the starting plot variable
C   --   TIMLIM - IN - the starting and ending times for a
C   --      variable-versus-variable curve
C   --   NAMES - IN - the variable names
C   --   PLTITL - OUT - the plot title describing the curves to be
C   --      plotted (e.g. "TIME vs SIGXX at ELEMENT 30" or
C   --      "LOAD vs SIGXX at ELEMENT 30 for times 0.000 to 15.000")
C   --   TXLAB, TYLAB - OUT - the X and Y axis labels, either the
C   --      user-input labels or the plot variable descriptions
C   --
C   --Common Variables:
C   --   Uses TIMPLT, ITVID, ITVNE of /TPVARS/
C   --   Uses XLAB, YLAB of /XYLAB/

      include 'params.blk'
      include 'tpvars.blk'
      include 'xylab.blk'

      REAL TIMLIM(2)
      CHARACTER*(*) NAMES(*)
      CHARACTER*(*) PLTITL
      CHARACTER*(*) TXLAB, TYLAB
      INTEGER MAPEL(*), MAPND(*)

      CHARACTER*(1024) PV1, PV2
      CHARACTER*20 RSTR(2)

C   --Get the plot legend

      N = IPVAR

      IF (TIMPLT) THEN
         PV1 = 'TIME'
         CALL TPLABV (-1, ITVID(N), NAMES(ITVID(N)), ITVNE(N), PV2,
     *     MAPEL, MAPND)
         PLTITL = PV1(:LENSTR(PV1)) // ' vs ' // PV2(:LENSTR(PV2))
         write (*,*) pltitl(:lenstr(pltitl))
      ELSE
         CALL TPLABV (-1, ITVID(N), NAMES(ITVID(N)), ITVNE(N), PV1,
     *    MAPEL, MAPND)
         N = N + 1
         CALL TPLABV (-1, ITVID(N), NAMES(ITVID(N)), ITVNE(N), PV2,
     *     MAPEL, MAPND)
         CALL NUMSTR (2, 4, TIMLIM, RSTR, LSTR)
         PLTITL = PV1(:LENSTR(PV1)) // ' vs ' // PV2(:LENSTR(PV2))
     &      // ' for times ' // RSTR(1)(:LENSTR(RSTR(1)))
     &      // ' to ' // RSTR(2)(:LSTR)
      END IF

C   --Get the axis labels

      IF (XLAB .NE. ' ') THEN
         TXLAB = XLAB
      ELSE
         TXLAB = PV1
      END IF

      IF (YLAB .NE. ' ') THEN
         TYLAB = YLAB
      ELSE
         TYLAB = PV2
      END IF

      RETURN
      END
