C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE SPLABN (IPVAR, TIME, NENUM, NAMES,
     &   PLTITL, TXLAB, TYLAB, MAPEL, MAPND)
C=======================================================================

C   --*** SPLABN *** (SPLOT) Get neutral file plot labels
C   --   Written by Amy Gilkey - revised 03/10/86
C   --
C   --SPLABN makes up the plot titles and labels for the neutral file.
C   --
C   --Parameters:
C   --   IPVAR - IN - the /SPVARS/ index of the starting plot variable
C   --   TIME - IN - the plot time
C   --   NENUM - IN - the node/element numbers
C   --   NAMES - IN - the variable names
C   --   PLTITL - OUT - the plot title describing the curves to be
C   --      plotted (e.g. "TIME vs SIGXX at ELEMENT 30")
C   --   TXLAB, TYLAB - OUT - the X and Y axis labels, either the
C   --      user-input labels or the plot variable descriptions
C   --
C   --Common Variables:
C   --   Uses NODVAR, NNENUM of /SELNE/
C   --   Uses ISVID of /SPVARS/
C   --   Uses XLAB, YLAB of /XYLAB/

      include 'params.blk'
      include 'dbnums.blk'
      include 'selne.blk'
      include 'spvars.blk'
      include 'xylab.blk'

      INTEGER NENUM(NNENUM)
      CHARACTER*(*) NAMES(*)
      CHARACTER*(*) PLTITL
      CHARACTER*(*) TXLAB, TYLAB
      INTEGER MAPEL(*), MAPND(*)

      CHARACTER*32 STRNUM
      CHARACTER*32 STRTIM
      CHARACTER*(MXNAME) NAM

C   --Get the plot legend

      NAM = NAMES(ISVID(IPVAR))

      if (nodvar) then
        WRITE (STRNUM, 10000, IOSTAT=IDUM)
     *    MAPND(NENUM(1)), MAPND(NENUM(NNENUM))
      else
        WRITE (STRNUM, 10000, IOSTAT=IDUM)
     *    MAPEL(NENUM(1)), MAPEL(NENUM(NNENUM))
      end if
10000  FORMAT (I9, '..', I9)
      CALL PCKSTR (1, STRNUM)

      CALL NUMSTR1(4, TIME, STRTIM, LSTR)

      IF (NODVAR) THEN
         PLTITL = 'DISTANCE vs ' // NAM(:LENSTR(NAM))
     &      // ' NODES ' // STRNUM(:LENSTR(STRNUM))
     &      // ' at TIME ' // STRTIM(:LSTR)
      ELSE
         PLTITL = 'DISTANCE vs ' // NAM(:LENSTR(NAM))
     &      // ' ELEMENTS ' // STRNUM(:LENSTR(STRNUM))
     &      // ' at TIME ' // STRTIM(:LSTR)
      END IF

C   --Get the axis labels

      IF (XLAB .NE. ' ') THEN
         TXLAB = XLAB
      ELSE
         TXLAB = 'DISTANCE'
      END IF

      IF (YLAB .NE. ' ') THEN
         TYLAB = YLAB
      ELSE
         TYLAB = NAM
      END IF

      RETURN
      END
