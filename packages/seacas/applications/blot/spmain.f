C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE SPMAIN (A, NEUTRL, NAMES, NPTIMS, IPTIMS,
     *  TIMES, NENUM,
     &   XN, YN, ZN, XE, YE, ZE,
     &   IE2ELB, ISEVOK, LIDSP, BLKCOL, IDELB,
     *  MAPEL, MAPND)
C=======================================================================

C   --*** SPMAIN *** (SPLOT) SPLOT main plot routine
C   --   Modified by John Glick - 11/9/88
C   --   Written by Amy Gilkey - revised 03/07/88
C   --
C   --SPMAIN calculates the distances, reads the variables, and displays
C   --the distance-versus-variable plots.
C   --
C   --Parameters:
C   --   A          - IN  - the dynamic memory base array
C   --   NEUTRL     - IN  - the type of neutral file to write.
C   --   NAMES      - IN  - the variable names
C   --   NPTIMS     - IN  - the number of selected time steps
C   --   IPTIMS     - IN  - the selected time steps
C   --   TIMES      - IN  - the database times
C   --   NENUM      - IN  - the selected node/element numbers
C   --   XN,YN,ZN   - IN  - the nodal mesh coordinates
C   --   XE,YE,ZE   - IN  - the element centroid mesh coordinates
C   --   IE2ELB     - IN  - the element block for each element
C   --   ISEVOK     - IN  - the element block variable truth table;
C   --                      variable i of block j exists iff ISEVOK(j,i)
C   --   LIDSP(0:*) - IN  - the indices of the selected variables
C   --                      whose values will be displayed on the plot legend.
C   --                      LIDSP(0) = the number of variables in the list.
C   --                      LIDSP(i) identifies the ith variable in the list.
C   --                      LIDSP(i) > 0, LIDSP(i) is the id of a history var.
C   --                      LIDSP(i) < 0, -LIDSP(i) is the id of a global var.
C   --                      LIDSP(i) = 0, TIME is displayed on the plot legend.
C   --   BLKCOL     - I/O - the user selected colors of the element blocks.
C   --                      BLKCOL(0) = 1 if the user defined material
C   --                                  colors should be used in mesh plots.
C   --                                = -1 if program selected colors should
C   --                                  be used.
C   --                      BLKCOL(i) = the user selected color of element
C   --                                  block i:
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

      include 'params.blk'

      PARAMETER (KLFT=1, KRGT=2, KBOT=3, KTOP=4)

      include 'dbnums.blk'
      include 'selne.blk'
      include 'spvars.blk'

      DIMENSION A(*)
      CHARACTER*(*) NAMES(*)
      INTEGER IPTIMS(NPTIMS)
      REAL TIMES(*)
      INTEGER NENUM(NNENUM)
      REAL XN(*), YN(*), ZN(*)
      REAL XE(*), YE(*), ZE(*)
      INTEGER IE2ELB(*)
      LOGICAL ISEVOK(*)
      INTEGER LIDSP(0:*)
      INTEGER BLKCOL(0:NELBLK)
      INTEGER IDELB(*)
      INTEGER MAPEL(*), MAPND(*)

C   --Use the selected color table
      CALL GRCOLU ('ALTERNATE')

C   --Reserve storage for plot variables

      LPLVAL = NNENUM * NSPVAR * NPTIMS
      CALL MDRSRV ('PLTVAL', KPLVAL, LPLVAL)
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 100

C   --Read plot variables

      CALL SPREAD (A, NPTIMS, IPTIMS, NENUM, A(KPLVAL))

      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 100

C   --Identify and index element variables with undefined values

      CALL MDRSRV ('IXSEGV', KXSEGV, NSPVAR)
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 100

      CALL SPDSEG (A, NENUM, IE2ELB, ISEVOK,
     &   NSEGV, A(KXSEGV), KSEGEL)

      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 100

C   --Calculate distances

      CALL MDRSRV ('DIST', KDIST, NNENUM)
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 100

      IF (NODVAR) THEN
         CALL CALDIS (NNENUM, NENUM, XN, YN, ZN, A(KDIST))
      ELSE
         CALL CALDIS (NNENUM, NENUM, XE, YE, ZE, A(KDIST))
      END IF

C   --Plot requested curves

      IF (NSEGV .GT. 0) THEN
         CALL MDRSRV ('XDIST', KXDIST, NNENUM)
         CALL MDSTAT (NERR, MEM)
         IF (NERR .GT. 0) GOTO 100
      END IF

      CALL SPPLOT (A, NEUTRL, NPTIMS, IPTIMS, TIMES, NENUM,
     &   A(KDIST), A(KPLVAL), NAMES,
     &   A(KXSEGV), A(KSEGEL), A(KXDIST), LIDSP, BLKCOL,
     *  MAPEL, MAPND)

      CALL MDDEL ('PLTVAL')
      CALL MDDEL ('DIST')
      IF (NSEGV .GT. 0) THEN
         CALL MDDEL ('XDIST')
      END IF
      CALL MDDEL ('IXSEGV')
      IF (NSEGV .GT. 0) THEN
         CALL MDDEL ('ISEGEL')
      END IF

  100 CONTINUE
      RETURN
      END
