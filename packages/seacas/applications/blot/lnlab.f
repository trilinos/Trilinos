C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE LNLAB (A, NEWSET, NPTIMS, IPTIMS, TIMES,
     &   NAMECO, NAMES, IELBST,
     &   ISSNPS, IDNPS, ISSESS, IDESS, ISCNPS, ISCESS,
     &   LIDSP, BLKCOL, *)
C=======================================================================

C   --*** LNLAB *** (PATHLN) Label plot
C   --   Modified by John Glick - 11/9/88
C   --   Written by Amy Gilkey - revised 05/31/88
C   --
C   --LNLAB calls MSLAB to draw the standard mesh plot label, then adds
C   --PATHLN-specific labeling.  MSLAB also draws the axes.
C   --
C   --Parameters:
C   --   A - IN - the dynamic memory base array
C   --   NEWSET - IN - true iff starting a new plot set
C   --   NPTIMS - IN - the number of times on this plot; 0 for no time
C   --   IPTIMS - IN - the plot time indices (starting with current plot)
C   --   TIMES - IN - the database times
C   --   NAMECO - IN - the coordinate names
C   --   NAMES - IN - the variable names
C   --   IELBST - IN - the element block status:
C   --      -1 = OFF, 0 = ON, but not selected, 1 = selected
C   --   ISSNPS - IN - the indices of the selected node sets
C   --   IDNPS - IN - the node set IDs
C   --   ISSESS - IN - the indices of the selected side sets
C   --   IDESS - IN - the side set IDs
C   --   ISCNPS - IN/OUT - size = NUMNPS, set iff NEWSET
C   --   ISCESS - IN/OUT - size = NUMESS, set iff NEWSET
C   --   LIDSP(0:*)  - IN - the indices of the selected variables
C   --          whose values will be displayed on the plot legend.
C   --          LIDSP(0) = the number of variables in the list.
C   --          LIDSP(i) identifies the ith variable in the list.
C   --          If LIDSP(i) > 0, LIDSP(i) is the id of a history variable.
C   --          If LIDSP(i) < 0, -LIDSP(i) is the id of a global variable.
C   --          If LIDSP(i) = 0, TIME is to be displayed on the plot legend.
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
C   --   Uses NDIM, NELBLK of /DBNUMS/
C   --   Uses IS3DIM of /D3NUMS/
C   --   Uses DOLEG of /LEGOPT/
C   --   Uses CHLSIZ of /LAYOUT/

      PARAMETER (KLFT=1, KRGT=2, KBOT=3, KTOP=4, KNEA=5, KFAR=6)

      include 'dbnums.blk'
      include 'dbnumgq.blk'
      include 'd3nums.blk'
      include 'legopt.blk'
      include 'lnvars.blk'
      include 'layout.blk'

      common /debugc/ cdebug
      common /debugn/ idebug
      character*8 cdebug

      DIMENSION A(*)
      LOGICAL NEWSET
      INTEGER IPTIMS(*)
      REAL TIMES(*)
      CHARACTER*(*) NAMECO(*)
      CHARACTER*(*) NAMES(*)
      INTEGER IELBST(NELBLK)
      INTEGER ISSNPS(NUMNPS,4)
      INTEGER ISSESS(NUMESS,4)
      INTEGER IDNPS(*)
      INTEGER IDESS(*)
      INTEGER ISCNPS(*)
      INTEGER ISCESS(*)
      INTEGER LIDSP(0:*)
      INTEGER BLKCOL(0:NELBLK)

      LOGICAL SOFTCH

      REAL DLEGND(KTOP)
      CHARACTER*80 STRING
      REAL RNUM(2)
      CHARACTER*20 RSTR(2)
      CHARACTER TYP

      REAL DYCRV, DYTIMS
      SAVE DYCRV, DYTIMS

C *** Standard plot labeling ***

      CALL MSLAB (A, NEWSET, 0, IDUM, RDUM,
     &   NAMECO, NAMES, IELBST,
     &   ISSNPS, IDNPS, ISSESS, IDESS, ISCNPS, ISCESS,
     &   DLEGND, LIDSP, BLKCOL, *130)
      CALL OUTLIN (BLKCOL, *130)

C *** Plot legend ***

      IF (.NOT. DOLEG(1)) GOTO 120

C   --Set up labeling

      IF (NEWSET) THEN

C      --Get the number of lines needed
         NCRVLN = 0
         DO 100 NP = 1, NLNCRV
            CALL DBVTYP_BL (ILVID(1,NP), TYP, IDUM)
            IF (TYP .EQ. 'H') THEN
               NCRVLN = NCRVLN + 1
            ELSE
               NCRVLN = NCRVLN + 2
            END IF
  100    CONTINUE

         NCENLN = NCRVLN + 1 + 2

         DYLTOP = DLEGND(KTOP) - 1.5*CHLSIZ
         DYLBOT = DLEGND(KBOT) + 1.5*CHLSIZ
         CALL GRYCEN (CHLSIZ, DYLTOP, DYLBOT, NCENLN, NOVER)
         DYCRV = DYLTOP
         DYLTOP = DYCRV - (NCRVLN+1) * CHLSIZ
         DYTIMS = DYLTOP
         DYLTOP = DYTIMS - (2+1) * CHLSIZ
      END IF

C   --Get software character flag for current device
      CALL GRGPARD ('SOFTCHAR', 0, SOFTCH, STRING)

C   --Display plot item(s) (variable names and number) for each pathline

      DY = DYCRV
      DO 110 NP = 1, NLNCRV
         CALL DBVTYP_BL (ILVID(1,NP), TYP, IDUM)
         WRITE (STRING, '(4 (A, 1X))')
     &      (NAMES(ILVID(I,NP)), I=1,NDIM)
         CALL SQZSTR (STRING, LSTR)
         CALL GRTEXT (DLEGND(KLFT), DY, STRING(:LSTR))
         DY = DY - CHLSIZ
         IF ((TYP .EQ. 'N') .OR. (TYP .EQ. 'E')) THEN
            CALL INTSTR (1, 0, ILVNE(NP), STRING, LSTR)
            IF (TYP .EQ. 'N') THEN
               CALL GRTEXT (DLEGND(KLFT), DY,
     &            'Node ' // STRING(:LSTR))
            ELSE
               CALL GRTEXT (DLEGND(KLFT), DY,
     &            'Element ' // STRING(:LSTR))
            END IF
            DY = DY - CHLSIZ
         END IF
  110 CONTINUE

      CALL UGRCOL (0, BLKCOL)
      DY = DY - CHLSIZ

C   --Display times for pathlines

      RNUM(1) = TIMES(IPTIMS(1))
      RNUM(2) = TIMES(IPTIMS(NPTIMS))
      CALL NUMSTR (2, 4, RNUM, RSTR, LSTR)
      CALL PCKSTR (2, RSTR)
      CALL GRTEXT (DLEGND(KLFT), DY, 'TIMES ' // RSTR(1))
      DY = DY - CHLSIZ
      CALL GRTEXT (DLEGND(KLFT), DY, '  TO ' // RSTR(2))
      DY = DY - CHLSIZ

  120 CONTINUE

C   --Flush buffer, so label is complete at this point
      CALL PLTFLU

      RETURN

  130 CONTINUE
      RETURN 1
      END
