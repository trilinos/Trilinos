C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE PLOTER (A, CURPRO, NEUTRL,
     &   NAMES, NPTIMS, IPTIMS, TIMES, WHOTIM, NLNKE, ISEVOK,
     &   NENUM, XN, YN, ZN, XE, YE, ZE,
     &   NAMECO, LENF, NLNKF, KLINKF, LENL, KLNSET,
     &   IE2ELB, NEWELB, IELBST,
     &   KNPSUR,
     &   ISSNPS, IDNPS, ISSESS, IDESS, LIDSP, BLKCOL,
     &   IDELB, NAMELB, NAMLEN, MAPEL, MAPND)
C=======================================================================

C   --*** PLOTER *** (BLOT) Main routine for plotting
C   --   Modified by John Glick - 11/29/88
C   --   Written by Amy Gilkey - revised 06/02/88
C   --
C   --PLOTER sets up the graphics and calls the appropriate plotting routine
C   --for the current subprogram.
C   --
C   --Parameters:
C   --   A - IN - the dynamic memory base array
C   --   CURPRO - IN - the current program name
C   --   NEUTRL - IN - the type of neutral file to write.
C   --   NAMES - IN - the variable names
C   --   NPTIMS - IN - the number of selected time steps
C   --   IPTIMS - IN - the selected time steps
C   --   TIMES - IN - the database times
C   --   WHOTIM - IN - true iff whole (versus history) time step
C   --   NLNKE - IN - the number of nodes per element
C   --   ISEVOK - IN - the element block variable truth table;
C   --      variable i of block j exists iff ISEVOK(j,i)
C   --   NENUM - IN - the selected node/element numbers
C   --   XN, YN, ZN - IN - the nodal mesh coordinates
C   --   XE, YE, ZE - IN - the element centroid mesh coordinates
C   --   NAMECO - IN - the coordinate names
C   --   LENF - IN/OUT - the cumulative face counts by element block
C   --   NLNKF - IN - the number of nodes per face
C   --   KLINKF - IN/OUT - the dynamic memory index of the connectivity
C   --      for all faces
C   --   LENL - IN/OUT - the cumulative line counts by element block
C   --   KLNSET - IN/OUT - the dynamic memory index of the sorted line set
C   --   IE2ELB - IN - the element block for each element
C   --   IELBST - IN - the element block status:
C   --      -1 = OFF, 0 = ON, but not selected, 1 = selected
C   --   KNPSUR - IN/OUT - the index of NPSURF -
C   --      the node numbers of the surface nodes or mesh boundary nodes (2D)
C   --   ISSNPS - IN - the indices of the selected node sets
C   --   IDNPS - IN - the node set ID for each set
C   --   ISSESS - IN - the indices of the selected side sets
C   --   IDESS - IN - the side set ID for each set
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
C   --
C   --Common Variables:
C   --   Uses DBORD0 of /LAYOUT/
C   --   Sets DVIEW0 of /LAYOUT/

      PARAMETER (KLFT=1, KRGT=2, KBOT=3, KTOP=4, KNEA=5, KFAR=6)
      include 'params.blk'
      include 'layout.blk'
      include 'dbnums.blk'

      DIMENSION A(*)
      CHARACTER*(*) CURPRO
      CHARACTER*(NAMLEN) NAMES(*)
      INTEGER IPTIMS(*)
      REAL TIMES(*)
      LOGICAL WHOTIM(*)
      INTEGER NLNKE(*)
      LOGICAL ISEVOK(*)
      INTEGER NENUM(*)
      REAL XN(*), YN(*), ZN(*)
      REAL XE(*), YE(*), ZE(*)
      CHARACTER*(*) NAMECO(*)
      INTEGER LENF(*)
      INTEGER NLNKF(*)
      INTEGER LENL(*)
      INTEGER IE2ELB(*)
      INTEGER NEWELB
      INTEGER IELBST(*)
      INTEGER ISSNPS(*)
      INTEGER IDNPS(*)
      INTEGER ISSESS(*)
      INTEGER IDESS(*)
      INTEGER LIDSP(0:*)
      INTEGER BLKCOL(0:NELBLK)
      INTEGER IDELB(*)
      CHARACTER*(*) NAMELB(*)
      INTEGER MAPEL(*), MAPND(*)

C   --Set up graphics

      IF ((CURPRO .EQ. 'MESH') .OR. (CURPRO .EQ. 'DETOUR')
     &   .OR. (CURPRO .EQ. 'PATHLINE')) THEN

C      --Set up graphics for mesh plot
         CALL PRESET (.FALSE., .TRUE., DBORD0, DVIEW0)

      ELSE IF ((CURPRO .EQ. 'TPLOT') .OR. (CURPRO .EQ. 'SPLOT')) THEN

C      --Set up graphics for line plot
         CALL PRESET (.TRUE., .TRUE., DBORD0, DVIEW0)

      ELSE
         CALL PRTERR ('PROGRAM', 'Undefined subprogram')
         GOTO 100
      END IF

C   --Use the standard color table
      CALL GRCOLU ('STANDARD')
C   --Wipe out the mesh plot pick
      CALL INPICK ('NONE')

C   --Call plotting routine for current subprogram

      IF (CURPRO .EQ. 'MESH') THEN
         CALL STCLST('MESH')
         CALL MSMAIN (A, NAMECO, NAMES,
     &      LENF, NLNKE, NLNKF, KLINKF, LENL, KLNSET,
     &      NEWELB, IELBST,
     &      KNPSUR,
     &      ISSNPS, IDNPS, ISSESS, IDESS, NENUM, LIDSP,
     &      BLKCOL, IDELB, NAMELB, MAPEL, MAPND)

      ELSE IF (CURPRO .EQ. 'DETOUR') THEN
         CALL STCLST('DETOUR')
         CALL DTMAIN (A, NAMECO, NAMES, NPTIMS, IPTIMS, TIMES,
     &      LENF, NLNKE, NLNKF, KLINKF, LENL, KLNSET,
     &      NEWELB, IELBST,
     &      KNPSUR, ISEVOK,
     &      ISSNPS, IDNPS, ISSESS, IDESS, LIDSP, BLKCOL,
     &      IDELB,  NAMELB, MAPEL, MAPND)

      ELSE IF (CURPRO .EQ. 'PATHLINE') THEN
         CALL STCLST('PATH')
         CALL LNMAIN (A, NAMECO, NAMES, NPTIMS, IPTIMS, TIMES, WHOTIM,
     &      LENF, NLNKE, NLNKF, KLINKF, LENL, KLNSET,
     &      NEWELB, IELBST,
     &      KNPSUR,
     &      ISSNPS, IDNPS, ISSESS, IDESS, LIDSP, BLKCOL,
     &      IDELB, NAMELB, MAPEL, MAPND)

      ELSE IF (CURPRO .EQ. 'TPLOT') THEN
         CALL STCLST('TPLOT')
         CALL TPMAIN (A, NEUTRL, NAMES,
     &      NPTIMS, IPTIMS, TIMES, WHOTIM, BLKCOL,
     &      IDELB, MAPEL, MAPND)

      ELSE IF (CURPRO .EQ. 'SPLOT') THEN
         CALL STCLST('SPLOT')
         CALL SPMAIN (A, NEUTRL, NAMES, NPTIMS, IPTIMS, TIMES,
     &      NENUM, XN, YN, ZN, XE, YE, ZE, IE2ELB, ISEVOK,
     &      LIDSP, BLKCOL, IDELB, MAPEL, MAPND)
      END IF

  100 CONTINUE
      RETURN
      END
