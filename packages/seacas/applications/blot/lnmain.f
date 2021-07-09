C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE LNMAIN (A, NAMECO, NAMES,
     &   NPTIMS, IPTIMS, TIMES, WHOTIM,
     &   LENF, NLNKE, NLNKF, KLINKF, LENL, KLNSET,
     &   NEWELB, IELBST,
     &   KNPSUR,
     &   ISSNPS, IDNPS, ISSESS, IDESS, LIDSP, BLKCOL,
     &   IDELB, NAMELB, MAPEL, MAPND)
C=======================================================================

C   --*** LNMAIN *** (PATHLN) PATHLINE main plot routine
C   --   Modified by John Glick - 11/29/88
C   --   Written by Amy Gilkey - revised 05/31/88
C   --   Modified version 1.1a  - November 1990 -  R.J. Meyers
C   --           added color coded sphere capability
C   --
C   --LNMAIN first determines the new face array, if it has changed.
C   --It then reads the pathline data for all selected times steps.
C   --The undeformed mesh is plotted with the pathline data.
C   --
C   --For 3D, elements are divided into faces and sorted as follows:
C   --   -1) Surface face
C   --    0) Interior element block boundary
C   --    n) Interior within element block 'n'
C   --
C   --Mesh lines are sorted as follows:
C   --   -1) Mesh boundary
C   --    0) Element block boundary
C   --    n) Interior within element block 'n'
C   --
C   --Dynamic memory arrays:
C   --   IN2ELB - the element block for each node;
C   --      <0 if not in any selected element block
C   --      =0 if in more than one selected element block
C   --   IDN2B - the element block for each dead node; dead if >= 0
C   --   IF2EL - the element number of each face
C   --   XN, YN, ZN - the nodal coordinates (ZN for 3D only)
C   --   HZ, VT, PD - the undeformed nodal coordinates
C   --      (rotated for 3D, PD for 3D only)
C   --   XF, YF, ZF - the undeformed face center coordinates
C   --      (rotated for 3D, ZF for 3D only)
C   --   HIDENP(i) - true iff node i is hidden (3D only)
C   --   HIDEF(i) - true iff face i is hidden (3D only)
C   --   XLN, YLN, ZLN - the pathline data
C   --
C   --Parameters:
C   --   A - IN - the dynamic memory base array
C   --   NAMECO - IN - the coordinates names
C   --   NAMES - IN - the variable names
C   --   NPTIMS - IN - the number of selected time steps
C   --   IPTIMS - IN - the selected time steps
C   --   TIMES - IN - the database times
C   --   WHOTIM - IN - true iff whole (versus history) time step
C   --   LENF - IN/OUT - the cumulative face counts by element block
C   --      LENF(0) is always 0
C   --      LENF(1..NELBLK) is the end of the surface faces of element block (i)
C   --      LENF(NELBLK+1) is the end of the interior faces
C   --      LENF(NELBLK+2) is the end of the faces that are dead
C   --      LENF(NELBLK+3) is the end of the faces outside a cut
C   --      LENF(NELBLK+4) is the end of the faces in a non-selected element
C   --         block
C   --   NLNKE - IN - the number of nodes per element
C   --   NLNKF - IN - the number of nodes per face
C   --   KLINKF - IN/OUT - the dynamic memory index of the connectivity
C   --      for all faces
C   --   LENL - IN/OUT - the cumulative line counts by element block
C   --   KLNSET - IN/OUT - the dynamic memory index of the sorted line set
C   --   NEWELB - IN/OUT - the new element blocks flag:
C   --      0 = no new element blocks
C   --      1 = new selected element blocks
C   --      2 = new displayed element blocks (implies new selected blocks)
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
C   --   Uses NDIM, NUMEL, NELBLK of /DBNUMS/
C   --   Uses IS3DIM, NNPSUR, NUMNPF, LLNSET of /D3NUMS/
C   --   Uses MSHDEF, MSHNUM, MSHLIN, MLNTYP, IHIDOP of /MSHOPT/
C   --   Uses MODDET, MODTYP, IDTVAR, NNDVAR, NEDVAR of /DETOPT/
C   --   Sets and uses ZMMESH of /MSHLIM/

      PARAMETER (MSHNON=0, MSHBOR=1, MSHDIV=2, MSHSEL=3, MSHALL=4)

      PARAMETER (KLFT=1, KRGT=2, KBOT=3, KTOP=4, KNEA=5, KFAR=6)

      include 'params.blk'
      include 'dbnums.blk'
      include 'dbnumgq.blk'
      include 'd3nums.blk'
      include 'mshopt.blk'
      include 'mshlim.blk'
      include 'lnvars.blk'

      CHARACTER*(MXSTLN) CDUM
      CHARACTER*(MXSTLN) NAMELB(*)

      DIMENSION A(*)
      CHARACTER*(*) NAMECO(*)
      CHARACTER*(*) NAMES(*)
      INTEGER IPTIMS(*)
      REAL TIMES(*)
      LOGICAL WHOTIM(*)
      INTEGER LENF(0:NELBLK+4)
      INTEGER NLNKE(NELBLK), NLNKF(NELBLK)
      INTEGER LENL(-2:NELBLK)
      INTEGER NEWELB
      INTEGER IELBST(NELBLK)
      INTEGER ISSNPS(NUMNPS,4)
      INTEGER ISSESS(NUMESS,4)
      INTEGER IDNPS(*)
      INTEGER IDESS(*)
      INTEGER LIDSP(0:*)
      INTEGER BLKCOL(0:NELBLK)
      INTEGER IDELB(*)
      INTEGER MAPEL(*), MAPND(*)

      INTEGER NDEFVW, IXVW
      LOGICAL NEWSET, NEWFAC, FIXFAC
      LOGICAL ANYDEF, ANYUND
      LOGICAL DOIXF, DON2B, DOELED, DOELEU, DODEAD, DONPS, DOESS, DOSCAL

************************************************************************
***                      Prepare for new plot set                    ***
************************************************************************

C   --Initialize flags

      CALL MSFLAG (ANYDEF, ANYUND,
     &   DOIXF, DON2B, DOELED, DOELEU, DODEAD, DONPS, DOESS, DOSCAL,
     &   MINMSH, MAXMSH, MAXHID)

      NEWSET = .TRUE.

C   --Reserve memory for plot set

      CALL MSMEMY (A, ANYDEF, ANYUND,
     &   DOIXF, DON2B, DOELED, DOELEU, DODEAD, DONPS, DOESS,
     &   LENF(NELBLK), KN2ELB, KDN2B, KIF2EL,
     &   KXN, KYN, KZN, KHZ, KVT, KPD, KDHZ, KDVT, KDPD,
     &   KXF, KYF, KZF, KDXF, KDYF, KDZF, KHIDEN, KHIDEF, KIXFAC,
     &   KSNPS, KSESS)
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 120

C   --Memory for variables is reserved at this time to prevent problems
C   --if the plot set is aborted; the use of memory is controlled with
C   --MDLONG
      CALL MDRSRV ('XLN', KXLN, 0)
      CALL MDRSRV ('YLN', KYLN, 0)
      IF (IS3DIM) CALL MDRSRV ('ZLN', KZLN, 0)
      IF (IS3DIM) CALL MDRSRV ('IXSCR', KIXSCR, 0)

      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 130

************************************************************************
***                         Plot sequence loop                       ***
************************************************************************

C   --Set the number of non-empty views for single time plots
      NVWPLT = NDEFVW (.FALSE.)

      ISTEP = 0

C   --Wipe out the mesh plot pick
      CALL INPICK ('NONE')

C   --Set up the faces for the new plot set

      CALL MSGEOM (A, 'PATHLINE', ISTEP,
     &   LENF, NLNKF, KLINKF, KXN, KYN, KZN,
     &   KIF2EL, NEWELB, IELBST, NEWFAC, FIXFAC)
      IF (LENF(NELBLK) .LE. 0) GOTO 120
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 120

      CALL MSSTEP (A, ISTEP, NEWSET, NEWFAC, ANYDEF, ANYUND,
     &   DOIXF, DON2B, DOELED, DOELEU, DODEAD, DONPS, DOESS,
     &   LENF, NLNKF, A(KLINKF), LENL, KLNSET,
     &   KHIDEN, KHIDEF, KIXFAC,
     &   KXN, KYN, KZN, KHZ, KVT, KPD, KDHZ, KDVT, KDPD,
     &   KXF, KYF, KZF, KDXF, KDYF, KDZF,
     &   NEWELB, IELBST, KN2ELB, KDN2B, A(KIF2EL), KNPSUR)
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 120

C   --Calculate the axis limits and set up graphics, if limits changed

      IF (DOSCAL) THEN

C      --Compute NPSURF nodes that determine the mesh limits

         IF ((.NOT. IS3DIM) .AND. (NNPSUR .LT. 0)) THEN
            CALL MDLONG ('NPSURF', KNPSUR, 0)
            CALL MDLONG ('NPSURF', KNPSUR, NUMNPF)
            CALL MDSTAT (NERR, MEM)
            IF (NERR .GT. 0) GOTO 120

            CALL MAKSU2 (LENL, A(KLNSET), MSHBOR,
     &         DODEAD, A(KDN2B), A(KNPSUR))

            CALL MDLONG ('NPSURF', KNPSUR, NNPSUR)
         END IF
      END IF

      CALL MSSCAL (DOSCAL, NNPSUR, A(KNPSUR),
     &   A(KHZ), A(KVT), A(KPD))

C   --Read variables for plot

      CALL MDLONG ('XLN', KXLN, NPTIMS * NLNCRV)
      CALL MDLONG ('YLN', KYLN, NPTIMS * NLNCRV)
      IF (IS3DIM) CALL MDLONG ('ZLN', KZLN, NPTIMS * NLNCRV)
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 130

      CALL LNREAD (A, NPTIMS, NPTIMW, IPTIMS, TIMES, WHOTIM,
     &   A(KXLN), A(KYLN), A(KZLN))

C   --Label plot frame

  100 CONTINUE
      CALL LNLAB (A, NEWSET, NPTIMS, IPTIMS, TIMES,
     &   NAMECO, NAMES, IELBST,
     &   ISSNPS, IDNPS, ISSESS, IDESS, A(KSNPS), A(KSESS),
     &   LIDSP, BLKCOL, *120)

      DO 110 IVW = 1, NVWPLT
         IVIEW = IXVW (.FALSE., IVW)

         KTXN = KHZ
         KTYN = KVT
         KTZN = KPD
         KTXF = KXF
         KTYF = KYF
         KTZF = KZF

C      --Identify the hidden faces and nodes on the surface

         IF (IS3DIM .AND. (NEWSET .OR. NEWFAC)) THEN
C         --NOTE: coordinates may be mangled by HIDDEN
            CALL HIDDEN (A, MAXHID, LENF, NLNKE, NLNKF, A(KLINKF),
     &         A(KTXN), A(KTYN), A(KTZN),
     &         LENL, A(KLNSET), MSCTYP, ZMMESH,
     &         MINMSH, MAXMSH, IELBST,
     &         DODEAD, A(KDN2B), A(KHIDEF), A(KHIDEN),
     &         A(KTZF), DOIXF, NXFAC, A(KIXFAC), NAMELB)

            CALL MDSTAT (NERR, MEM)
            IF (NERR .GT. 0) GOTO 120
         END IF

C      --Set up to plot view

         CALL SETVW (IVIEW, *120)

C      --Plot undeformed mesh

C      --Set up the mesh plot pick
         CALL INPICK ('UNDEFORM')
         CDUM = ' '
         CALL MSPLT1 (A, .FALSE.,
     &      MSHNUM(IVIEW), MSHLIN(IVIEW), MLNTYP(-1,IVIEW),
     &      LENF, NLNKF, A(KLINKF), LENL, A(KLNSET),
     &      A(KHIDEN), A(KHIDEF),
     &      A(KTXN), A(KTYN), A(KTZN), A(KTXF), A(KTYF), A(KTZF),
     &      IELBST, A(KN2ELB), DODEAD, A(KDN2B),
     &      .FALSE., 0, IDUM,
     &      NNPSET(IVIEW), ISSNPS(1,IVIEW),
     &      NESSET(IVIEW), ISSESS(1,IVIEW), BLKCOL,
     &      IDELB, VARNP, CDUM, IHIDOP, MAPEL, MAPND, *120)

         IF (IS3DIM) CALL MDLONG ('IXSCR', KIXSCR, NPTIMS)
         CALL MDSTAT (NERR, MEM)
         IF (NERR .GT. 0) GOTO 120

         CALL LNPLOT (NPTIMS, NPTIMW,
     &      A(KXLN), A(KYLN), A(KZLN), A(KIXSCR), *120)

         IF (IS3DIM) CALL MDLONG ('IXSCR', KIXSCR, 0)

         IF (NEWSET) THEN
            NEWSET = .FALSE.
            NEWFAC = .FALSE.
            FIXFAC = .FALSE.
            DOSCAL = .FALSE.
         END IF

  110 CONTINUE

C   --Check if user wants to quit or get hardcopy of plot

C   --Set color in case text is requested
      CALL UGRCOL (0, BLKCOL)
      N = 1
      CALL GRPEND (.TRUE., .TRUE., N, N, .FALSE., *100, *120)

  120 CONTINUE

C   --Release the memory reserved for plot set

      CALL MSDONE (ANYDEF, ANYUND,
     &   DOIXF, DON2B, DOELED, DOELEU, DODEAD, DONPS, DOESS)

      CALL MDDEL ('XLN')
      CALL MDDEL ('YLN')
      IF (IS3DIM) CALL MDDEL ('ZLN')
      IF (IS3DIM) CALL MDDEL ('IXSCR')

  130 CONTINUE
      RETURN
      END
