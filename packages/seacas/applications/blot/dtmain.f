C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE DTMAIN (A, NAMECO, NAMES, NPTIMS, IPTIMS, TIMES,
     &   LENF, NLNKE, NLNKF, KLINKF, LENL, KLNSET,
     &   NEWELB, IELBST, KNPSUR, ISEVOK,
     &   ISSNPS, IDNPS, ISSESS, IDESS, LIDSP, BLKCOL, IDELB, NAMELB,
     *  MAPEL, MAPND)
C=======================================================================

C   --*** DTMAIN *** (DETOUR) DETOUR main plot routine
C   --   Modified by John Glick - 11/29/88
C   --   Written by Amy Gilkey - revised 05/31/88
C   --   D. P. Flanagan, 11/17/82
C   --
C   --DTMAIN first determines the new face array, if it has changed.
C   --It then loops for each selected time step.  For each step, the
C   --deformed coordinates are calculated, the needed database variables
C   --are input, and the plots are displayed.
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
C   --DTMAIN calls DTREAD to pre-process variables for modes involving
C   --variables.  This includes obtaining the variables and converting
C   --from element to nodal variables or vice versa.
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
C   --   DHZ, DVT, DPD - the deformed nodal coordinates
C   --      (rotated for 3D, DPD for 3D only)
C   --   XF, YF, ZF - the undeformed face center coordinates
C   --      (rotated for 3D, ZF for 3D only)
C   --   DXF, DYF, DZF - the deformed face center coordinates
C   --      (rotated for 3D, DZF for 3D only)
C   --   HIDENP(i) - true iff node i is hidden (3D only)
C   --   HIDEF(i) - true iff face i is hidden (3D only)
C   --   VARNP - the nodal variable values
C   --   VARFAC - the face variable values
C   --   VAREL - sized to hold an element variable
C   --
C   --Parameters:
C   --   A - IN - the dynamic memory base array
C   --   NAMECO - IN - the coordinates names
C   --   NAMES - IN - the variable names
C   --   NPTIMS - IN - the number of selected time steps
C   --   IPTIMS - IN - the selected time steps
C   --   TIMES - IN - the database times
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
C   --   ISEVOK - IN - the element block variable truth table;
C   --      variable i of block j exists iff ISEVOK(j,i)
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
C   --   Uses NDIM, NUMEL, NELBLK, NVARNP, NVAREL of /DBNUMS/
C   --   Uses IS3DIM, NNPSUR, NUMNPF, LLNSET of /D3NUMS/
C   --   Uses DEFOK, DFAC of /DEFORM/
C   --   Uses MSHDEF, MSHNUM, MSHLIN, MLNTYP, IHIDOP, NALVAR, DEADNP of /MSHOPT/
C   --   Uses MODDET, MODTYP, IDTVAR, NNDVAR, NEDVAR of /DETOPT/
C   --   Uses MULTIM of /VIEWS/
C   --   Sets and uses ZMMESH of /MSHLIM/

      PARAMETER (MSHNON=0, MSHBOR=1, MSHDIV=2, MSHSEL=3, MSHALL=4)

      PARAMETER (KLFT=1, KRGT=2, KBOT=3, KTOP=4, KNEA=5, KFAR=6)

      EXTERNAL BLKDAT

      include 'params.blk'
      include 'dbnums.blk'
      include 'dbnumgq.blk'
      include 'd3nums.blk'
      include 'deform.blk'
      include 'mshopt.blk'
      include 'detopt.blk'
      include 'views.blk'
      include 'mshlim.blk'
      include 'nodzom.blk'
      include 'axsplt.blk'
      include 'cntr.blk'
      include 'dbase.blk'
      include 'icexct.blk'

      LOGICAL ISLINE, GOBCK

      DIMENSION A(*)
      CHARACTER*(MXSTLN) NAMECO(*)
      CHARACTER*(*) NAMES(*)
      INTEGER IPTIMS(*)
      REAL TIMES(*)
      INTEGER LENF(0:NELBLK+4)
      INTEGER NLNKE(NELBLK), NLNKF(NELBLK)
      INTEGER LENL(-2:NELBLK)
      INTEGER NEWELB
      INTEGER IELBST(NELBLK)
      LOGICAL ISEVOK(NELBLK,*)
      INTEGER ISSNPS(NUMNPS,4)
      INTEGER ISSESS(NUMESS,4)
      INTEGER IDNPS(*)
      INTEGER IDESS(*)
      INTEGER LIDSP(0:*)
      INTEGER BLKCOL(0:NELBLK)
      INTEGER IDELB(*)
      CHARACTER*(MXSTLN) NAMELB(*)
      INTEGER MAPEL(*), MAPND(*)

      INTEGER NUMMOD, NDEFVW, IXVW
      CHARACTER TYP
      LOGICAL NEWSET, NEWFAC, FIXFAC
      LOGICAL ANYDEF, ANYUND
      LOGICAL DOIXF, DON2B, DOELED, DOELEU, DODEAD, DONPS, DOESS, DOSCAL
      LOGICAL DOVN2B
      LOGICAL ONESHO
      LOGICAL DEFORM, LSTDEF

************************************************************************
***                      Prepare for new plot set                    ***
************************************************************************

C  Save contour values
      CMINSV=CMIN
      CMAXSV=CMAX
      DELCSV=DELC
      NCSAV=NCNTR

C   --Initialize flags

      CALL MSFLAG (ANYDEF, ANYUND,
     &   DOIXF, DON2B, DOELED, DOELEU, DODEAD, DONPS, DOESS, DOSCAL,
     &   MINMSH, MAXMSH, MAXHID)

      IF (IS3DIM .AND. (IHIDOP .GE. 4)) THEN
         IF ((NUMMOD (MODDET, MODTYP, 'CONTOUR', 'PAINT') .GE. 1)
     &      .OR. (NUMMOD (MODDET, MODTYP, 'ELEMCONT', 'PAINT') .GE. 1)
     &      .OR. (NUMMOD (MODDET, ' ', 'SOLID', ' ') .GE. 1))
     &      DOIXF = .TRUE.
      END IF

      IF (NNDVAR .GT. 0) DON2B = .TRUE.

      IF ((NEDVAR .GT. 0) .OR. DOIXF
     &   .OR. (NUMMOD (MODDET, ' ', 'CONTOUR', ' ') .GE. 1)
     &   .OR. (NUMMOD (MODDET, ' ', 'ELEMCONT', ' ') .GE. 1)) THEN
         IF (ANYDEF) DOELED = .TRUE.
         IF (ANYUND) DOELEU = .TRUE.
      END IF

      NVEC = 0
      IF (NUMMOD (MODDET, ' ', 'VECTOR', ' ') .GE. 1) THEN
         IF ((NUMMOD (MODDET, MODTYP, 'VECTOR', 'NODE') .GE. 1)
     &      .OR. (NUMMOD (MODDET, MODTYP, 'VECTOR', 'ELEMENT') .GE. 1))
     &      NVEC = MAX (NVEC, NDIM)
         IF ((NUMMOD (MODDET, MODTYP, 'VECTOR', 'SIGMAX') .GE. 1)
     &      .OR. (NUMMOD (MODDET, MODTYP, 'VECTOR', 'SIGMIN') .GE. 1))
     &      NVEC = MAX (NVEC, 3)
      END IF

      ONESHO = (.NOT. ANYDEF) .AND. ((NNDVAR + NEDVAR) .LE. 0)
     &   .AND. (NALVAR .LE. 0)

      NEWSET = .TRUE.

C   --Reserve memory for plot set

      CALL MSMEMY (A, ANYDEF, ANYUND,
     &   DOIXF, DON2B, DOELED, DOELEU, DODEAD, DONPS, DOESS,
     &   LENF(NELBLK), KN2ELB, KDN2B, KIF2EL,
     &   KXN, KYN, KZN, KHZ, KVT, KPD, KDHZ, KDVT, KDPD,
     &   KXF, KYF, KZF, KDXF, KDYF, KDZF, KHIDEN, KHIDEF, KIXFAC,
     &   KSNPS, KSESS)
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 160

C   --Memory for variables is reserved at this time to prevent problems
C   --if the plot set is aborted; the use of memory is controlled with
C   --MDLONG
      IF ((NNDVAR + NEDVAR) .GT. 0) THEN
         CALL MDRSRV ('VARNP', KVARNP, 0)
C      --Note that VARFAC must be resized after element birth/death
         CALL MDRSRV ('VARFAC', KVARF, 0)
         CALL MDRSRV ('ISVOK', KISVOK, NELBLK)
         DOVN2B = .TRUE.
         CALL MDRSRV ('IVN2B', KIVN2B, 0)
         CALL MDRSRV ('IFVCNT', KFVCNT, 0)
      END IF

      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 170

************************************************************************
***                         Plot sequence loop                       ***
************************************************************************

C   --Set the number of non-empty views for multiple and single time plots
      IF (MULTIM) THEN
         NVWTIM = NDEFVW (.FALSE.)
         NVWPLT = 1
      ELSE
         NVWPLT = NDEFVW (.FALSE.)
         NVWTIM = 1
      END IF

      IF (ONESHO) THEN
         NPTIMX = 1
      ELSE
         NPTIMX = NPTIMS
      END IF

      NPLOT1 = 1
 150  CONTINUE

C      --Set number of times for last plot
         NVWTIM = MIN (NVWTIM, NPTIMX-NPLOT1+1)

  100    CONTINUE
         DO 140 IVWTIM = 1, NVWTIM

            NPLOT = NPLOT1 + IVWTIM-1
            ISTEP = IPTIMS(NPLOT)

C         --Wipe out the mesh plot pick
            CALL INPICK ('NONE')

C         --Set up the faces for the new plot set

            IF (NEWSET .OR. (NALVAR .GT. 0)) THEN
               CALL MSGEOM (A, 'DETOUR', ISTEP,
     &            LENF, NLNKF, KLINKF, KXN, KYN, KZN,
     &            KIF2EL, NEWELB, IELBST, NEWFAC, FIXFAC)
               IF (LENF(NELBLK) .LE. 0) GOTO 160
               CALL MDSTAT (NERR, MEM)
               IF (NERR .GT. 0) GOTO 160

               IF (FIXFAC) DOVN2B = .TRUE.
            END IF

            CALL MSSTEP (A, ISTEP, NEWSET, NEWFAC, ANYDEF, ANYUND,
     &         DOIXF, DON2B, DOELED, DOELEU, DODEAD, DONPS, DOESS,
     &         LENF, NLNKF, A(KLINKF), LENL, KLNSET,
     &         KHIDEN, KHIDEF, KIXFAC,
     &         KXN, KYN, KZN, KHZ, KVT, KPD, KDHZ, KDVT, KDPD,
     &         KXF, KYF, KZF, KDXF, KDYF, KDZF,
     &         NEWELB, IELBST, KN2ELB, KDN2B, A(KIF2EL), KNPSUR)
            CALL MDSTAT (NERR, MEM)
            IF (NERR .GT. 0) GOTO 160

C         --Calculate the axis limits and set up graphics, if limits changed

            IF (NEWSET .OR. ((IVWTIM .EQ. 1) .AND. DOSCAL)) THEN
               IF (DOSCAL) THEN

C               --Compute NPSURF nodes that determine the mesh limits

                  IF ((.NOT. IS3DIM) .AND. (NNPSUR .LT. 0)) THEN
                     CALL MDLONG ('NPSURF', KNPSUR, 0)
                     CALL MDLONG ('NPSURF', KNPSUR, NUMNPF)
                     CALL MDSTAT (NERR, MEM)
                     IF (NERR .GT. 0) GOTO 160

                     CALL MAKSU2 (LENL, A(KLNSET), MSHBOR,
     &                  DODEAD, A(KDN2B), A(KNPSUR))

                     CALL MDLONG ('NPSURF', KNPSUR, NNPSUR)
                  END IF
               END IF

               IF (ANYDEF) THEN
                  CALL MSSCAL (DOSCAL, NNPSUR, A(KNPSUR),
     &               A(KDHZ), A(KDVT), A(KDPD))
               ELSE
                  CALL MSSCAL (DOSCAL, NNPSUR, A(KNPSUR),
     &               A(KHZ), A(KVT), A(KPD))
               END IF

C            --Set up scaling for the vectors, if needed

               IF (NVEC .GT. 0) THEN
                  IF (NEWSET) THEN
                     VVMAX = -1.0E-36
                     DO 110 I = 1, NVEC
                        IF (IDTVAR(I) .GT. 0) THEN
C                        --Note that SCALER does not print result, and that
C                        --min/max has been calculated in DTCOMD
                          CALL SCALER (A, A,
     &                      0, NAMES(IDTVAR(I)), IDTVAR(I), .TRUE.,
     *                      IELBST, NALVAR, FMIN, FMAX, MAPEL, MAPND)
                          VVMAX = MAX (VVMAX, ABS(FMIN), ABS(FMAX))
                        END IF
  110                CONTINUE
                     IF (VVMAX .EQ. 0.0) VVMAX = 1.0
                  END IF
                  IF (NEWSET .OR. DOSCAL) THEN
                     VECMAX = MAX (ZMMESH(KRGT)-ZMMESH(KLFT)
     &                  ,          ZMMESH(KTOP)-ZMMESH(KBOT)) / VVMAX
                  END IF
               END IF
            END IF

C         --Read variables for plot; done here so min/max can be calculated
C         --for the label

            IF ((NNDVAR + NEDVAR) .GT. 0) THEN
               CALL MDLONG ('VARNP', KVARNP, NNDVAR * NUMNPF)
               CALL MDLONG ('VARFAC', KVARF, NEDVAR * LENF(NELBLK))
               CALL MDRSRV ('VAREL', KVAREL, MIN (NEDVAR, 1) * NUMEL)
               CALL MDSTAT (NERR, MEM)
               IF (NERR .GT. 0) GOTO 170

               LVARF = LENF(NELBLK)

               CALL DTREAD (A, ISTEP, IDTVAR, NNDVAR, NEDVAR,
     &            LENF, A(KIF2EL), IELBST, ISEVOK,
     &            A(KVARNP), A(KVARF), A(KVAREL), LVARF)

               CALL MDDEL ('VAREL')

C            --If contour mode, convert element variable to nodal variable,
C            --if needed, and count min/max

               IF (NUMMOD (MODDET, ' ', 'CONTOUR', ' ') .GE. 1) THEN
                  CALL DBVTYP_BL (IDTVAR(1), TYP, ID)
                  IF (TYP .EQ. 'E') THEN
                     IF (DOVN2B) THEN
                        CALL MDLONG ('IVN2B', KIVN2B, 0)
                        CALL MDLONG ('IVN2B', KIVN2B, NUMNPF)
                        CALL MDLONG ('IFVCNT', KFVCNT, 0)
                        CALL MDLONG ('IFVCNT', KFVCNT, NUMNPF)
                        CALL MDSTAT (NERR, MEM)
                        IF (NERR .GT. 0) GOTO 170
                     END IF
                     CALL CFV2NV (LENF, NLNKF, A(KLINKF),
     &                  IELBST, ISEVOK(1,ID), A(KVARF),
     &                  DOVN2B, A(KIVN2B), A(KFVCNT),
     &                  A(KVARNP))
                  ELSE
                     KIVN2B = KN2ELB
                  END IF
               END IF
            END IF

C         --If nodal or element contour mode, count min/max

            IF (NUMMOD (MODDET, ' ', 'CONTOUR', ' ') .GE. 1) THEN
               CALL CNVMAX (NUMNPF, A(KVARNP), A(KIVN2B),
     &            NMIN, NMAX, FMIN, FMAX)
C ... Echo contour minimum and maximum values to log file.
               if (nlog .gt. 0) then
                  WRITE (nlog, 9000) 'Nodal', NAMES(IDTVAR(1)),
     $                 FMIN, FMAX
               end if
             ELSE IF (NUMMOD (MODDET, ' ', 'ELEMCONT', ' ') .GE. 1)
     &           THEN
               CALL CFVMAX (LENF(NELBLK), A(KVARF),
     &           NMIN, NMAX, FMIN, FMAX)
C ... Echo contour minimum and maximum values to log file.
               if (nlog .gt. 0) then
                  WRITE (nlog, 9000) 'Element', NAMES(IDTVAR(1)),
     $                 FMIN, FMAX
               end if
             ELSE
               NMIN = 0
               NMAX = 0
               FMIN = 0.0
               FMAX = 0.0
            END IF
 9000       FORMAT ('$$$ ',A,'= ',A,
     *        ' FMIN=',1pe13.5,' FMAX=',1pe13.5)

            IF(FMIN.NE.FMAX.AND.NCNTR.NE.0.AND.IEXCON.EQ.1)THEN
C  RESET CONTOUR VALUES TO THE MIN AND MAX FOR THIS PLOT
             IF (NUMMOD (MODTYP, ' ', 'LINE', ' ') .GE. 1) THEN
              ISLINE=.TRUE.
              CALL CONRNG(ISLINE,FMIN,FMAX,NCNTR,DELC,CMIN,CMAX)
             ELSE
              ESPRNG=.0001*(FMAX-FMIN)
              CMIN=FMIN-ESPRNG
              CMAX=FMAX+ESPRNG
              DELC=(CMAX-CMIN)/NCNTR
             ENDIF
            ENDIF

C         --Label plot frame

  120       CONTINUE
            IF (IVWTIM .EQ. 1) THEN
               N = NVWTIM
               IF (ONESHO) N = 0
               CALL DTLAB (A, NEWSET, N, IPTIMS(NPLOT), TIMES,
     &            NAMECO, NAMES, IELBST,
     &            NMIN, NMAX, FMIN, FMAX,
     &            ISSNPS, IDNPS, ISSESS, IDESS, A(KSNPS), A(KSESS),
     &            LIDSP, BLKCOL, *160)
            END IF

C -- IF PLOTTING AXIS ONLY, SKIP TO END

            IF(AXONLY) THEN
                AXONLY = .FALSE.
                GO TO 180
            END IF

            DO 130 IVW = 1, NVWPLT
               IVIEW = IXVW (.FALSE., MAX (IVWTIM, IVW))

               DEFORM = (ANYDEF .AND. (MSHDEF(IVIEW) .NE. 'UNDEFORM'))

               IF (DEFORM) THEN
                  KTXN = KDHZ
                  KTYN = KDVT
                  KTZN = KDPD
                  KTXF = KDXF
                  KTYF = KDYF
                  KTZF = KDZF
               ELSE
                  KTXN = KHZ
                  KTYN = KVT
                  KTZN = KPD
                  KTXF = KXF
                  KTYF = KYF
                  KTZF = KZF
               END IF

C FOR "ZOOM NODE" MODE, RECALCULATE THE ZOOM WINDOW

               IF(NZMON) THEN
                  CALL ZOOMND(A(KTXN), A(KTYN), A(KTZN), RDMESH)
                  CALL EXPLIM(2, RDMESH, RDMESH)
                  CALL ADJLIM(MSHDEF, XISSYM, YISSYM, LFTSYM, BOTSYM,
     &                        XAXSYM, YAXSYM, SQMESH, RDMESH, ZMMESH)
                  CALL SETUP(MSHDEF,ZMMESH)
               END IF

C            --Identify the hidden faces and nodes on the surface

               IF (IS3DIM .AND. (NEWSET .OR. NEWFAC
     &            .OR. DEFORM .OR. LSTDEF)) THEN
C               --NOTE: coordinates may be mangled by HIDDEN
                  CALL HIDDEN (A, MAXHID, LENF, NLNKE, NLNKF, A(KLINKF),
     &               A(KTXN), A(KTYN), A(KTZN),
     &               LENL, A(KLNSET), MSCTYP, ZMMESH,
     &               MINMSH, MAXMSH, IELBST,
     &               DODEAD, A(KDN2B), A(KHIDEF), A(KHIDEN),
     &               A(KTZF), DOIXF, NXFAC, A(KIXFAC), NAMELB)
                  LSTDEF = DEFORM

                  CALL MDSTAT (NERR, MEM)
                  IF (NERR .GT. 0) GOTO 160
               END IF

C            --Set up to plot view

               CALL SETVW (IVIEW, *160)

C            --Plot view, with deformed or undeformed mesh

C            --Set up the mesh plot pick
               IF (DEFORM) THEN
                  CALL INPICK ('DEFORM')
               ELSE
                  CALL INPICK ('UNDEFORM')
               END IF

               CALL DTPLT1 (A,
     &            MSHNUM(IVIEW), MSHLIN(IVIEW), MLNTYP(-1,IVIEW),
     &            MODDET(IVIEW), MODTYP(IVIEW),
     &            LENF, NLNKF, A(KLINKF), LENL, A(KLNSET),
     &            A(KHIDEN), A(KHIDEF), DOIXF, NXFAC, A(KIXFAC),
     &            A(KTXN), A(KTYN), A(KTZN), A(KTXF), A(KTYF), A(KTZF),
     &            IELBST, ISEVOK,
     &            A(KN2ELB), A(KIVN2B), DODEAD, A(KDN2B),
     &            NNPSET(IVIEW), ISSNPS(1,IVIEW),
     &            NESSET(IVIEW), ISSESS(1,IVIEW),
     &            IDTVAR, A(KVARNP), A(KVARF), NMIN, NMAX, FMIN, FMAX,
     &            VECMAX, A(KISVOK), BLKCOL, IDELB, MAXHID,
     *            MAPEL, MAPND, *160)

               IF (NEWSET) THEN
                  NEWSET = .FALSE.
                  NEWFAC = .FALSE.
                  FIXFAC = .FALSE.
                  IF (.NOT. ANYDEF) DOSCAL = .FALSE.
               END IF
  130       CONTINUE

C              Outline the viewport windows

            CALL OUTLIN (BLKCOL, *160)

C         --Check if user wants to quit or get hardcopy of single-time plot

            IF (.NOT. MULTIM) THEN
C            --Set color in case text is requested
               CALL UGRCOL (0, BLKCOL)
               GOBCK = .TRUE.
               CALL GRPEND (.TRUE., .TRUE., NPLOT1, NPTIMX, GOBCK,
     $              *120, *160)
            END IF

C         --Release variable memory
            IF ((NNDVAR + NEDVAR) .GT. 0) THEN
               CALL MDLONG ('VARNP', KVARNP, 0)
               CALL MDLONG ('VARFAC', KVARF, 0)
            END IF
            CALL MDSTAT (NERR, MEM)
            IF (NERR .GT. 0) GOTO 170

  140    CONTINUE

C      --Check if user wants to quit or get hardcopy for multiple time plot

         IF (MULTIM) THEN
C         --Set color in case text is requested
            CALL UGRCOL (0, BLKCOL)
            I = INT ((NPTIMX-1) / DBLE(NVWTIM)) + 1
            GOBCK = .FALSE.
            CALL GRPEND (.TRUE., .TRUE., NPLOT1, I, GOBCK, *100, *160)
         END IF
         if (gobck) then
            nplot1 = nplot1 - nvwtim
            if (nplot1 .lt. 1) nplot1 = 1
         else
            NPLOT1 = NPLOT1 + NVWTIM
         end if
         if (nplot1 .le. nptimx) go to 150

  160 CONTINUE

C   --Release the memory reserved for plot set

180   CONTINUE

      CALL MSDONE (ANYDEF, ANYUND,
     &   DOIXF, DON2B, DOELED, DOELEU, DODEAD, DONPS, DOESS)

      IF ((NNDVAR + NEDVAR) .GT. 0) THEN
         CALL MDDEL ('VARNP')
         CALL MDDEL ('VARFAC')
         CALL MDDEL ('ISVOK')
         DOVN2B = .FALSE.
         CALL MDDEL ('IVN2B')
         CALL MDDEL ('IFVCNT')
      END IF

  170 CONTINUE

C  Restore contour values
      CMIN=CMINSV
      CMAX=CMAXSV
      DELC=DELCSV
      NCNTR=NCSAV

      RETURN
      END
