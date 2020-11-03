C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE MSPLT1 (A, WIDLIN, MSHNUM, MSHLIN, MLNTYP,
     &   LENF, NLNKF, LINKF, LENL, LINSET,
     &   HIDENP, HIDEF, XN, YN, ZN, XF, YF, ZF,
     &   IELBST, IN2ELB, DODEAD, IDN2B,
     &   NODSEL, NNESEL, NESEL,
     &   NNPSET, ISSNPS, NESSET, ISSESS, BLKCOL,
     &   IDELB, VARNP, MODDET, IHIDOP, MAPEL, MAPND, *)
C=======================================================================

C   --*** MSPLT1 *** (MESH) Plot one view
C   --   Modified by John Glick - 11/29/88
C   --   Written by Amy Gilkey - revised 04/11/88
C   --   Modified version 1.1a  November 1990  - R.J. Meyers
C   --           added color coded sphere capability
C   --
C   --MSPLT1 displays the mesh for a single view.  The labeling of the
C   --view is done elsewhere.
C   --
C   --SSMEMY is called to get the node set and side set
C   --information.
C   --
C   --This routine uses MDFIND to find the following dynamic memory arrays:
C   --   IX2NP - the node number for each mesh index
C   --   IF2EL - the element number of each face
C   --   IE2ELB - the element block for each element
C   --
C   --Parameters:
C   --   A - IN - the dynamic memory array
C   --   WIDLIN - IN - true iff mesh lines should be wide versus narrow
C   --   MSHNUM - IN - the mesh numbering (as in /MSHOPT/)
C   --   MSHLIN - IN - the display type for the mesh lines (as in /MSHOPT/)
C   --   MLNTYP - IN - the line type of lines (as in /MSHOPT/)
C   --   LENF - IN - the cumulative face counts by element block
C   --   NLNKF - IN - the number of nodes per face
C   --   LINKF - IN - the connectivity for all faces
C   --   LENL - IN - the cumulative line counts by element block
C   --   LINSET - IN - the sorted line set
C   --   HIDENP(i) - IN - true iff node i is hidden (3D only)
C   --   HIDEF(i) - IN - true iff face i is hidden (3D only)
C   --   XN, YN, ZN - IN - the nodal coordinates
C   --   XF, YF, ZF - IN - the face center coordinates
C   --   IELBST - IN - the element block status (>0 if selected)
C   --   IN2ELB - IN - the element block for each node;
C   --      <0 if not in any selected element block
C   --      =0 if in more than one selected element block
C   --   DODEAD - IN - true iff dead nodes are to be displayed
C   --   IDN2B - IN - the element block for each dead node; dead if >= 0
C   --   NODSEL - IN - true iff nodes (versus elements) are selected
C   --   NNESEL - IN - the number of selected nodes/elements
C   --   NESEL - IN - the selected nodes/elements
C   --   NNPSET - IN - the number of selected node sets
C   --   ISSNPS - IN - the indices of the selected node sets
C   --   NESSET - IN - the number of selected side sets
C   --   ISSESS - IN - the indices of the selected side sets
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
C   --   VARNP - IN - the function value for spheres for coloring
C   --   MODDET - IN - DETOUR mode to check for painting
C   --   * - return statement if cancel function active
C   --
C   --Common Variables:
C   --   Uses NELBLK of /DBNUMS/
C   --   Uses IS3DIM, NUMNPF of /D3NUMS/
C   --   Uses SELOK of /SELNE/

      PARAMETER (MSHNON=0, MSHBOR=1, MSHDIV=2, MSHSEL=3, MSHALL=4)

      include 'dbnums.blk'
      include 'd3nums.blk'
      include 'sphele.blk'

      DIMENSION A(*)
      LOGICAL WIDLIN
      CHARACTER*(*) MSHNUM
      INTEGER MLNTYP(-1:1)
      INTEGER LENF(0:NELBLK)
      INTEGER NLNKF(NELBLK)
      INTEGER LINKF(*)
      INTEGER LENL(-2:NELBLK), LINSET(LLNSET,*)
      LOGICAL HIDENP(NUMNPF)
      LOGICAL HIDEF(*)
      REAL XN(NUMNPF), YN(NUMNPF), ZN(NUMNPF)
      REAL XF(*), YF(*), ZF(*)
      INTEGER IELBST(NELBLK)
      INTEGER IN2ELB(NUMNPF)
      LOGICAL DODEAD
      INTEGER IDN2B(NUMNPF)
      LOGICAL NODSEL
      INTEGER NESEL(*)
      INTEGER ISSNPS(*)
      INTEGER ISSESS(*)
      INTEGER BLKCOL(0:NELBLK)
      INTEGER IDELB(*)
      REAL VARNP(*)
      INTEGER MAPEL(*), MAPND(*)

      CHARACTER*8 MODDET

      CHARACTER*6 FFLAG

      CHARACTER*8 INUM

      LOGICAL FIRST
      DATA FIRST /.TRUE./

C   --Display mesh

      CALL MESHUP (WIDLIN, MSHLIN, MLNTYP, IELBST, LENL, LINSET,
     &   BLKCOL, XN, YN, ZN, IDELB, *100)

C   --Display "dead" nodes, if requested

      IF (DODEAD) THEN
         CALL DEADUP (HIDENP, XN, YN, ZN, IDN2B, *100)
      END IF

C   --Number nodes or elements, if requested

      IF ((MSHNUM .NE. 'NONE') .OR. (NNESEL .GT. 0)) THEN
         IF (NNESEL .GT. 0) THEN
            IF (NODSEL) THEN
               INUM = 'NODE'
            ELSE
               INUM = 'ELEMENT'
            END IF
         ELSE
            INUM = 'NONE'
         END IF
         CALL MDFIND ('IX2NP', KIX2NP, IDUM)
         CALL MDFIND ('IF2EL', KIF2EL, IDUM)
         CALL MDSTAT (NERR, MEM)
         IF (NERR .GT. 0) GOTO 100

         CALL NUMBER (MSHNUM, LENF, NLNKF, A(KIX2NP), A(KIF2EL),
     &     HIDENP, HIDEF, XN, YN, ZN, XF, YF, ZF,
     &     IELBST, IN2ELB, DODEAD, IDN2B,
     &     INUM, NNESEL, NESEL, BLKCOL, IDELB,
     *     MAPEL, MAPND, *100)
      END IF

C   --Mark node sets, if selected

      IF (NNPSET .GT. 0) THEN
         CALL SSMEMY (A,
     &      .TRUE., KIDNS, KNNNS, KIXNNS, KLTNNS,
     &      .FALSE., KIDSS, KNESS, KNNSS,
     &      KIXESS, KIXNSS, KLTESS, KLTNSS)
         CALL MDFIND ('IX2NP', KIX2NP, IDUM)
         CALL MDSTAT (NERR, MEM)
         IF (NERR .GT. 0) GOTO 100

         CALL MRKNPS (HIDENP, XN, YN, ZN,
     &      A(KIX2NP), IN2ELB, DODEAD, IDN2B,
     &      NNPSET, ISSNPS, A(KIDNS), A(KNNNS), A(KIXNNS), A(KLTNNS),
     &      *100)
      END IF

C   --Mark side sets, if selected

      IF (NESSET .GT. 0) THEN
         CALL SSMEMY (A,
     &      .FALSE., KIDNS, KNNNS, KIXNNS, KLTNNS,
     &      .TRUE., KIDSS, KNESS, KNNSS,
     &      KIXESS, KIXNSS, KLTESS, KLTNSS)
         CALL MDFIND ('IX2NP', KIX2NP, IDUM)
         CALL MDFIND ('IF2EL', KIF2EL, IDUM)
         CALL MDFIND ('IE2ELB', KE2ELB, IDUM)
         CALL MDSTAT (NERR, MEM)
         IF (NERR .GT. 0) GOTO 100

         CALL MRKESS (LENF, NLNKF, LINKF, HIDENP, HIDEF,
     &      XN, YN, ZN, XF, YF, ZF,
     &      IELBST, A(KIX2NP), IN2ELB, DODEAD, IDN2B,
     &      A(KIF2EL), A(KE2ELB),
     &      NESSET, ISSESS, A(KIDSS), A(KNESS), A(KNNSS),
     &      A(KIXESS), A(KIXNSS), A(KLTESS), A(KLTNSS), *100)
      END IF

C        Draw elements as spheres, if requested

      IF (SPHPLT .NE. 0) THEN
         CALL MDFIND( 'LENE',   KLENE,  LDUM)
         CALL MDFIND( 'LINK',   KLINK,  LDUM)
         CALL MDFIND( 'NUMLNK', KNUMLN, LDUM)
         CALL MDFIND( 'NUMATR', KNUMAT, LDUM)
         CALL MDFIND( 'ATRIB',  KATRIB, LDUM)
         IF (FIRST) THEN
           CALL MDRSRV( 'ISPSOR', KSPSOR, NUMEL)
           CALL MDRSRV( 'SPRAD',  KSPRAD, NUMNP)
           CALL MDRSRV( 'ISPBLK', KSPBLK, NUMNP)
           FIRST = .FALSE.
         ELSE
           CALL MDFIND( 'ISPSOR', KSPSOR, LDUM)
           CALL MDFIND( 'SPRAD',  KSPRAD, LDUM)
           CALL MDFIND( 'ISPBLK', KSPBLK, LDUM)
         END IF
         CALL MDSTAT (NERR, MEM)
         IF (NERR .GT. 0) GOTO 100
         IF (SPHPLT .GE. 1) THEN
            FFLAG = 'NOFILL'
         ELSE
            FFLAG = 'FILL  '
         ENDIF
         IF (is3dim .and. IHIDOP .GE. 5 .AND.
     &        MODDET .NE. 'CONTOUR') then
           CALL MDFIND('SHDCOL', KSHDCL, LEN)
           CALL MDFIND('ISHDCL', KISHCL, LEN)
           CALL SHDSPH ( A(KLENE), A(KLINK), A(KNUMLN), A(KNUMAT),
     &       XN, YN, ZN, A(KATRIB), BLKCOL, IDELB,
     &       A(KSPSOR), A(KSPRAD), IELBST,
     *       A(KSPBLK), A(KSHDCL), A(KISHCL), HIDENP, *100)
         ELSE
           CALL SYMSPH ( A(KLENE), A(KLINK), A(KNUMLN), A(KNUMAT),
     &       XN, YN, ZN, A(KATRIB), BLKCOL, FFLAG, IDELB, VARNP,
     &       MODDET, A(KSPSOR), A(KSPRAD), IELBST,
     *       A(KSPBLK), HIDENP, *100)
         ENDIF
       END IF
      RETURN

  100 CONTINUE
      RETURN 1
      END
