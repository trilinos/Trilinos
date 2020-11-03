C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE DTPLT1 (A, MSHNUM, MSHLIN, MLNTYP, MODDET, MODTYP,
     &   LENF, NLNKF, LINKF, LENL, LINSET,
     &   HIDENP, HIDEF, DOIXF, NXFAC, IXFAC,
     &   XN, YN, ZN, XF, YF, ZF,
     &   IELBST, ISEVOK, IN2ELB, IVN2B, DODEAD, IDN2B,
     &   NNPSET, ISSNPS, NESSET, ISSESS,
     &   IDTVAR, VARNP, VARFAC, NMIN, NMAX, FMIN, FMAX, VECMAX, ISVOK,
     &   BLKCOL, IDELB, IHIDOP, MAPEL, MAPND, *)
C=======================================================================

C   --*** DTPLT1 *** (DETOUR) Plot one view
C   --   Modified by John Glick - 11/28/88
C   --   Written by Amy Gilkey - revised 04/11/88
C   --   Modified version 1.1a  - November 1990 - R.J. Meyers
C   --           added color coded sphere capability
C   --
C   --DTPLT1 does the plotting for a single view.  The labeling of the view
C   --is done elsewhere.
C   --
C   --Parameters:
C   --   A - IN - the dynamic memory array
C   --   MSHNUM - IN - the mesh numbering (as in /MSHOPT/)
C   --   MSHLIN - IN - the display type for the mesh lines (as in /MSHOPT/)
C   --   MLNTYP - IN - the line type of lines (as in /MSHOPT/)
C   --   MODDET - IN - the display mode (as in /DETOPT/)
C   --   MODTYP - IN - the display mode type (as in /DETOPT/)
C   --   LENF - IN - the cumulative face counts by element block
C   --   NLNKF - IN - the number of nodes per face
C   --   LINKF - IN - the connectivity for all faces
C   --   LENL - IN - the cumulative line counts by element block
C   --   LINSET - IN - the sorted line set
C   --   HIDENP(i) - IN - true iff node i is hidden (3D only)
C   --   HIDEF(i) - IN - true iff face i is hidden (3D only)
C   --   DOIXF - IN - IXFAC valid iff true
C   --   NXFAC - IN - the number of ordered faces (if DOIXF)
C   --   IXFAC - IN - the indices of the ordered faces (if DOIXF)
C   --   XN, YN, ZN - IN - the nodal coordinates
C   --   XF, YF, ZF - IN - the face center coordinates
C   --   IELBST - IN - the element block status (>0 if selected)
C   --   ISEVOK - IN - the element block variable truth table;
C   --      variable i of block j exists iff ISEVOK(j,i)
C   --   IN2ELB - IN - the element block for each node;
C   --      <0 if not in any selected element block
C   --      =0 if in more than one selected element block
C   --   IVN2B - IN - the element block for each node for the contour variable;
C   --      <0 if not in any selected element block
C   --      =0 if in more than one selected element block
C   --   DODEAD - IN - true iff dead nodes are to be displayed
C   --   IDN2B - IN - the element block for each dead node; dead if >= 0
C   --   NNPSET - IN - the number of selected node sets
C   --   ISSNPS - IN - the indices of the selected node sets
C   --   NESSET - IN - the number of selected side sets
C   --   ISSESS - IN - the indices of the selected side sets
C   --   IDTVAR - IN - the variable numbers, if any
C   --   VARNP - IN - the nodal variable values
C   --   VARFAC - IN - the face variable values
C   --   NMIN, NMAX - IN - the number of variables values matching the
C   --      minimum and the maximum (for contour plots only)
C   --   FMIN, FMAX - IN - the nodal variable minimum and maximum
C   --      (for contour modes only)
C   --   VECMAX - IN - the vector maximum, scaled (for vector only)
C   --   ISVOK - SCRATCH - size = NELBLK, only if any variables
C   --   BLKCOL - IN - the user selected colors of the element blocks.
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
C   --   * - return statement if cancel function active
C   --
C   --Common Variables:
C   --   Uses NELBLK of /DBNUMS/
C   --   Uses IS3DIM, NUMNPF of /D3NUMS/

      include 'dbnums.blk'
      include 'dbnumgq.blk'
      include 'd3nums.blk'

      DIMENSION A(*)
      CHARACTER*(*) MSHNUM
      INTEGER MLNTYP(-1:1)
      CHARACTER*8 MODDET, MODTYP
      INTEGER LENF(0:NELBLK)
      INTEGER NLNKF(NELBLK)
      INTEGER LINKF(*)
      INTEGER LENL(-2:NELBLK), LINSET(LLNSET,*)
      LOGICAL HIDENP(NUMNPF)
      LOGICAL HIDEF(*)
      LOGICAL DOIXF
      INTEGER IXFAC(*)
      REAL XN(NUMNPF), YN(NUMNPF), ZN(NUMNPF)
      REAL XF(*), YF(*), ZF(*)
      INTEGER IELBST(NELBLK)
      LOGICAL ISEVOK(NELBLK,*)
      INTEGER IN2ELB(NUMNPF)
      INTEGER IVN2B(NUMNPF)
      LOGICAL DODEAD
      INTEGER IDN2B(NUMNPF)
      INTEGER ISSNPS(*)
      INTEGER ISSESS(*)
      INTEGER IDTVAR(*)
      REAL VARNP(*), VARFAC(*)
      REAL VECMAX
      LOGICAL ISVOK(NELBLK)
      INTEGER BLKCOL(0:NELBLK)
      INTEGER IDELB(*)
      INTEGER MAPEL(*), MAPND(*)

      LOGICAL WIDLIN
      LOGICAL LDUM

C   --Call appropriate paint function routine, if any; mesh is drawn
C   --over paint

      IF (MODDET .EQ. 'SOLID') THEN
        IF (DOIXF) THEN
          IF (is3dim .and. IHIDOP .EQ. 5) then
            CALL MDFIND('SHDCOL', KSHDCL, LEN)
            CALL MDFIND('ISHDCL', KISHCL, LEN)
            CALL SHADE (LENF, NLNKF, LINKF, NXFAC, IXFAC,
     &        XN, YN, ZN, IELBST, BLKCOL, IDELB,
     *        A(KSHDCL), A(KISHCL), IHIDOP, *100)
          else if (is3dim .and. ihidop .ge. 6) then
            CALL WRTRAY (LENF, NLNKF, LINKF, NXFAC, IXFAC,
     &        XN, YN, ZN, IELBST, BLKCOL, IDELB, *100)
          ELSE
            CALL SOLID (LENF, NLNKF, LINKF, NXFAC, IXFAC,
     &        XN, YN, ZN, IELBST, BLKCOL, IDELB, *100)
          END IF
        ELSE
          CALL QSOLID (LENF, NLNKF, LINKF, HIDEF,
     &      XN, YN, ZN, IELBST, BLKCOL, IDELB, *100)
        END IF

      ELSE IF (MODDET .EQ. 'CONTOUR') THEN
         IF (MODTYP .EQ. 'PAINT') THEN
C         --Use the selected color table
            CALL GRCOLU ('ALTERNATE')

            CALL EVAROK (1, IDTVAR, NELBLK, IELBST, ISEVOK, ISVOK)
            IF (DOIXF) THEN
               CALL PAINT (VARNP, LENF, NLNKF, LINKF, NXFAC, IXFAC,
     &            XN, YN, ZN, XF, YF, ZF, ISVOK, FMIN, FMAX, *100)
            ELSE
               CALL QPAINT (VARNP, LENF, NLNKF, LINKF, HIDEF,
     &            XN, YN, ZN, XF, YF, ZF, ISVOK, FMIN, FMAX, *100)
            END IF

C         --Use the standard color table
            CALL GRCOLU ('STANDARD')
         END IF

      ELSE IF (MODDET .EQ. 'ELEMCONT') THEN
         IF (MODTYP .EQ. 'PAINT') THEN
C         --Use the selected color table
            CALL GRCOLU ('ALTERNATE')

            CALL EVAROK (1, IDTVAR, NELBLK, IELBST, ISEVOK, ISVOK)
            IF (DOIXF) THEN
               CALL EPAINT (VARFAC, LENF, NLNKF, LINKF, NXFAC, IXFAC,
     &            XN, YN, ZN, ISVOK, FMIN, FMAX, *100)
            ELSE
               CALL QEPAIN (VARFAC, LENF, NLNKF, LINKF, HIDEF,
     &            XN, YN, ZN, ISVOK, FMIN, FMAX, *100)
            END IF

C         --Use the standard color table
            CALL GRCOLU ('STANDARD')
         END IF
      END IF

C   --Draw mesh (with numbering, etc.)

      WIDLIN = (MODDET .EQ. 'WIREFRAM')

      CALL MSPLT1 (A, WIDLIN, MSHNUM, MSHLIN, MLNTYP,
     &   LENF, NLNKF, LINKF, LENL, LINSET,
     &   HIDENP, HIDEF, XN, YN, ZN, XF, YF, ZF,
     &   IELBST, IN2ELB, DODEAD, IDN2B,
     &   .FALSE., 0, IDUM,
     &   NNPSET, ISSNPS, NESSET, ISSESS, BLKCOL,
     &   IDELB, VARNP, MODDET, IHIDOP, MAPEL, MAPND, *100)

C   --Call appropriate function routine

      IF (MODDET .EQ. 'CONTOUR') THEN
         IF (MODTYP .EQ. 'LINE') THEN
C         --Use the selected color table
            CALL GRCOLU ('ALTERNATE')

            CALL EVAROK (1, IDTVAR, NELBLK, IELBST, ISEVOK, ISVOK)
            CALL CONTOR (VARNP, LENF, NLNKF, LINKF, HIDEF,
     &         XN, YN, ZN, XF, YF, ZF, LENL, LINSET,
     &         IVN2B, ISVOK, FMIN, FMAX, *100)

C         --Use the standard color table
            CALL GRCOLU ('STANDARD')
         END IF

         CALL MRKNOD (VARNP, HIDENP,
     &      XN, YN, ZN, IVN2B, NMIN, NMAX, FMIN, FMAX,
     &      BLKCOL,  *100)

      ELSE IF (MODDET .EQ. 'ELEMCONT') THEN
         CALL MRKFAC (LENF(NELBLK), VARFAC, HIDEF, XF, YF, ZF,
     &      NMIN, NMAX, FMIN, FMAX, BLKCOL, *100)

      ELSE IF (MODDET .EQ. 'VECTOR') THEN
         IF ((MODTYP .EQ. 'NODE') .OR. (MODTYP .EQ. 'ELEMENT')) THEN
            N = NDIM
         ELSE
            N = 3
         END IF
         IF (MODTYP .EQ. 'NODE') THEN
            CALL VECTORN (MODTYP, VARNP, NUMNPF,
     &         HIDENP, XN, YN, ZN, IN2ELB, VECMAX,
     &         BLKCOL,  IDELB, *100)
         ELSE
            CALL EVAROK (N, IDTVAR, NELBLK, IELBST, ISEVOK, ISVOK)
            CALL VECTOR (MODTYP, VARFAC, LENF(NELBLK), LENF, NLNKF,
     &         HIDEF, XF, YF, ZF, IDUM, ISVOK, VECMAX, BLKCOL,
     &         IDELB, *100)
         END IF

      ELSE IF (MODDET .EQ. 'SYMBOL') THEN
         CALL EVAROK (1, IDTVAR, NELBLK, IELBST, ISEVOK, ISVOK)
         IF (MODTYP .NE. 'STATE') THEN
            CALL SYMBOL_BL (MODTYP, VARFAC, LENF, NLNKF, HIDEF,
     &         XF, YF, ZF, ISVOK, BLKCOL, IDELB, *100)
         ELSE
            CALL ELESTA (MODTYP, VARFAC, LENF, NLNKF, LINKF, HIDEF,
     &         XN, YN, ZN, ISVOK, *100)
         END IF

      ELSE IF (MODDET .EQ. 'GAUSS') THEN
         CALL EVAROK (4, IDTVAR, NELBLK, IELBST, ISEVOK, ISVOK)
         LVARF = LENF(NELBLK)
         CALL GAUSS_BL (MODTYP, VARFAC, LENF, NLNKF, LINKF, HIDEF,
     &      XN, YN, ZN, ISVOK, LVARF, *100)
      END IF

      RETURN

  100 CONTINUE
      RETURN 1
      END
