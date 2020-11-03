C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE DBLIST (A, INLINE,
     &   VERB, IFLD, INTYP, CFIELD, IFIELD, RFIELD,
     &   NAMECO, NAMELB, NAMES, QAREC, INFREC,
     &   TIMES, WHOTIM, NPTIMS, IPTIMS, XN, YN, ZN,
     &   IDELB, LENE, NLNKE, LINKE, ISEVOK,
     &   IDNPS, NNNPS, IDESS, NEESS, NNESS,
     &   NCSTEP, LISNP, NLISEL, LISEL, LISNPS, LISESS,
     &   LISHV, LISGV, LISNV, LISEV,
     $   EBNAME, NSNAME, SSNAME, NAMLEN, MAPEL, MAPND)
C=======================================================================

C   --*** DBLIST *** (BLOT) Process LIST/PRINT commands
C   --   Written by Amy Gilkey - revised 05/12/88
C   --
C   --DBLIST processes a LIST or PRINT command.
C   --
C   --This routine uses MDFIND to find the following dynamic memory arrays:
C   --   MAPEL - the element order map
C   --   NUMATR - the number of attributes for each element block
C   --   ATRIB - the attributes for each element block
C   --   IXNNPS - the index of the first node for each node set
C   --   LTNNPS - the nodes for all node sets
C   --   FACNPS - the distribution factors for all node sets
C   --   IXEESS - the index of the first element for each side set
C   --   IXNESS - the index of the first node for each side set
C   --   LTEESS - the elements for all side sets
C   --   LTNESS - the nodes for all side sets
C   --   FACESS - the distribution factors for all side sets
C   --
C   --Parameters:
C   --   A - IN - the dynamic memory base array
C   --   INLINE - IN/OUT - the parsed input lines for the log file
C   --   VERB - IN/OUT - the command verb
C   --   IFLD, INTYP, CFIELD, IFIELD, RFIELD - IN/OUT - the free-field
C   --      reader index and fields
C   --   NAMECO - IN - the names of the coordinates
C   --   NAMELB - IN - the names of the element block types
C   --   NAMES - IN - the variable names
C   --   QAREC - IN - the QA records containing:
C   --      (1) - the analysis code name
C   --      (2) - the analysis code QA descriptor
C   --      (3) - the analysis date
C   --      (4) - the analysis time
C   --   INFREC - IN - the information records
C   --   TIMES - IN - the times for all time steps
C   --   WHOTIM - IN - true iff whole (versus history) time step
C   --   NPTIMS - IN - the number of selected time steps
C   --   IPTIMS - IN - the selected time steps
C   --   XN, YN, ZN - IN - the nodal coordinates
C   --   IDELB - IN - the element block ID for each block
C   --   LENE - IN - the cumulative element counts by element block
C   --   NLNKE - IN - the number of nodes per element for each block
C   --   LINKE - IN - the connectivity array for all blocks
C   --   ISEVOK - IN - the element block variable truth table
C   --   IDNPS - IN - the node set ID for each set
C   --   NNNPS - IN - the number of nodes for each set
C   --   IDESS - IN - the side set ID for each set
C   --   NEESS - IN - the number of elements for each set
C   --   NNESS - IN - the number of nodes for each set
C   --   NCSTEP - IN/OUT - the current step number for display
C   --   LISNP - IN/OUT - the indices of the selected nodes
C   --   NLISEL - IN/OUT - the number of selected elements for each block
C   --   LISEL - IN/OUT - the indices of the selected elements (by block)
C   --   LISNPS - IN/OUT - the indices of the selected node sets
C   --   LISESS - IN/OUT - the indices of the selected side sets
C   --   LISHV - IN/OUT - the indices of the selected history variables
C   --   LISGV - IN/OUT - the indices of the selected global variables
C   --   LISNV - IN/OUT - the indices of the selected nodal variables
C   --   LISEV - IN/OUT - the indices of the selected element variables
C   --
C   --Common Variables:
C   --   Uses NDB of /DBASE/
C   --   Uses TITLE of /DBTITL/
C   --   Uses NDIM, NUMNP, NUMEL, NELBLK, NVARHI, NVARGL, NVARNP, NVAREL,
C   --      NSTEPS of /DBNUMS/
C   --   Uses NUMNPS, LNPSNL, NUMESS, LESSEL, LESSNL of /DBNUMG/
C   --   Uses NQAREC, NINFO of /DBNUMQ/
C   --   Uses NOUT, NCRT, NPRT, ANYPRT of /OUTFIL/

      include 'dbname.blk'
      include 'params.blk'
      include 'progqa.blk'
      include 'dbase.blk'
      include 'dbtitl.blk'
      include 'dbnums.blk'
      include 'dbnumgq.blk'
      include 'outfil.blk'

      DIMENSION A(*)
      CHARACTER*(*) INLINE(*)
      CHARACTER*(*) VERB
      INTEGER     INTYP(*)
      CHARACTER*(*) CFIELD(*)
      INTEGER     IFIELD(*)
      REAL        RFIELD(*)
      CHARACTER*(MXSTLN) NAMECO(*), NAMELB(*)
      CHARACTER*(NAMLEN) NAMES(*)
      CHARACTER*(NAMLEN) EBNAME(*), NSNAME(*), SSNAME(*)
      CHARACTER*(MXSTLN) QAREC(4,*)
      CHARACTER*(MXLNLN) INFREC(*)
      CHARACTER*2048 FILNAM, ERRMSG

      REAL TIMES(*)
      LOGICAL WHOTIM(*)
      INTEGER IPTIMS(*)
      REAL XN(*), YN(*), ZN(*)
      INTEGER IDELB(*), NLNKE(*)
      INTEGER LENE(0:*), LINKE(*)
      LOGICAL ISEVOK(*)
      INTEGER IDNPS(*), NNNPS(*)
      INTEGER IDESS(*), NEESS(*), NNESS(*)
      INTEGER LISNP(0:*)
      INTEGER NLISEL(0:*), LISEL(0:*)
      INTEGER LISNPS(0:*), LISESS(0:*)
      INTEGER LISHV(0:*), LISGV(0:*), LISNV(0:*), LISEV(0:*)
      INTEGER MAPEL(*)
      INTEGER MAPND(*)

      LOGICAL FFEXST

      CHARACTER*(MXNAME) WORD, LISTYP
      CHARACTER*(MXNAME) DUMLIN
      CHARACTER*(MXSTLN) OPT

      LOGICAL FIRST
      SAVE FIRST
C      --FIRST - true iff first time through routine

      CHARACTER*(MXSTLN) CMDTBL(32)
      SAVE CMDTBL
C      --CMDTBL - the valid commands table

      DATA FIRST / .TRUE. /

C   --Command table follows.  Remember to change the dimensioned size when
C   --changing the table.
      DATA CMDTBL /
     1   'TITLE   ', 'VARS    ',
     2   'COORDINA', 'MAP     ', 'NODEMAP ', 'NMAP    ',
     3   'BLOCKS  ', 'MATERIAL', 'LINK    ', 'CONNECTI', 'ATTRIBUT',
     4   'NSETS   ', 'NNODES  ', 'NFACTORS',
     5   'SSETS   ', 'SELEMS  ', 'SNODES  ', 'SFACTORS',
     6   'QA      ', 'NAMES   ',
     7   'HVARS   ', 'GVARS   ', 'NVARS   ', 'EVARS   ',
     8   'HISTORY ', 'GLOBALS ', 'NODALS  ', 'ELEMENTS',
     9   'STEPS   ', 'TIMES   ', 'MINMAX  ',
     1   '        ' /

C *** Initialization ***

C   --Initialize parameters first time through, then reset

      IF (FIRST) THEN

C      --Set up the print file

         NCRT = -1
         NOUT = NCRT
         NPRT = 21
         ANYPRT = .FALSE.

C      --Reset selection if not done by SELECT command

         CALL DBSEL (A, A, INLINE,
     &      'reset', IFLD, INTYP, CFIELD, IFIELD, RFIELD,
     &      NAMES, TIMES, WHOTIM, NPTIMS, IPTIMS,
     &      IDELB, LENE, IDNPS, IDESS,
     &      NCSTEP, LISNP, NLISEL, LISEL, LISNPS, LISESS,
     &      LISHV, LISGV, LISNV, LISEV, MAPEL, MAPND)

         FIRST = .FALSE.
      END IF

C   --Get the command verb

      IF (VERB .EQ. 'LIST') THEN
         NOUT = NCRT

      ELSE IF (VERB .EQ. 'PRINT') THEN
         IF (NPRT .LE. 0) THEN
            CALL PRTERR ('CMDERR', 'Print file cannot be opened')
            GOTO 150
         END IF

         IF (.NOT. ANYPRT) THEN
           filnam = basenam(:lenstr(basenam)) // '.lis'
           open (unit=nprt, file=filnam(:lenstr(filnam)),
     *       form='formatted', status='unknown', iostat=ierr)
           IF (IERR .NE. 0) THEN
             ERRMSG = 'Print file "'//FILNAM(:LENSTR(FILNAM))//
     *         '" could not be opened.'
             CALL PRTERR ('CMDERR', ERRMSG(:LENSTR(ERRMSG)))
             NPRT = -1
             GOTO 150
           END IF

           CALL BANNER (NPRT, QAINFO,
     &       ' ', ' ', ' ')

           CALL PRINIT ('N', NPRT, NDB, DBNAME, TITLE,
     &       NDIM, NUMNP, NUMEL, NELBLK,
     &       NUMNPS, LNPSNL, LNPSDF, NUMESS, LESSEL, LESSNL,
     &       NVARGL, NVARNP, NVAREL, NVARNS, NVARSS)

           ANYPRT = .TRUE.
         END IF

         NOUT = NPRT

      END IF

      CALL FFADDC (VERB, INLINE(1))

      CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)
      CALL ABRSTR (LISTYP, WORD, CMDTBL)
      IF (LISTYP .EQ. ' ') LISTYP = WORD

C *** GENESIS Print Commands ***

      IF (FFEXST (IFLD, INTYP)) THEN

         IF (LISTYP .EQ. 'COORDINA') THEN
            CALL PRTERR ('CMDREQ',
     &         'Please use a SELECT NODES command')
            CALL DBSEL (A, A, DUMLIN,
     &         'NODES', IFLD, INTYP, CFIELD, IFIELD, RFIELD,
     &         NAMES, TIMES, WHOTIM, NPTIMS, IPTIMS,
     &         IDELB, LENE, IDNPS, IDESS,
     &         NCSTEP, LISNP, NLISEL, LISEL, LISNPS, LISESS,
     &         LISHV, LISGV, LISNV, LISEV, MAPEL, MAPND)

         ELSE IF ((LISTYP .EQ. 'BLOCKS')
     &      .OR. (LISTYP .EQ. 'MATERIAL')) THEN
            CALL PRTERR ('CMDREQ',
     &         'Please use a SELECT BLOCKS command')
            CALL DBSEL (A, A, DUMLIN,
     &         'BLOCKS', IFLD, INTYP, CFIELD, IFIELD, RFIELD,
     &         NAMES, TIMES, WHOTIM, NPTIMS, IPTIMS,
     &         IDELB, LENE, IDNPS, IDESS,
     &         NCSTEP, LISNP, NLISEL, LISEL, LISNPS, LISESS,
     &         LISHV, LISGV, LISNV, LISEV, MAPEL, MAPND)

         ELSE IF ((LISTYP .EQ. 'LINK') .OR. (LISTYP .EQ. 'CONNECTI')
     &      .OR. (LISTYP .EQ. 'ATTRIBU')) THEN
            CALL PRTERR ('CMDREQ',
     &         'Please use a SELECT ELEMENTS command')
            CALL DBSEL (A, A, DUMLIN,
     &         'ELEMENTS', IFLD, INTYP, CFIELD, IFIELD, RFIELD,
     &         NAMES, TIMES, WHOTIM, NPTIMS, IPTIMS,
     &         IDELB, LENE, IDNPS, IDESS,
     &         NCSTEP, LISNP, NLISEL, LISEL, LISNPS, LISESS,
     &         LISHV, LISGV, LISNV, LISEV, MAPEL, MAPND)

         ELSE IF (LISTYP .EQ. 'NSETS') THEN
            CALL PRTERR ('CMDREQ',
     &         'Please use a SELECT NSETS command')
            CALL DBSEL (A, A, DUMLIN,
     &         'NSETS', IFLD, INTYP, CFIELD, IFIELD, RFIELD,
     &         NAMES, TIMES, WHOTIM, NPTIMS, IPTIMS,
     &         IDELB, LENE, IDNPS, IDESS,
     &         NCSTEP, LISNP, NLISEL, LISEL, LISNPS, LISESS,
     &         LISHV, LISGV, LISNV, LISEV, MAPEL, MAPND)

         ELSE IF (LISTYP .EQ. 'SSETS') THEN
            CALL PRTERR ('CMDREQ',
     &         'Please use a SELECT SSETS command')
            CALL DBSEL (A, A, DUMLIN,
     &         'SSETS', IFLD, INTYP, CFIELD, IFIELD, RFIELD,
     &         NAMES, TIMES, WHOTIM, NPTIMS, IPTIMS,
     &         IDELB, LENE, IDNPS, IDESS,
     &         NCSTEP, LISNP, NLISEL, LISEL, LISNPS, LISESS,
     &         LISHV, LISGV, LISNV, LISEV, MAPEL, MAPND)

         ELSE IF (LISTYP .EQ. 'HVARS') THEN
            CALL PRTERR ('CMDREQ',
     &         'Please use a SELECT HVARS command')
            CALL DBSEL (A, A, DUMLIN,
     &         LISTYP, IFLD, INTYP, CFIELD, IFIELD, RFIELD,
     &         NAMES, TIMES, WHOTIM, NPTIMS, IPTIMS,
     &         IDELB, LENE, IDNPS, IDESS,
     &         NCSTEP, LISNP, NLISEL, LISEL, LISNPS, LISESS,
     &         LISHV, LISGV, LISNV, LISEV, MAPEL, MAPND)

         ELSE IF (LISTYP .EQ. 'GVARS') THEN
            CALL PRTERR ('CMDREQ',
     &         'Please use a SELECT GVARS command')
            CALL DBSEL (A, A, DUMLIN,
     &         LISTYP, IFLD, INTYP, CFIELD, IFIELD, RFIELD,
     &         NAMES, TIMES, WHOTIM, NPTIMS, IPTIMS,
     &         IDELB, LENE, IDNPS, IDESS,
     &         NCSTEP, LISNP, NLISEL, LISEL, LISNPS, LISESS,
     &         LISHV, LISGV, LISNV, LISEV, MAPEL, MAPND)

         ELSE IF (LISTYP .EQ. 'NVARS') THEN
            CALL PRTERR ('CMDREQ',
     &         'Please use a SELECT NVARS command')
            CALL DBSEL (A, A, DUMLIN,
     &         LISTYP, IFLD, INTYP, CFIELD, IFIELD, RFIELD,
     &         NAMES, TIMES, WHOTIM, NPTIMS, IPTIMS,
     &         IDELB, LENE, IDNPS, IDESS,
     &         NCSTEP, LISNP, NLISEL, LISEL, LISNPS, LISESS,
     &         LISHV, LISGV, LISNV, LISEV, MAPEL, MAPND)

         ELSE IF (LISTYP .EQ. 'EVARS') THEN
            CALL PRTERR ('CMDREQ',
     &         'Please use a SELECT EVARS command')
            CALL DBSEL (A, A, DUMLIN,
     &         LISTYP, IFLD, INTYP, CFIELD, IFIELD, RFIELD,
     &         NAMES, TIMES, WHOTIM, NPTIMS, IPTIMS,
     &         IDELB, LENE, IDNPS, IDESS,
     &         NCSTEP, LISNP, NLISEL, LISEL, LISNPS, LISESS,
     &         LISHV, LISGV, LISNV, LISEV, MAPEL, MAPND)

         ELSE IF (LISTYP .EQ. 'NODALS') THEN
            CALL PRTERR ('CMDREQ',
     &         'Please use a SELECT NODES command')
            CALL DBSEL (A, A, DUMLIN,
     &         'NODES', IFLD, INTYP, CFIELD, IFIELD, RFIELD,
     &         NAMES, TIMES, WHOTIM, NPTIMS, IPTIMS,
     &         IDELB, LENE, IDNPS, IDESS,
     &         NCSTEP, LISNP, NLISEL, LISEL, LISNPS, LISESS,
     &         LISHV, LISGV, LISNV, LISEV, MAPEL, MAPND)

         ELSE IF (LISTYP .EQ. 'ELEMENTS') THEN
            CALL PRTERR ('CMDREQ',
     &         'Please use a SELECT ELEMENTS command')
            CALL DBSEL (A, A, DUMLIN,
     &         'XELEM', IFLD, INTYP, CFIELD, IFIELD, RFIELD,
     &         NAMES, TIMES, WHOTIM, NPTIMS, IPTIMS,
     &         IDELB, LENE, IDNPS, IDESS,
     &         NCSTEP, LISNP, NLISEL, LISEL, LISNPS, LISESS,
     &         LISHV, LISGV, LISNV, LISEV, MAPEL, MAPND)
         END IF
      END IF

      IF ((LISTYP .EQ. 'TITLE') .OR. (LISTYP .EQ. 'VARS')) THEN
         CALL FFADDC (LISTYP, INLINE(1))

         IF (EXODUS) THEN
            CALL PRINIT ('NTSIV', NOUT, NDB, DBNAME, TITLE,
     &         NDIM, NUMNP, NUMEL, NELBLK,
     &         NUMNPS, LNPSNL, LNPSDF, NUMESS, LESSEL, LESSNL,
     &         NVARGL, NVARNP, NVAREL, NVARNS, NVARSS)
         ELSE
            CALL PRINIT ('NTSI', NOUT, NDB, DBNAME, TITLE,
     &         NDIM, NUMNP, NUMEL, NELBLK,
     &         NUMNPS, LNPSNL, LNPSDF, NUMESS, LESSEL, LESSNL,
     &         NVARGL, NVARNP, NVAREL, NVARNS, NVARSS)
         END IF

      ELSE IF (LISTYP .EQ. 'COORDINA') THEN
         CALL FFADDC (LISTYP, INLINE(1))
         CALL CKNONE (NUMNP, .FALSE., 'nodes', *150)
         CALL CKNONE (NDIM, .FALSE., 'coordinates', *150)
         CALL CKNONE (LISNP(0), .TRUE., 'coordinates', *150)

         CALL PRXYZ ('*', NOUT, NDIM, NAMECO, NUMNP, LISNP, XN, YN, ZN,
     &     MAPND)

      ELSE IF (LISTYP .EQ. 'MAP') THEN
         CALL FFADDC (LISTYP, INLINE(1))
         CALL CKNONE (NUMEL, .FALSE., 'elements', *150)

         CALL PRMAP ('*', NOUT, NUMEL, MAPEL, 'Element')

      ELSE IF ((LISTYP .EQ. 'NODEMAP') .OR. (LISTYP .EQ. 'NMAP')) THEN
         CALL FFADDC (LISTYP, INLINE(1))
         CALL CKNONE (NUMNP, .FALSE., 'nodes', *150)

         CALL PRMAP ('*', NOUT, NUMNP, MAPND, 'Node')

      ELSE IF ((LISTYP .EQ. 'BLOCKS') .OR. (LISTYP .EQ. 'MATERIAL')
     &   .OR. (LISTYP .EQ. 'LINK') .OR. (LISTYP .EQ. 'CONNECTI')
     &   .OR. (LISTYP .EQ. 'ATTRIBUT')) THEN
         CALL FFADDC (LISTYP, INLINE(1))
         CALL CKNONE (NELBLK, .FALSE., 'element blocks', *150)
         CALL CKNONE (NLISEL(0), .TRUE., 'element blocks', *150)

         CALL MDFIND ('NUMATR', KNATR, IDUM)
         CALL MDSTAT (NERR, MEM)
         IF (NERR .GT. 0) GOTO 150

         IF ((LISTYP .EQ. 'BLOCKS') .OR. (LISTYP .EQ. 'MATERIAL')) THEN
            IF (EXODUS) THEN
               OPT = 'NV'
            ELSE
               OPT = 'N'
            END IF
         ELSE IF ((LISTYP .EQ. 'LINK')
     &      .OR. (LISTYP .EQ. 'CONNECTI')) THEN
            OPT = 'C'
         ELSE IF (LISTYP .EQ. 'ATTRIBUT') THEN
            OPT = 'A'
            CALL MDFIND ('ATRIB', KATRIB, N)
            CALL MDSTAT (NERR, MEM)
            IF (NERR .GT. 0) GOTO 150
         END IF

         IF ((LISTYP .EQ. 'LINK')
     &      .OR. (LISTYP .EQ. 'CONNECTI')) THEN
            IF (NDIM .EQ. 2) THEN
               DO 100 I = 1, NELBLK
                  IF (NLNKE(I) .EQ. 8) THEN
                     CALL PRTERR ('WARNING',
     &                  'Connectivity has been reordered')
                     GOTO 110
                  END IF
  100          CONTINUE
  110          CONTINUE
            END IF
         END IF

         CALL DBVIX_BL ('E', 1, IXEV)
         IF ((OPT .EQ. '*') .OR. (INDEX (OPT, 'V') .GT. 0)) THEN
            CALL MDRSRV ('XLISEV', KXLSEV, NVAREL)
            CALL MDSTAT (NERR, MEM)
            IF (NERR .GT. 0) GOTO 150
         END IF

         CALL PRELB (OPT, NOUT, NELBLK, NLISEL, LISEL,
     &      IDELB, LENE, NLNKE, A(KNATR), LINKE, A(KATRIB),
     &      NAMELB, EBNAME, NVAREL, NAMES(IXEV), ISEVOK, A(KXLSEV),
     *      MAPEL, MAPND)

         IF ((OPT .EQ. '*') .OR. (INDEX (OPT, 'V') .GT. 0)) THEN
            CALL MDDEL ('XLISEV')
         END IF

      ELSE IF ((LISTYP .EQ. 'NSETS')
     &   .OR. (LISTYP .EQ. 'NNODES')
     &   .OR. (LISTYP .EQ. 'NFACTORS')) THEN
         CALL FFADDC (LISTYP, INLINE(1))
         CALL CKNONE (NUMNPS, .FALSE., 'node sets', *150)
         CALL CKNONE (LISNPS(0), .TRUE., 'node sets', *150)

         CALL MDFIND ('IXNNPS', KIXNNS, IDUM)
         CALL MDSTAT (NERR, MEM)
         IF (NERR .GT. 0) GOTO 150

         IF (LISTYP .EQ. 'NSETS') THEN
            OPT = ' '
         ELSE IF (LISTYP .EQ. 'NNODES') THEN
            OPT = 'N'
            CALL MDFIND ('LTNNPS', KLTNNS, N)
            CALL MDSTAT (NERR, MEM)
            IF (NERR .GT. 0) GOTO 150
         ELSE IF (LISTYP .EQ. 'NFACTORS') THEN
            OPT = 'F'
            CALL MDFIND ('FACNPS', KFACNS, N)
            CALL MDSTAT (NERR, MEM)
            IF (NERR .GT. 0) GOTO 150
         END IF

         CALL PRNPS (OPT, NOUT, NUMNPS, LISNPS, LNPSNL,
     &      IDNPS, NNNPS, A(KIXNNS), A(KLTNNS), A(KFACNS),
     $      NSNAME, MAPND)

      ELSE IF ((LISTYP .EQ. 'SSETS')
     &   .OR. (LISTYP .EQ. 'SELEMS')
     &   .OR. (LISTYP .EQ. 'SNODES')
     &   .OR. (LISTYP .EQ. 'SFACTORS')) THEN
         CALL FFADDC (LISTYP, INLINE(1))
         CALL CKNONE (NUMESS, .FALSE., 'side sets', *150)
         CALL CKNONE (LISESS(0), .TRUE., 'side sets', *150)

         CALL MDFIND ('IXEESS', KIXESS, IDUM)
         CALL MDFIND ('IXNESS', KIXNSS, IDUM)
         CALL MDSTAT (NERR, MEM)
         IF (NERR .GT. 0) GOTO 150

         IF (LISTYP .EQ. 'SSETS') THEN
            OPT = ' '
         ELSE IF (LISTYP .EQ. 'SELEMS') THEN
            OPT = 'E'
            CALL MDFIND ('LTEESS', KLTESS, N)
            CALL MDSTAT (NERR, MEM)
            IF (NERR .GT. 0) GOTO 150
         ELSE IF (LISTYP .EQ. 'SNODES') THEN
            OPT = 'N'
            CALL MDFIND ('LTNESS', KLTNSS, N)
            CALL MDSTAT (NERR, MEM)
            IF (NERR .GT. 0) GOTO 150
         ELSE IF (LISTYP .EQ. 'SFACTORS') THEN
            OPT = 'F'
            CALL MDFIND ('FACESS', KFACSS, N)
            CALL MDSTAT (NERR, MEM)
            IF (NERR .GT. 0) GOTO 150
         END IF

         if (nness(1) .lt. 0) then
            call getssn(a, ierr)
         end if

         CALL PRESS (OPT, NOUT, NUMESS, LISESS, LESSEL, LESSNL,
     &      IDESS, NEESS, NNESS, A(KIXESS), A(KIXNSS),
     &      A(KLTESS), A(KLTNSS), A(KFACSS), SSNAME, MAPEL, MAPND)

      ELSE IF (LISTYP .EQ. 'QA') THEN
         CALL FFADDC (LISTYP, INLINE(1))

         CALL PRQA ('*', NOUT, NQAREC, QAREC, NINFO, INFREC)

C *** EXODUS Print Commands ***

      ELSE IF (LISTYP .EQ. 'NAMES') THEN
         CALL FFADDC (LISTYP, INLINE(1))
         CALL CKEXOD (EXODUS, *150)

         CALL DBVIX_BL ('H', 1, IXHV)
         CALL DBVIX_BL ('G', 1, IXGV)
         CALL DBVIX_BL ('N', 1, IXNV)
         CALL DBVIX_BL ('E', 1, IXEV)
         CALL DBVIX_BL ('M', 1, IXNS)
         CALL DBVIX_BL ('S', 1, IXSS)
         CALL PRNAME (NOUT, NAMLEN,
     *     NVARGL, NVARNP, NVAREL, NVARNS, NVARSS,
     &     NAMES(IXGV), NAMES(IXNV), NAMES(IXEV),
     *     NAMES(IXNS), NAMES(IXSS))

      ELSE IF ((LISTYP .EQ. 'HVARS') .OR. (LISTYP .EQ. 'HISTORY')) THEN
         CALL FFADDC (LISTYP, INLINE(1))
         CALL CKEXOD (EXODUS, *150)
         CALL CKNONE (NVARHI, .FALSE., 'history variables', *150)
         CALL CKNONE (LISHV(0), .TRUE., 'history variables', *150)

         CALL DBVIX_BL ('H', 1, IXHV)

         CALL MDRSRV ('SCRVAR', KVAR, NVARHI)
         CALL MDSTAT (NERR, MEM)
         IF (NERR .GT. 0) GOTO 150

         CALL GETVAR (A, IXHV, -999, NCSTEP, NVARHI, A(KVAR))

         CALL PRSTEP ('*', NOUT,
     &      TIMES(NCSTEP), WHOTIM(NCSTEP), NCSTEP, NSTEPS)

         CALL PRHIST ('*', NOUT, NVARHI, LISHV, NAMES(IXHV), A(KVAR))

         CALL MDDEL ('SCRVAR')

      ELSE IF ((LISTYP .EQ. 'GVARS') .OR. (LISTYP .EQ. 'GLOBALS')) THEN
         CALL FFADDC (LISTYP, INLINE(1))
         CALL CKEXOD (EXODUS, *150)
         CALL CKNONE (NVARGL, .FALSE., 'global variables', *150)
         CALL CKWHOL (WHOTIM(NCSTEP), *150)
         CALL CKNONE (LISGV(0), .TRUE., 'global variables', *150)

         CALL DBVIX_BL ('G', 1, IXGV)

         CALL MDRSRV ('SCRVAR', KVAR, NVARGL)
         CALL MDSTAT (NERR, MEM)
         IF (NERR .GT. 0) GOTO 150

         CALL GETVAR (A, IXGV, -999, NCSTEP, NVARGL, A(KVAR))

         CALL PRSTEP ('*', NOUT,
     &      TIMES(NCSTEP), WHOTIM(NCSTEP), NCSTEP, NSTEPS)

         CALL PRGLOB ('*', NOUT, NVARGL, LISGV, NAMES(IXGV), A(KVAR))

         CALL MDDEL ('SCRVAR')

      ELSE IF ((LISTYP .EQ. 'NVARS') .OR. (LISTYP .EQ. 'NODALS')) THEN
         CALL FFADDC (LISTYP, INLINE(1))
         CALL CKEXOD (EXODUS, *150)
         CALL CKNONE (NUMNP, .FALSE., 'nodes', *150)
         CALL CKNONE (NVARNP, .FALSE., 'nodal variables', *150)
         CALL CKWHOL (WHOTIM(NCSTEP), *150)
         CALL CKNONE (LISNP(0), .TRUE., 'nodes', *150)
         CALL CKNONE (LISNV(0), .TRUE., 'nodal variables', *150)

         CALL DBVIX_BL ('N', 1, IXNV)

         CALL MDRSRV ('SCRVAR', KVAR, NUMNP * LISNV(0))
         CALL MDSTAT (NERR, MEM)
         IF (NERR .GT. 0) GOTO 150

         DO 120 I = 1, NVARNP
            IX = LOCINT (I, LISNV(0), LISNV(1))
            IF (IX .GT. 0) THEN
               IXV = NUMNP * (IX-1)
               CALL GETVAR (A, IXNV+I-1, -999, NCSTEP, NUMNP,
     &            A(KVAR+IXV))
            END IF
  120    CONTINUE

         CALL PRSTEP ('*', NOUT,
     &      TIMES(NCSTEP), WHOTIM(NCSTEP), NCSTEP, NSTEPS)

         CALL PRNODE ('*', NOUT, NUMNP, LISNP, NVARNP, LISNV,
     &      NAMES(IXNV), A(KVAR), MAPND)

         CALL MDDEL ('SCRVAR')

      ELSE IF ((LISTYP .EQ. 'EVARS') .OR. (LISTYP .EQ. 'ELEMENTS')) THEN
         CALL FFADDC (LISTYP, INLINE(1))
         CALL CKEXOD (EXODUS, *150)
         CALL CKNONE (NUMEL, .FALSE., 'elements', *150)
         CALL CKNONE (NVAREL, .FALSE., 'element variables', *150)
         CALL CKWHOL (WHOTIM(NCSTEP), *150)
         CALL CKNONE (LISEL(0), .TRUE., 'elements', *150)
         CALL CKNONE (LISEV(0), .TRUE., 'element variables', *150)

         CALL DBVIX_BL ('E', 1, IXEV)

         CALL MDRSRV ('SCRVAR', KVAR, NUMEL * LISEV(0))
         CALL MDSTAT (NERR, MEM)
         IF (NERR .GT. 0) GOTO 150

C      --Transfer needed variables to random disk
         DO 130 I = 1, NVAREL
            IX = LOCINT (I, LISEV(0), LISEV(1))
            IF (IX .GT. 0) THEN
               CALL GETVAR (A, IXEV+I-1, -1, -NCSTEP, NUMEL, A(KVAR))
            END IF
  130    CONTINUE

         DO 140 I = 1, NVAREL
            IX = LOCINT (I, LISEV(0), LISEV(1))
            IF (IX .GT. 0) THEN
               IXV = NUMEL * (IX-1)
               CALL GETVAR (A, IXEV+I-1, -1, NCSTEP, NUMEL,
     &            A(KVAR+IXV))
            END IF
  140    CONTINUE

         CALL PRSTEP ('*', NOUT,
     &      TIMES(NCSTEP), WHOTIM(NCSTEP), NCSTEP, NSTEPS)

         CALL PRELEM ('*', NOUT, NELBLK, NUMEL, NLISEL, LISEL, LENE,
     &      NVAREL, LISEV, NAMES(IXEV), ISEVOK, A(KVAR), MAPEL)

         CALL MDDEL ('SCRVAR')

C *** Miscellaneous Commands ***

      ELSE IF (LISTYP .EQ. 'STEPS') THEN
         CALL PRTIMS ('NM', NOUT, .TRUE., .TRUE.,
     &      NSTEPS, TIMES, WHOTIM)

      ELSE IF (LISTYP .EQ. 'TIMES') THEN
         CALL PRTIMS ('NT', NOUT, .TRUE., .TRUE.,
     &      NSTEPS, TIMES, WHOTIM)

      ELSE IF (LISTYP .EQ. 'MINMAX') THEN
         CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)
         CALL FFADDC (WORD, INLINE(1))
         IVAR = LOCSTR (WORD, NVARHI+NVARGL+NVARNP+NVAREL, NAMES)
         IF (IVAR .LE. 0) THEN
            CALL PRTERR ('CMDERR', 'Expected variable name')
            GOTO 150
         END IF

         CALL SCALER (A, A, 2, NAMES(IVAR), IVAR,
     &      .FALSE., IDUM, 0, DUMMIN, DUMMAX, MAPEL, MAPND)
         CALL MDSTAT (NERR, MEM)
         IF (NERR .GT. 0) GOTO 150

      ELSE
         CALL SHOCMD ('LIST/PRINT Options:', CMDTBL)
         GOTO 150
      END IF

      GOTO 160

  150 CONTINUE
      INLINE(1) = ' '

  160 CONTINUE
      IF (VERB .NE. 'PRINT') INLINE(1) = ' '
      RETURN
      END
