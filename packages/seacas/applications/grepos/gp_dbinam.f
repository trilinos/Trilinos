C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE DBINAM (NDB, OPTION, NDIM, NELBLK, NUMNPS, NUMESS,
     *  NNDIM, NNELB, NVARGL, NVARNP, NVAREL, NVARNS, NVARSS,
     &  IXGV, IXNV, IXEV, IXNSV, IXSSV, A, IA,
     *  KIEVOK, KNSVOK, KSSVOK,
     *  C, KNAMES, EXODUS, *)
C=======================================================================
C   --*** DBINAM *** (EXOLIB) Read database names
C   --   Written by Amy Gilkey - revised 02/08/88
C   --
C   --DBINAM reads the names of the coordinates, the element block types,
C   --and the database variables from the database.  All names are converted
C   --to uppercase and all embedded blanks within a name are removed.
C   --The element block variable truth table is also read.
C   --
C   --Note that the numbers of variables are read in this routine.
C   --
C   --Dynamic memory is reserved in this routine.  If there is a problem,
C   --the routine returns normally without printing an error message.
C   --
C   --Parameters:
C   --   NDB    - IN  - the database number
C   --   OPTION - IN  - ' ' to not store, '*' to store all, else store options:
C   --                  'C' to store coordinate names
C   --                  'B' to store element block names
C   --                  'V' to store variables names
C   --                  'T' to store element block variable truth table
C   --   NDIM   - IN  - the number of coordinates per node
C   --   NELBLK - IN  - the number of element blocks
C   --   NNDIM  - OUT - the number of coordinates per node; <0 if end-of-file
C   --   NNELB  - OUT - the number of element blocks; <0 if end-of-file
C   --   NVARGL - OUT - the number of global variables; <0 if end-of-file
C   --   NVARNP - OUT - the number of nodal variables; <0 if end-of-file
C   --   NVAREL - OUT - the number of element variables; <0 if end-of-file
C   --   NVARNS - OUT - the number of nodeset variables; <0 if end-of-file
C   --   NVARSS - OUT - the number of sideset variables; <0 if end-of-file
C   --   IXGV   - OUT - the VNAMES index of the global var names (if OPTION)
C   --   IXNV   - OUT - the VNAMES index of the nodal var names (if OPTION)
C   --   IXEV   - OUT - the VNAMES index of the element var names (if OPTION)
C   --   IXNSV  - OUT - the VNAMES index of the nodeset var names (if OPTION)
C   --   IXSSV  - OUT - the VNAMES index of the sideset var names (if OPTION)
C   --   A      - OUT - the dynamic memory base array
C   --   KIEVOK - OUT - the dynamic memory index of the element block variable
C   --                  truth table; (if OPTION)
C   --                  variable i of block j exists iff ISEVOK(j,i)
C   --   KNSVOK - OUT - the dynamic memory index of the nodeset variable
C   --                  truth table; (if OPTION)
C   --                  variable i of block j exists iff ISEVOK(j,i)
C   --   KSSVOK - OUT - the dynamic memory index of the sideset variable
C   --                  truth table; (if OPTION)
C   --                  variable i of block j exists iff ISEVOK(j,i)
C   --   C      - OUT - the dynamic memory base array (character)
C   --   KNAMES - OUT - the dynamic memory index of the variable names.
C   --   EXODUS - OUT - false if GENESIS file, true if EXODUS file so far
C   --   *      - OUT - return statement if error encountered
C   --                  NOT used if valid GENESIS file; message is printed
C   --
C   --Routines Called:
C   --   EXUPCS - (SUPES) Convert to uppercase and blank non-standard
C   --   MDRSRV - (SUPES) Reserve dynamic memory
C   --   PCKSTR - (STRLIB) Remove embedded blanks

      include 'exodusII.inc'
      include 'gp_namlen.blk'

      PARAMETER (MAXDIM=6)

      INTEGER NDB
      CHARACTER*(*) OPTION
      INTEGER NDIM, NELBLK
      INTEGER NNDIM, NNELB
      INTEGER NVARGL, NVARNP, NVAREL, NVARNS, NVARSS
      INTEGER IXGV, IXNV, IXEV, IXNSV, IXSSV
      DIMENSION A(*), IA(*)
      INTEGER KIEVOK, KNSVOK, KSSVOK
      CHARACTER*1 C(*)
      INTEGER KNAMES
      LOGICAL EXODUS

      CHARACTER*80 ERRMSG
      EXODUS = .FALSE.
      NNDIM = -999
      NNELB = -999
      NVARGL = -999
      NVARNP = -999
      NVAREL = -999

      NNDIM = NDIM

C   --Read the number of variables

      call exgvp(ndb, 'G', nvargl, ierr)
      call exgvp(ndb, 'N', nvarnp, ierr)
      call exgvp(ndb, 'E', nvarel, ierr)
      call exgvp(ndb, 'M', nvarns, ierr)
      call exgvp(ndb, 'S', nvarss, ierr)
      if ((nvargl + nvarnp + nvarel + nvarns + nvarss) .eq. 0) return

      EXODUS = .TRUE.

      IF ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'V') .GT. 0)) THEN

C --Get the name indices
        IXGV = 1
        IXGVE = IXGV + NVARGL - 1
        IXNV = IXGVE + 1
        IXNVE = IXNV + NVARNP - 1
        IXEV = IXNVE + 1
        IXEVE = IXEV + NVAREL - 1
        IXNSV= IXEVE + 1
        IXNSE = IXNSV + NVARNS - 1
        IXSSV = IXNSE + 1
        IXSSE = IXSSV + NVARSS - 1

C ... Allocate space for variable names
        nname = nvargl + nvarnp + nvarel + nvarns + nvarss
        call mcrsrv('NAMES', knames, nname*maxnam)
        call mcstat(nerr, nused)
        if (nerr .gt. 0) go to 180

C ... Read variable names (wrapper used to get character length correct)
        ioff = 0
        call rdnam2(ndb, c(knames+ioff), nvargl, 'g', ierr)
        if (ierr .ne. 0) go to 180

        ioff = ioff + (nvargl * maxnam)
        call rdnam2(ndb, c(knames+ioff), nvarnp, 'n', ierr)
        if (ierr .ne. 0) go to 180

        ioff = ioff + (nvarnp * maxnam)
        call rdnam2(ndb, c(knames+ioff), nvarel, 'e', ierr)
        if (ierr .ne. 0) go to 180

        ioff = ioff + (nvarel * maxnam)
        call rdnam2(ndb, c(knames+ioff), nvarns, 'm', ierr)
        if (ierr .ne. 0) go to 180

        ioff = ioff + (nvarns * maxnam)
        call rdnam2(ndb, c(knames+ioff), nvarss, 's', ierr)
        if (ierr .ne. 0) go to 180

      END IF
C ... Read the element block variable truth table
      IF ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'T') .GT. 0)) THEN
        kmax = max(nelblk*nvarel, numnps*nvarns, numess*nvarss)

        CALL MDRSRV ('ITMP',    KITMP,  KMAX)
        CALL MDRSRV ('ISEVOK',  KIEVOK, NELBLK * NVAREL)
        CALL MDRSRV ('ISNSVOK', KNSVOK, NUMNPS * NVARNS)
        CALL MDRSRV ('ISSSVOK', KSSVOK, NUMESS * NVARSS)
        CALL MDSTAT (NERR, MEM)
        IF (NERR .GT. 0) GOTO 210

        CALL DBINM1 (NDB, 'E', OPTION, NELBLK, NVAREL, IA(KIEVOK),
     &    IA(KIEVOK), IA(KITMP), IERR, MAX(NELBLK,1), *210)

        CALL DBINM1 (NDB, 'M', OPTION, NUMNPS, NVARNS, IA(KNSVOK),
     &    IA(KNSVOK), IA(KITMP), IERR, MAX(NUMNPS,1), *210)

        CALL DBINM1 (NDB, 'S', OPTION, NUMESS, NVARSS, IA(KSSVOK),
     &    IA(KSSVOK), IA(KITMP), IERR, MAX(NUMESS,1), *210)

        call mddel('ITMP')
      END IF

      RETURN

 180  CONTINUE
      ERRMSG = 'ELEMENT BLOCK NAMES'
      GOTO 220
 210  CONTINUE
      ERRMSG = 'ELEMENT BLOCK VARIABLE TRUTH TABLE'
      GOTO 220
 220  CONTINUE
      CALL DBERR (IERR, ERRMSG)
      RETURN 1
      END

      subroutine rdnam2(ndb, varnam, nvar, vtype, ierr)
      include 'gp_namlen.blk'
      character*(maxnam) varnam(*)
      character*1 vtype

      ierr = 0
C ... Initialize in case there is error
      DO 100 I = 1, nvar
        varnam(i) = '-'
 100  CONTINUE

      if (nvar .gt. 0) then
        call exgvan (ndb, vtype, nvar, varnam, ierr)
      end if
      do 200 i = 1, nvar
        call exupcs(varnam(i))
 200  continue
      call pckstr(nvar, varnam)
      return
      end

