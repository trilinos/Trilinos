C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE DBINAM (NDB, OPTION, NDIM, NELBLK, NNDIM, NNELB,
     &                   NVARHI, NVARGL, NVARNP, NVAREL, NVARNS,
     &                   NVARSS, NAMECO,
     &                   IXHV, IXGV, IXNV, IXEV, IXNS, IXSS,
     &                   A, IA, KIEVOK, C, KNAMES, EXODUS, IDELB,
     &                   ISHEX, KHEXID, NAMLEN, *)
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
C   --This routine calls DBVINI and uses DBVIX_BL to get the variable name
C   --indices.
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
C   --   NVARHI - OUT - the number of history variables; <0 if end-of-file
C   --   NVARGL - OUT - the number of global variables; <0 if end-of-file
C   --   NVARNP - OUT - the number of nodal variables; <0 if end-of-file
C   --   NVAREL - OUT - the number of element variables; <0 if end-of-file
C   --   NAMECO - OUT - the names of the coordinates; max size = 6 (if OPTION)
C   --   IXHV   - OUT - the VNAMES index of the history var names (if OPTION)
C   --   IXGV   - OUT - the VNAMES index of the global var names (if OPTION)
C   --   IXNV   - OUT - the VNAMES index of the nodal var names (if OPTION)
C   --   IXEV   - OUT - the VNAMES index of the element var names (if OPTION)
C   --   IXNS   - OUT - the VNAMES index of the nodeset var names (if OPTION)
C   --   IXSS   - OUT - the VNAMES index of the sideset var names (if OPTION)
C   --   A      - OUT - the dynamic memory base array
C   --   KIEVOK - OUT - the dynamic memory index of the element block variable
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
      PARAMETER (MAXDIM=6)

      INTEGER NDB
      CHARACTER*(*) OPTION
      INTEGER NDIM, NELBLK
      INTEGER NNDIM, NNELB
      CHARACTER*(MXSTLN) NAMECO(*)
      INTEGER NVARHI, NVARGL, NVARNP, NVAREL
      INTEGER IXHV, IXGV, IXNV, IXEV, IXNS, IXSS
      DIMENSION A(*), IA(*)
      INTEGER KIEVOK
      CHARACTER*1 C(*)
      INTEGER KNAMES
      LOGICAL EXODUS
      INTEGER IDELB(*)

      CHARACTER*80 ERRMSG
      EXODUS = .FALSE.
      NNDIM = -999
      NNELB = -999
      NVARHI = -999
      NVARGL = -999
      NVARNP = -999
      NVAREL = -999
      NVARNS = -999
      NVARSS = -999
C   --Read and pack coordinate names

      IF ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'C') .GT. 0)) THEN
         IF (NDIM .GT. MAXDIM) CALL PRTERR ('WARNING',
     &      'Too many coordinate names in the database')

C ... Easier to just hardwire the coordinate names...
         nameco(1) = 'X'
         nameco(2) = 'Y'
         if (ndim .eq. 3) nameco(3) = 'Z'

         DO 100 I = 1, MIN(NDIM,MAXDIM)
           CALL EXUPCS (NAMECO(I))
  100    CONTINUE
         CALL PCKSTR (MIN(NDIM,MAXDIM), NAMECO)
      END IF
      NNDIM = NDIM

C   --Read the number of variables

      nvarhi = 0
      call exgvp(ndb, 'G', nvargl, ierr)
      call exgvp(ndb, 'N', nvarnp, ierr)
      call exgvp(ndb, 'E', nvarel, ierr)
      call exgvp(ndb, 'M', nvarns, ierr)
      call exgvp(ndb, 'S', nvarss, ierr)

C   --Initialize for DBVTYP_BL and DBVIX_BL

      CALL DBVINI_BL (NVARGL, NVARNP, NVAREL, NVARNS, NVARSS)

      if ((nvarhi + nvargl + nvarnp + nvarel + nvarns + nvarss) .eq. 0)
     *  go to 160

      EXODUS = .TRUE.
      IF ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'V') .GT. 0)) THEN

C      --Get the name indices

         CALL MCRSRV('NAMES', KNAMES,
     *    NAMLEN*(NVARHI+NVARGL+NVARNP+NVAREL+NVARNS+NVARSS))
         CALL MCSTAT (NERR, MEM)
         IF (NERR .GT. 0) GOTO 180
         CALL DBVIX_BL ('H', 1, IXHV)
         CALL DBVIX_BL ('H', NVARHI, IXHVE)

         CALL DBVIX_BL ('G', 1, IXGV)
         CALL DBVIX_BL ('G', NVARGL, IXGVE)

         CALL DBVIX_BL ('N', 1, IXNV)
         CALL DBVIX_BL ('N', NVARNP, IXNVE)

         CALL DBVIX_BL ('E', 1, IXEV)
         CALL DBVIX_BL ('E', NVAREL, IXEVE)

         CALL DBVIX_BL ('M', 1, IXNS)
         CALL DBVIX_BL ('M', NVARNS, IXNSE)

         CALL DBVIX_BL ('S', 1, IXSS)
         CALL DBVIX_BL ('S', NVARSS, IXSSE)

C      --Read and pack variable names

         call dbinm2 (ndb, nvargl, nvarnp, nvarel, nvarns, nvarss,
     &                C(KNAMES+NAMLEN*(IXGV-1)),
     &                C(KNAMES+NAMLEN*(IXNV-1)),
     &                C(KNAMES+NAMLEN*(IXEV-1)),
     &                C(KNAMES+NAMLEN*(IXNS-1)),
     &                C(KNAMES+NAMLEN*(IXSS-1)), namlen)
      END IF

C   --Read the element block variable truth table

      IF ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'T') .GT. 0)) THEN
         CALL MDRSRV ('ISEVOK', KIEVOK, NELBLK * NVAREL)
         CALL MDRSRV ('ITMP',   KITMP,  NELBLK * NVAREL)

         CALL MDSTAT (NERR, MEM)
         IF (NERR .GT. 0) GOTO 160

         CALL DBINM1 (NDB, OPTION, NELBLK, NVAREL, IA(KIEVOK),
     &                IA(KIEVOK), IA(KITMP), IERR, MAX(NELBLK,1),
     &                IDELB, ISHEX, IA(KHEXID), A, IA, *210)
         call mddel('ITMP')
      END IF

  160 CONTINUE
      RETURN

  180 CONTINUE
      ERRMSG = 'ELEMENT BLOCK NAMES'
      GOTO 220
  210 CONTINUE
      ERRMSG = 'ELEMENT BLOCK VARIABLE TRUTH TABLE'
      GOTO 220
  220 CONTINUE
      CALL DBERR (IERR, ERRMSG)
      RETURN 1
      END

      subroutine dbinm2 (ndb, nvargl, nvarnp, nvarel, nvarns, nvarss,
     &                   namgv, namnv, namev, namns, namss, namlen)

      include 'params.blk'
      character*(namlen) namgv(*), namnv(*), namev(*), namns(*),
     *  namss(*)

      if (nvargl .gt. 0) then
        call exgvan(ndb, 'G', nvargl, namgv, ierr)
      end if
      if (nvarnp .gt. 0) then
        call exgvan(ndb, 'N', nvarnp, namnv, ierr)
      end if
      if (nvarel .gt. 0) then
        call exgvan(ndb, 'E', nvarel, namev, ierr)
      end if
      if (nvarns .gt. 0) then
        call exgvan(ndb, 'M', nvarns, namns, ierr)
      end if
      if (nvarss .gt. 0) then
        call exgvan(ndb, 'S', nvarss, namss, ierr)
      end if

      DO 130 I = 1, nvargl
        CALL EXUPCS (namgv(i))
 130  CONTINUE
      CALL PCKSTR (NVARGL, namgv)

      DO 140 I = 1, nvarnp
        CALL EXUPCS (namnv(i))
 140  CONTINUE
      CALL PCKSTR (NVARNP, namnv)

      DO 150 I = 1, nvarel
        CALL EXUPCS (namev(i))
 150  CONTINUE
      CALL PCKSTR (NVAREL, namev)

      DO I = 1, nvarns
        CALL EXUPCS (namns(i))
      END DO
      CALL PCKSTR (NVARNS, namns)

      DO I = 1, nvarss
        CALL EXUPCS (namss(i))
      END DO
      CALL PCKSTR (NVARSS, namss)

      return
      end
