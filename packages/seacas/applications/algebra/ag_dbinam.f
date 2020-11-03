C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details
C=======================================================================
      SUBROUTINE DBINAM (NDB, C, KNAMES, NVARGL, NVARNP, NVAREL,
     &                   IXGV, IXNV, IXEV, MAXVAR)
C=======================================================================

C   --*** DBINAM *** Read the names of the of global, node, and
C   --               element result variables
C   --   Modified for ExodusIIV2 database 8/26/95
C   --This code originated from the subroutine
C   --*** DBINAM *** (EXOLIB) Read database names
C   --   Written by Amy Gilkey - revised 02/08/88
C   --
C   --DBINAM performed a number of different input file read base
C   --on the passed in option argument.  DBINAM was split up
C   --into a number of different subroutins
C   --
C   --This routine calls DBVINI and uses DBVIX to get the variable name
C   --indices.  This subroutine includes and call DBINM2.
C   --
C   --Parameters:
C   --   NDB    - IN  - the database number
C   --   C      - OUT - the dynamic memory base array (character)
C   --   NVARGL - IN  - the number of global variables; <0 if end-of-file
C   --   NVARNP - IN  - the number of nodal variables; <0 if end-of-file
C   --   NVAREL - IN  - the number of element variables; <0 if end-of-file
C   --   IXGV   - OUT - the VNAMES index of the global variable names
C   --   IXNV   - OUT - the VNAMES index of the nodal variable names
C   --   IXEV   - OUT - the VNAMES index of the element variable names
C   --   MAXVAR - IN  - maximum number of all variables
C   --Routines Called:
C   --   EXUPCS - (SUPES) Convert to uppercase and blank non-standard
C   --   MDRSRV - (SUPES) Reserve dynamic memory
C   --   PCKSTR - (STRLIB) Remove embedded blanks

      include 'ag_namlen.blk'

      INTEGER NDB
      CHARACTER*1 C(*)
      CHARACTER*(namlen) KNAMES(*)
      INTEGER NVARGL, NVARNP, NVAREL
      INTEGER IXGV, IXNV, IXEV

C     There are no global, nodal, or element variables
      IF ((nvargl + nvarnp + nvarel) .eq. 0) RETURN

C     DBVINI Initializes the indices for DBVTYP and DBVIX
C     DBVINI(num_global_var, num_nodal_var, num_elem_var)
      CALL DBVINI (NVARGL, NVARNP, NVAREL)

C     call dbvix(var_index, ret_var_type, ret_var_num)
      CALL DBVIX ('G', 1, IXGV)
      CALL DBVIX ('G', NVARGL, IXGVE)
      IXGVE = MIN(IXGVE, MAXVAR)
      CALL DBVIX ('N', 1, IXNV)
      CALL DBVIX ('N', NVARNP, IXNVE)
      IXNVE = MIN(IXNVE, MAXVAR)
      CALL DBVIX ('E', 1, IXEV)
      CALL DBVIX ('E', NVAREL, IXEVE)
      IXEVE = MIN(IXEVE, MAXVAR)

C      --Read and pack variable names
      call dbinm2 (ndb, nvargl, nvarnp, nvarel, KNAMES(IXGV),
     &             KNAMES(IXNV), KNAMES(IXEV))

      do 10 i=ixgv, ixgve
        call exupcs(knames(i))
 10   continue
      do 20 i=ixnv, ixnve
        call exupcs(knames(i))
 20   continue
      do 30 i=ixev, ixeve
        call exupcs(knames(i))
 30   continue

      RETURN
      END

C     *********************************************************
      SUBROUTINE DBINM2 (ndb, nvargl, nvarnp, nvarel,
     &                   namgv, namnv, namev)
C     *********************************************************
C     NDB    - IN  - File ID
C     NVARGL - IN  - Number of global variables read
C     NVARNP - IN  - Number of nodal variables read
C     NVAREL - IN  - Number of element variables rea
C     NAMGV  - OUT - Returned character array with global variable names
C     NAMNV  - OUT - Returned character array with nodal variable names
C     NAMEV  - OUT - Returned character array with element variable names

      include 'ag_namlen.blk'

      character*(namlen) namgv(*), namnv(*), namev(*)

C     Reads the names of the results variables from the database
C     call exgvan(fileid, vartyp, num_var, array_var_names(num_var), errorid)
      if (nvargl .gt. 0)
     &   call exgvan(ndb, 'G', nvargl, namgv, ierr)

      if (nvarnp .gt. 0)
     &   call exgvan(ndb, 'N', nvarnp, namnv, ierr)

      if (nvarel .gt. 0)
     &   call exgvan(ndb, 'E', nvarel, namev, ierr)

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

      return
      end
