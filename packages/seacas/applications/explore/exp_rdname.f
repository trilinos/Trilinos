C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details
C=======================================================================
      SUBROUTINE RDNAME (A, C, NDB, KVNAMI, KVNAMO,
     &     IXGV, IXNV, IXEV, IXNS, IXSS, KIEVOK, KNSVOK, KSSVOK)
C=======================================================================

C   --*** RDNAME *** (EXPLORE) Read database names
C   --
C   --RDNAME reads the names of the coordinates, the element block types,
C   --and the database variables.  The element block variable truth table
C   --is also read.  An error message is displayed if the end of file is read.
C   --
C   --Dynamic memory is reserved in this routine.  If there is a problem,
C   --the routine returns normally without printing an error message.
C   --
C   --Parameters:
C   --   A - OUT - the dynamic memory base array
C   --   NDB - IN - the database number
C   --   NDIM - IN - the number of coordinates per node
C   --   NELBLK - IN - the number of element blocks
C   --   MAXELB - IN - the maximum number of element block type names to store
C   --   MAXVAR - IN - the maximum number of variable names to store
C   --   NAMECO - OUT - the names of the coordinates
C   --   EBTYPE - OUT - the names of the element block types
C   --   VNAMEI - OUT - the names of the global, nodal, and element variables
C   --      exactly as input (including blanks and lowercase)
C   --   VNAMEO - OUT - the names of the global, nodal, and element variables
C   --      with blanks deleted and converted to uppercase
C   --   NVARGL - OUT - the number of global variables
C   --   IXGV - OUT - the VNAME index of the global variable names
C   --   NVARNP - OUT - the number of nodal variables
C   --   IXNV - OUT - the VNAME index of the nodal variable names
C   --   NVAREL - OUT - the number of element variables
C   --   IXEV - OUT - the VNAME index of the element variable names
C   --   KIEVOK - OUT - the dynamic memory index of the element block variable
C   --      truth table;
C   --      variable i of block j exists iff ISEVOK(i,j) is NOT 0

      include 'exodusII.inc'
      include 'exp_dbnums.blk'

      DIMENSION A(*)
      CHARACTER*1 C(*)

      CHARACTER*80 ERRMSG

      NVARGL = 0
      NVARNP = 0
      NVAREL = 0
      NVARNS = 0
      NVARSS = 0

      IXGV = 1
      IXNV = 1
      IXEV = 1
      IXNS = 1
      IXSS = 1

      KIEVOK = 1
      KNSVOK = 1
      KSSVOK = 1

C ... Read number of variables of each type
      call exgvp (ndb, 'g', nvargl, ierr)
      if (ierr .ne. 0) go to 210
      call exgvp (ndb, 'n', nvarnp, ierr)
      if (ierr .ne. 0) go to 211
      call exgvp (ndb, 'e', nvarel, ierr)
      if (ierr .ne. 0) go to 212
      call exgvp (ndb, 'm', nvarns, ierr)
      if (ierr .ne. 0) go to 212
      call exgvp (ndb, 's', nvarss, ierr)
      if (ierr .ne. 0) go to 212

C   --Get the name indices
      IXGV  = 1
      IXGVE = IXGV  + NVARGL - 1
      IXNV  = IXGVE + 1
      IXNVE = IXNV  + NVARNP - 1
      IXEV  = IXNVE + 1
      IXEVE = IXEV  + NVAREL - 1
      IXNS  = IXEVE + 1
      IXNSE = IXNS  + NVARNS - 1
      IXSS  = IXNSE + 1
      IXSSE = IXSS  + NVARSS - 1

C ... Allocate space for variable names
      nname = nvargl + nvarnp + nvarel + nvarns + nvarss
      call mcrsrv('VNAMEI', kvnami, nname*namlen)
      call mcrsrv('VNAMEO', kvnamo, nname*namlen)
      call mcstat(nerr, nused)
      if (nerr .gt. 0) go to 240

C ... Read variable names (wrapper used to get character length correct)
      ioff = 0
      call rdnam2(ndb, c(kvnami+ioff), nvargl, 'g', ierr, namlen)
      if (ierr .ne. 0) go to 220
      call fixnam(c(kvnami+ioff), c(kvnamo+ioff), nvargl, namlen)

      ioff = ioff + (nvargl * namlen)
      call rdnam2(ndb, c(kvnami+ioff), nvarnp, 'n', ierr, namlen)
      if (ierr .ne. 0) go to 220
      call fixnam(c(kvnami+ioff), c(kvnamo+ioff), nvarnp, namlen)

      ioff = ioff + (nvarnp * namlen)
      call rdnam2(ndb, c(kvnami+ioff), nvarel, 'e', ierr, namlen)
      if (ierr .ne. 0) go to 220
      call fixnam(c(kvnami+ioff), c(kvnamo+ioff), nvarel, namlen)

      ioff = ioff + (nvarel * namlen)
      call rdnam2(ndb, c(kvnami+ioff), nvarns, 'm', ierr, namlen)
      if (ierr .ne. 0) go to 220
      call fixnam(c(kvnami+ioff), c(kvnamo+ioff), nvarns, namlen)

      ioff = ioff + (nvarns * namlen)
      call rdnam2(ndb, c(kvnami+ioff), nvarss, 's', ierr, namlen)
      if (ierr .ne. 0) go to 220
      call fixnam(c(kvnami+ioff), c(kvnamo+ioff), nvarss, namlen)

C ... Read the element block variable truth table
      CALL MDRSRV ('ISEVOK',  KIEVOK, NELBLK * NVAREL)
      CALL MDRSRV ('ISNSVOK', KNSVOK, NUMNPS * NVARNS)
      CALL MDRSRV ('ISSSVOK', KSSVOK, NUMESS * NVARSS)
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 240

      CALL RDNAM1 (NDB, 'E', NELBLK, NVAREL, A(KIEVOK))
      CALL RDNAM1 (NDB, 'M', NUMNPS, NVARNS, A(KNSVOK))
      CALL RDNAM1 (NDB, 'S', NUMESS, NVARSS, A(KSSVOK))
      RETURN

 210  CONTINUE
      WRITE (ERRMSG, 10000) 'NUMBER OF GLOBAL VARIABLES'
      GOTO 230
 211  CONTINUE
      WRITE (ERRMSG, 10000) 'NUMBER OF NODAL VARIABLES'
      GOTO 230
 212  CONTINUE
      WRITE (ERRMSG, 10000) 'NUMBER OF ELEMENT VARIABLES'
      GOTO 230
 220  CONTINUE
      WRITE (ERRMSG, 10000) 'GLOBAL VARIABLE NAMES'
      GOTO 230
 230  CONTINUE
      CALL WDBERR (IERR, ERRMSG)
 240  CONTINUE
      RETURN

10000 FORMAT (A)
      END

      subroutine rdnam2(ndb, varnam, nvar, vtype, ierr, namlen)
      character*(namlen) varnam(*)
      character*1 vtype

      ierr = 0
C ... Initialize in case there is error
      DO 100 I = 1, nvar
        varnam(i) = '--------------------------------'
  100 CONTINUE

      if (nvar .gt. 0) then
         call exgvan (ndb, vtype, nvar, varnam, ierr)
      end if
      return
      end

      subroutine fixnam(namin, namout, nvar, namlen)
C ... Convert the variable names to uppercase and compress blanks
      character*(namlen) namin(*), namout(*)

      do 100 i = 1, nvar
        namout(i) = namin(i)
        call exupcs(namout(i))
 100  continue
      call pckstr(nvar, namout)

      return
      end
