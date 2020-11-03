C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details
C=======================================================================
      SUBROUTINE DBICON (NDB, NDIM, NAMECO)
C=======================================================================

C   --*** DBICON *** Read and pack coordinate names
C   --   Modified from DBINAM for ExodusIIV2 database format 8/26/95
C   --*** DBINAM *** (EXOLIB) Read database names
C   --   Written by Amy Gilkey - revised 02/08/88
C   --
C   --DBINAM performed a number of different input file read base
C   --on the passed in option argument.  DBINAM was split up
C   --into a number of different subroutins
C   --
C   --Parameters:
C   --   NDB    - IN  - the database number
C   --   NDIM   - IN  - the number of coordinates per node
C   --   NAMECO - OUT - the names of the coordinates; max size = 6 (if OPTION)

C   --Routines Called:
C   --   EXUPCS - (SUPES) Convert to uppercase and blank non-standard
C   --   PCKSTR - (STRLIB) Remove embedded blanks

      include 'ag_namlen.blk'
      PARAMETER (MAXDIM=6)

      INTEGER NDB
      INTEGER NDIM
      CHARACTER*(namlen) NAMECO(*)

C     Read and pack coordinate names
      IF (NDIM .GT. MAXDIM) THEN
         CALL PRTERR ('WARNING',
     &   'Too many coordinate names in the database')
         RETURN
      END IF

C     Read the name of the coordinate arrays from the database
      call exgcon(ndb, nameco, ierr)

C     Make upper case & remove blanks
      DO 100 I = 1, MIN(NDIM,MAXDIM)
         CALL EXUPCS (NAMECO(I))
  100 CONTINUE
      CALL PCKSTR (MIN(NDIM,MAXDIM), NAMECO)

      RETURN
      END
