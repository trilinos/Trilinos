C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE DBONAM (NDB,
     &   NDIM, NELBLK, NVARHI, NVARGL, NVARNP, NVAREL,
     &   NAMECO, NAMELB, NAMEHV, NAMEGV, NAMENV, NAMEEV, ISEVOK)
C=======================================================================
C   --*** DBONAM *** (EXOLIB) Write database names
C   --
C   --DBONAM writes the names of the coordinates, the element block types,
C   --and the database variables to the database.  The element block variable
C   --truth table is also written.
C   --
C   --Parameters:
C   --   NDB - IN - the database number
C   --   NDIM - IN - the number of coordinates per node; written only if >= 0
C   --   NELBLK - IN - the number of element blocks; written only if >= 0
C   --   NVARHI - IN - the number of history variables; written only if >= 0
C   --   NVARGL - IN - the number of global variables
C   --   NVARNP - IN - the number of nodal variables
C   --   NVAREL - IN - the number of element variables
C   --   NAMECO - IN - the names of the coordinates
C   --   NAMELB - IN - the names of the element block types
C   --   NAMEHV - IN - the names of the history variables
C   --   NAMEGV - IN - the names of the global variables
C   --   NAMENV - IN - the names of the nodal variables
C   --   NAMEEV - IN - the names of the element variables
C   --   ISEVOK - IN - the element block variable truth table;
C   --      variable i of block j exists iff ISEVOK(i,j)
C   --
C   --Database must be positioned in front of coordinate names upon entry;
C   --upon exit positioned after element block variable truth table.

      INTEGER NDB
      INTEGER NDIM, NELBLK, NVARHI, NVARGL, NVARNP, NVAREL
      CHARACTER*8 NAMECO(*)
      CHARACTER*8 NAMELB(*)
      CHARACTER*8 NAMEHV(*)
      CHARACTER*8 NAMEGV(*)
      CHARACTER*8 NAMENV(*)
      CHARACTER*8 NAMEEV(*)
c      LOGICAL ISEVOK(*)
      integer ISEVOK(*)

C   --Write coordinate names

      IF (NDIM .LT. 0) GOTO 100

      IF (NDIM .GT. 0) THEN
         WRITE (NDB) (NAMECO(I), I=1,NDIM)
      ELSE
         WRITE (NDB) 0
      END IF

C   --Write element block type names

      IF (NELBLK .LT. 0) GOTO 100

      IF (NELBLK .GT. 0) THEN
         WRITE (NDB) (NAMELB(I), I=1,NELBLK)
      ELSE
         WRITE (NDB) 0
      END IF

C   --Write the variable names

      IF (NVARHI .LT. 0) GOTO 100

      WRITE (NDB) NVARHI, NVARGL, NVARNP, NVAREL

      IF (NVARHI + NVARGL + NVARNP + NVAREL .GT. 0) THEN
         WRITE (NDB)
     &      (NAMEHV(I), I=1,NVARHI),
     &      (NAMEGV(I), I=1,NVARGL),
     &      (NAMENV(I), I=1,NVARNP),
     &      (NAMEEV(I), I=1,NVAREL)
      ELSE
         WRITE (NDB) 0
      END IF

C   --Write the element block variable truth table

      CALL DBONM1 (NDB, NELBLK, NVAREL, ISEVOK, ISEVOK, MAX(NVAREL,1))

  100 CONTINUE
      RETURN
      END
