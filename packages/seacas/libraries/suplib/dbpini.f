C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE DBPINI (OPTION, NDB, TITLE, NDIM, NUMNP, NUMEL, NELBLK,
     &   NUMNPS, LNPSNL, NUMESS, LESSEL, LESSNL,
     &   NVARHI, NVARGL, NVARNP, NVAREL)
C=======================================================================

C   --*** DBPINI *** (EXOLIB) Display database title and initial variables
C   --   Written by Amy Gilkey - revised 10/15/87
C   --
C   --DBPINI displays the database filename (optional) and title and
C   --the database initial variables.
C   --
C   --Parameters:
C   --   OPTION - IN - '*' to print all, else print options:
C   --      'N' to print database name
C   --      'T' to print title
C   --      'I' to print number of nodes, etc.
C   --      'S' to print node set and side set information
C   --      'V' to print number of variables
C   --   NDB - IN - the database file (if OPTION)
C   --   TITLE - IN - the database title (if OPTION)
C   --   NDIM - IN - the number of coordinates per node (if OPTION)
C   --   NUMNP - IN - the number of nodes (if OPTION)
C   --   NUMEL - IN - the number of elements (if OPTION)
C   --   NELBLK - IN - the number of element blocks (if OPTION)
C   --   NUMNPS - IN - the number of node sets (if OPTION)
C   --   LNPSNL - IN - the length of the node sets node list (if OPTION)
C   --   NUMESS - IN - the number of side sets (if OPTION)
C   --   LESSEL - IN - the length of the side sets element list
C   --      (if OPTION)
C   --   LESSNL - IN - the length of the side sets node list (if OPTION)
C   --   NVARHI - IN - the number of history variables (if OPTION)
C   --   NVARGL - IN - the number of global variables (if OPTION)
C   --   NVARNP - IN - the number of nodal variables (if OPTION)
C   --   NVAREL - IN - the number of element variables (if OPTION)

      CHARACTER*(*) OPTION
      INTEGER NDB
      CHARACTER*80 TITLE
      INTEGER NDIM, NUMNP, NUMEL, NELBLK,
     &   NUMNPS, LNPSNL, NUMESS, LESSEL, LESSNL

      CHARACTER*132 FILNAM

      IF ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'N') .GT. 0)) THEN
         INQUIRE (NDB, NAME=FILNAM)
         WRITE (*, *)
         WRITE (*, 10040) 'Database:  ', FILNAM(:LENSTR(FILNAM))
      END IF

      IF ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'T') .GT. 0)) THEN
         WRITE (*, *)
         WRITE (*, 10040) TITLE(1:LENSTR(TITLE))
      END IF

      IF ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'I') .GT. 0)) THEN
         WRITE (*, 10000, IOSTAT=IDUM)
     &      NDIM, NUMNP, NUMEL, NELBLK
10000     FORMAT (
     &      /, 1X, 'Number of coordinates per node       =', I12
     &      /, 1X, 'Number of nodes                      =', I12
     &      /, 1X, 'Number of elements                   =', I12
     &      /, 1X, 'Number of element blocks             =', I12
     &      )
      END IF

      IF ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'S') .GT. 0)) THEN
         IF (NUMNPS .GT. 0) THEN
            WRITE (*, 10010, IOSTAT=IDUM)
     &         NUMNPS, LNPSNL
         ELSE
            WRITE (*, 10010, IOSTAT=IDUM)
     &         NUMNPS
         END IF
10010     FORMAT (
     &      /, 1X, 'Number of node sets                  =', I12, :
     &      /, 1X, '   Length of node list               =', I12
     &      )
         IF (NUMESS .GT. 0) THEN
            WRITE (*, 10020, IOSTAT=IDUM)
     &         NUMESS, LESSEL, LESSNL
         ELSE
            WRITE (*, 10020, IOSTAT=IDUM)
     &         NUMESS
         END IF
10020     FORMAT
     &      (  1X, 'Number of side sets                  =', I12 :
     &      /, 1X, '   Length of element list            =', I12
     &      /, 1X, '   Length of node list               =', I12
     &      )
      END IF

      IF ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'V') .GT. 0)) THEN
         WRITE (*, 10030, IOSTAT=IDUM) NVARHI, NVARGL, NVARNP, NVAREL
10030     FORMAT (
     &      /, 1X, 'Number of history variables          =', I12
     &      /, 1X, 'Number of global variables           =', I12
     &      /, 1X, 'Number of variables at each node     =', I12
     &      /, 1X, 'Number of variables at each element  =', I12
     &      )
      END IF

      RETURN
10040  FORMAT (1X, 5A)
      END
