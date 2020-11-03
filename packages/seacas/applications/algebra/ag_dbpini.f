C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE DBPINI (OPTION, NDB, TITLE, NDIM, NUMNP, NUMEL, NELBLK,
     &                   NUMNPS, LNPSNL, LNPSDF, NUMESS, LESSEL,
     &                   LESSDF, NVARGL, NVARNP, NVAREL, FILNAM)
C=======================================================================
C   --*** DBPINI *** (EXOLIB) Display database title and initial variables
C   --   Written by Amy Gilkey - revised 10/15/87
C   -- Revised 8/5/95 to work with EXODUSIIV2 database format
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
C   --   NDB    - IN - the database file (if OPTION)
C   --   TITLE  - IN - the database title (if OPTION)
C   --   NDIM   - IN - the number of coordinates per node (if OPTION)
C   --   NUMNP  - IN - the number of nodes (if OPTION)
C   --   NUMEL  - IN - the number of elements (if OPTION)
C   --   NELBLK - IN - the number of element blocks (if OPTION)
C   --   NUMNPS - IN - the number of node sets (if OPTION)
C   --   LNPSNL - IN - the length of the node sets node list (if OPTION)
C   --   LNPSDF - IN - the length of the node set distribution list
C   --   NUMESS - IN - the number of side sets (if OPTION)
C   --   LESSEL - IN - the length of the side sets element list
C   --                 (if OPTION)
C   --   NVARGL - IN - the number of global variables (if OPTION)
C   --   NVARNP - IN - the number of nodal variables (if OPTION)
C   --   NVAREL - IN - the number of element variables (if OPTION)

      include 'exodusII.inc'

      CHARACTER*(*) FILNAM
      CHARACTER*(*) OPTION
      CHARACTER*(MXLNLN) TITLE
      INTEGER NDB
      INTEGER NDIM, NUMNP, NUMEL, NELBLK, NUMNPS,
     &        LNPSNL, LNPSDF, NUMESS, LESSEL, LESSDF
      LOGICAL ALL
      ALL = (OPTION .EQ. '*')

      IF (ALL .OR. (INDEX (OPTION, 'I') .GT. 0)) THEN
         WRITE(*, 10000)
      ENDIF

      IF (ALL .OR. (INDEX (OPTION, 'N') .GT. 0)) THEN
         WRITE (*, 10010) FILNAM(:LENSTR(FILNAM))
      END IF

      IF (ALL .OR. (INDEX (OPTION, 'T') .GT. 0)) THEN
         WRITE (*, 10020) TITLE(1:LENSTR(TITLE))
      END IF

      IF (ALL .OR. (INDEX (OPTION, 'I') .GT. 0)) THEN
         WRITE (*, 10030, IOSTAT=IDUM)
     &      NDIM, NUMNP, NUMEL, NELBLK
      END IF

      IF (ALL .OR. (INDEX (OPTION, 'S') .GT. 0)) THEN
         If (NUMNPS .LE. 0) THEN
            WRITE (*, 10040, IOSTAT=IDUM) NUMNPS
         ELSE
            WRITE (*, 10040, IOSTAT=IDUM) NUMNPS, LNPSNL, LNPSDF
         END IF
         IF (NUMESS .LE. 0) THEN
            WRITE (*, 10050, IOSTAT=IDUM) NUMESS
         ELSE
            WRITE (*, 10050, IOSTAT=IDUM)
     &             NUMESS, LESSEL, LESSDF
         END IF
      END IF

      IF (ALL .OR. (INDEX (OPTION, 'V') .GT. 0)) THEN
         WRITE (*, 10060, IOSTAT=IDUM) NVARGL, NVARNP, NVAREL
      END IF

10000 FORMAT (/, 1X, 'DATABASE INITIAL VARIABLES')
10010 FORMAT (/, 1X, 'Database:  ', A)
10020 FORMAT (/, 1X, A)
10030 FORMAT (
     &      /, 1X, 'Number of coordinates per node         =', I10
     &      /, 1X, 'Number of nodes                        =', I10
     &      /, 1X, 'Number of elements                     =', I10
     &      /, 1X, 'Number of element blocks               =', I10
     &      )
10040 FORMAT (
     &      /, 1X, 'Number of node sets                    =', I10, :
     &      /, 1X, '   Length of node list                 =', I10
     &      /, 1X, '   Length of distribution factors list =', I10
     &      )
10050 FORMAT
     &      (  1X, 'Number of side sets                    =', I10, :
     &      /, 1X, '   Length of element list              =', I10
     &      /, 1X, '   Length of distribution factors list =', I10
     &      )
10060 FORMAT (
     &      /, 1X, 'Number of global variables             =', I10
     &      /, 1X, 'Number of variables at each node       =', I10
     &      /, 1X, 'Number of variables at each element    =', I10
     &      )

      RETURN
      END
