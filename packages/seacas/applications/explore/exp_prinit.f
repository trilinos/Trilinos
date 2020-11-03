C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details
C=======================================================================
      SUBROUTINE PRINIT (OPTION, NOUT, DBNAME, TITLE,
     &   NDIM, NUMNP, NUMEL, NELBLK,
     &   NUMNPS, LNPSNL, lnpsdf, NUMESS, LESSEL, LESSNL, LESSDF,
     &   NVARGL, NVARNP, NVAREL, NVARNS, NVARSS)
C=======================================================================

C   --*** PRINIT *** (EXPLORE) Display database initial variables
C   --
C   --PRINIT displays the database initial variables.
C   --
C   --Parameters:
C   --   OPTION - IN - '*' to print all, else print options:
C   --      'N' to print database name
C   --      'T' to print title
C   --      'I' to print number of nodes, etc.
C   --      'S' to print nodal point and element side set information
C   --      'V' to print number of variables
C   --      'C' to print number of coordinate frames

C   --   NOUT - IN - the output file, <=0 for standard
C   --   DBNAME - IN - the database name
C   --   TITLE - IN - the database title
C   --   NUMNP - IN - the number of nodes
C   --   NDIM - IN - the number of coordinates per node
C   --   NUMEL - IN - the number of elements
C   --   NELBLK - IN - the number of element blocks
C   --   NUMNPS - IN - the number of nodal points sets
C   --   LNPSNL - IN - the length of the nodal point sets node list
C   --   NUMESS - IN - the number of side sets
C   --   LESSEL - IN - the length of the element side sets element list
C   --   LESSNL - IN - the length of the element side sets node list
C   --   LESSDF - IN - the length of the side sets distribution list
C   --   NVARGL - IN - the number of global variables (if OPTION)
C   --   NVARNP - IN - the number of nodal variables (if OPTION)
C   --   NVAREL - IN - the number of element variables (if OPTION)
C   --   NVARNS - IN - the number of nodeset variables (if OPTION)
C   --   NVARSS - IN - the number of sideset variables (if OPTION)

      include 'exp_dbase.blk'
      include 'exodusII.inc'

      CHARACTER*(*) OPTION
      CHARACTER*(*) DBNAME
      CHARACTER*80 TITLE

      IF ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'I') .GT. 0)) THEN
         IF (NOUT .GT. 0) WRITE (NOUT, 10000)
      END IF

      IF ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'N') .GT. 0)) THEN
         IF (NOUT .GT. 0) THEN
            WRITE (NOUT, 10010) DBNAME(:LENSTR(DBNAME))
         ELSE
            WRITE (*, 10010) DBNAME(:LENSTR(DBNAME))
         END IF
      END IF

      IF ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'T') .GT. 0)) THEN
         IF (NOUT .GT. 0) THEN
            WRITE (NOUT, 10020) TITLE(1:LENSTR(TITLE))
         ELSE
            WRITE (*, 10020) TITLE(1:LENSTR(TITLE))
         END IF
      END IF

      IF ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'I') .GT. 0)) THEN
         IF (NOUT .GT. 0) THEN
            WRITE (NOUT, 10030, IOSTAT=IDUM)
     &         NDIM, NUMNP, NUMEL, NELBLK
         ELSE
            WRITE (*, 10030, IOSTAT=IDUM)
     &         NDIM, NUMNP, NUMEL, NELBLK
         END IF
      END IF

      IF ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'S') .GT. 0)) THEN
         IF (NUMNPS .LE. 0) THEN
            IF (NOUT .GT. 0) THEN
               WRITE (NOUT, 10040, IOSTAT=IDUM) NUMNPS
            ELSE
               WRITE (*, 10040, IOSTAT=IDUM) NUMNPS
            END IF
         ELSE
            IF (NOUT .GT. 0) THEN
               WRITE (NOUT, 10040, IOSTAT=IDUM) NUMNPS, LNPSNL, LNPSDF
            ELSE
               WRITE (*, 10040, IOSTAT=IDUM) NUMNPS, LNPSNL, LNPSDF
            END IF
         END IF
         IF (NUMESS .LE. 0) THEN
            IF (NOUT .GT. 0) THEN
               WRITE (NOUT, 10050, IOSTAT=IDUM) NUMESS
            ELSE
               WRITE (*, 10050, IOSTAT=IDUM) NUMESS
            END IF
         ELSE
            IF (NOUT .GT. 0) THEN
               WRITE (NOUT, 10050, IOSTAT=IDUM) NUMESS, LESSEL, LESSNL,
     $              LESSDF
            ELSE
               WRITE (*, 10050, IOSTAT=IDUM) NUMESS, LESSEL, LESSNL,
     $              LESSDF
            END IF
         END IF
      END IF

      IF ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'C') .GT. 0)) THEN
        call exinq(ndb, EXNCF, ncf, rdum, cdum, ierr)
        if (INDEX (OPTION, 'C') .GT. 0 .or. ncf .gt. 0) then
          IF (NOUT .GT. 0) THEN
            WRITE (NOUT, 10055, IOSTAT=IDUM) NCF
          ELSE
            WRITE (*, 10055, IOSTAT=IDUM) NCF
          END IF
        END IF
      END IF

      IF ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'V') .GT. 0)) THEN
         IF (NOUT .GT. 0) THEN
            WRITE (NOUT, 10060, IOSTAT=IDUM)
     &         NVARGL, NVARNP, NVAREL, NVARNS, NVARSS
         ELSE
            WRITE (*, 10060, IOSTAT=IDUM)
     &         NVARGL, NVARNP, NVAREL, NVARNS, NVARSS
         END IF
      END IF

      RETURN

10000  FORMAT (/, 1X, 'DATABASE INITIAL VARIABLES')
10010  FORMAT (/, 1X, 'Database:  ', A)
10020  FORMAT (/, 1X, A)
10030  FORMAT (
     &   /, 1X, 'Number of coordinates per node       =', I12
     &   /, 1X, 'Number of nodes                      =', I12
     &   /, 1X, 'Number of elements                   =', I12
     &   /, 1X, 'Number of element blocks             =', I12
     &   )
10040  FORMAT (
     &   /, 1X, 'Number of nodal point sets           =', I12, :
     &   /, 1X, '   Length of node list               =', I12
     &   /, 1X, '   Length of distribution list       =', I12
     &   )
10050  FORMAT
     &   (  1X, 'Number of element side sets          =', I12, :
     &   /, 1X, '   Length of element list            =', I12
     &   /, 1X, '   Length of node list               =', I12
     &   /, 1X, '   Length of distribution list       =', I12
     &   )
10055  FORMAT (
     &   /, 1X, 'Number of coordinate frames          =', I12
     &   )

10060  FORMAT (
     &   /, 1X, 'Number of global variables           =', I12
     &   /, 1X, 'Number of variables at each node     =', I12
     &   /, 1X, 'Number of variables at each element  =', I12
     &   /, 1X, 'Number of variables at each nodeset  =', I12
     &   /, 1X, 'Number of variables at each sideset  =', I12
     &   )
      END
