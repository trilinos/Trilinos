C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE WRNPS (NTXT, NUMNPS, LNPSNL, LNPSDF, IDNPS, NNNPS,
     &                  NDNPS, IXNNPS, IXDNPS, LTNNPS, FACNPS)
C=======================================================================

C   --*** WRNPS *** (EXOTXT) Write database nodal point sets
C   --   Written by Amy Gilkey - revised 02/27/86
C   --   Modified for ExodusIIv2 database format 10/12/95
C   --
C   --WRNPS writes the nodal point set information from the database.
C   --
C   --Parameters:
C   --   LNPSNL - IN - the length of the nodal point sets node list
C   --   IDNPS - IN - the nodal point set ID for each set
C   --   NNNPS - IN - the number of nodes for each set
C   --   IXNNPS - IN - the index of the first node for each set
C   --   LSTNPS - IN - the nodes for all sets
C   --   FACNPS - IN - the distribution factors for all sets

C   --   NTXT   - IN - the text file
C   --   NUMNPS - IN - the number of nodal point sets
C   --   LNPSNL - IN - length of node set node list
C   --   LNPSDF - IN - length of node set distribution factor list
C   --   IDNPS  - IN - array containing the node set ID's for each node set
C   --   NNNPS  - IN - array containing the number of nodes for each node set
C   --   NDNPS  - IN - array containing number of dist. fact for each node set
C   --   IXNNPS - IN - array containing indices into the LTNNPS array which
C                      are the location of the 1st nodes for each set
C   --   IXDNPS - IN - array containing indices into the FACNPS array which
C                      are the location of the 1st dist factor for each set
C   --   LTNNPS - IN - array containing the nodes for all node sets
C                      Internal node IDs
C   --   FACNPS - IN - array containing the distribution factors
C                      for all node sets

      IMPLICIT NONE
      INTEGER NTXT, NUMNPS, LNPSNL, LNPSDF, IDS, OFF
      INTEGER INS, INE, INPS, NL
      INTEGER IDNPS(*)
      INTEGER NNNPS(*)
      INTEGER NDNPS(*)
      INTEGER IXNNPS(*)
      INTEGER IXDNPS(*)
      INTEGER LTNNPS(*)
      REAL    FACNPS(*)

      WRITE (NTXT, 10005) '! Node sets                  ', NUMNPS
      WRITE (NTXT, 10005) '! Len node set node list     ', LNPSNL
      WRITE (NTXT, 10005) '! Len node set dist fact list', LNPSDF
      DO 100 INPS = 1, NUMNPS
         WRITE (NTXT, '(A, I10)') '! Nodal point set', INPS
         WRITE (NTXT, 10000) IDNPS(INPS), NNNPS(INPS), NDNPS(INPS),
     &      '! ID, nodes, dist factors'

         INS = IXNNPS(INPS)
         INE = IXNNPS(INPS) + NNNPS(INPS) - 1

         if (ndnps(inps) .ne. 0) then
            IDS = IXDNPS(INPS)
            OFF = IDS - INS
            if (NNNPS(INPS) .NE. NDNPS(INPS)) then
               call prterr('FATAL',
     $             'Number of df in nodeset does not match node count')
            end if
            if (ins .le. ine)
     &        WRITE (NTXT, 10020)
     $           (LTNNPS(NL), FACNPS(NL+OFF), NL=INS,INE)
         else
            if (ins .le. ine)
     &           WRITE (NTXT, 10010) (LTNNPS(NL), NL=INS,INE)
         endif

  100 CONTINUE

      RETURN
10000 FORMAT (3I10, 6X, A)
10005 FORMAT (A, I10)
10010 FORMAT (I10)
10020 FORMAT (I10, 2X, 1pE16.7)
      END
