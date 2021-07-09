C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE DBONPS (NDB, NUMNPS, LNPSNL,
     &   IDNPS, NNNPS, IXNNPS, LTNNPS, FACNPS)
C=======================================================================

C   --*** DBONPS *** (EXOLIB) Write database node sets
C   --   Written by Amy Gilkey - revised 02/27/86
C   --
C   --DBONPS writes the node set information to the database.
C   --
C   --Parameters:
C   --   NDB - IN - the database file
C   --   NUMNPS - IN - the number of node sets
C   --   LNPSNL - IN - the length of the node sets node list
C   --   IDNPS - IN - the node set ID for each set
C   --   NNNPS - IN - the number of nodes for each set
C   --   IXNNPS - IN - the index of the first node for each set
C   --   LTNNPS - IN - the nodes for all sets
C   --   FACNPS - IN - the distribution factors for all sets
C   --
C   --Database must be positioned at start of node set information
C   --upon entry; upon exit at end of node set information.

      INTEGER NDB
      INTEGER NUMNPS, LNPSNL
      INTEGER IDNPS(*)
      INTEGER NNNPS(*)
      INTEGER IXNNPS(*)
      INTEGER LTNNPS(*)
      REAL FACNPS(*)

      IF (NUMNPS .GT. 0) THEN
         WRITE (NDB) (IDNPS(INPS), INPS=1,NUMNPS)
         WRITE (NDB) (NNNPS(INPS), INPS=1,NUMNPS)
         WRITE (NDB) (IXNNPS(INPS), INPS=1,NUMNPS)
      ELSE
         WRITE (NDB) 0
         WRITE (NDB) 0
         WRITE (NDB) 0
      END IF
      IF (LNPSNL .GT. 0) THEN
         WRITE (NDB) (LTNNPS(NL), NL=1,LNPSNL)
         WRITE (NDB) (FACNPS(NL), NL=1,LNPSNL)
      ELSE
         WRITE (NDB) 0
         WRITE (NDB) 0
      END IF

      RETURN
      END
