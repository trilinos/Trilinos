C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE DBOESS (NDB, NUMESS, LESSEL, LESSNL,
     &   IDESS, NEESS, NNESS, IXEESS, IXNESS, LTEESS, LTNESS, FACESS)
C=======================================================================

C   --*** DBOESS *** (EXOLIB) Write database node sets
C   --   Written by Amy Gilkey - revised 02/27/86
C   --
C   --DBOESS writes the side set information to the database.
C   --
C   --Parameters:
C   --   NDB - IN - the database file
C   --   NUMESS - IN - the number of side sets
C   --   LESSEL - IN - the length of the side sets element list
C   --   LESSNL - IN - the length of the side sets node list
C   --   IDESS - IN - the side set ID for each set
C   --   NEESS - IN - the number of elements for each set
C   --   NNESS - IN - the number of nodes for each set
C   --   IXEESS - IN - the index of the first element for each set
C   --   IXNESS - IN - the index of the first node for each set
C   --   LTEESS - IN - the elements for all sets
C   --   LTNESS - IN - the nodes for all sets
C   --   FACESS - IN - the distribution factors for all sets
C   --
C   --Database must be positioned at start of side set information
C   --upon entry; upon exit at end of side set information.

      INTEGER NDB
      INTEGER NUMESS, LESSEL, LESSNL
      INTEGER IDESS(*)
      INTEGER NEESS(*)
      INTEGER NNESS(*)
      INTEGER IXEESS(*)
      INTEGER IXNESS(*)
      INTEGER LTEESS(*)
      INTEGER LTNESS(*)
      REAL FACESS(*)

      IF (NUMESS .GT. 0) THEN
         WRITE (NDB) (IDESS(IESS), IESS=1,NUMESS)
         WRITE (NDB) (NEESS(IESS), IESS=1,NUMESS)
         WRITE (NDB) (NNESS(IESS), IESS=1,NUMESS)
         WRITE (NDB) (IXEESS(IESS), IESS=1,NUMESS)
         WRITE (NDB) (IXNESS(IESS), IESS=1,NUMESS)
      ELSE
         WRITE (NDB) 0
         WRITE (NDB) 0
         WRITE (NDB) 0
         WRITE (NDB) 0
         WRITE (NDB) 0
      END IF
      IF (LESSEL .GT. 0) THEN
         WRITE (NDB) (LTEESS(NL), NL=1,LESSEL)
      ELSE
         WRITE (NDB) 0
      END IF
      IF (LESSNL .GT. 0) THEN
         WRITE (NDB) (LTNESS(NL), NL=1,LESSNL)
         WRITE (NDB) (FACESS(NL), NL=1,LESSNL)
      ELSE
         WRITE (NDB) 0
         WRITE (NDB) 0
      END IF

      RETURN
      END
