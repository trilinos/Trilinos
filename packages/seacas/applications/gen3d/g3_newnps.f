C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE NEWNPS (IDFRO, IDBCK,
     &   IDNPS, NNNPS, NNNP3, IXNNPS, IXNNP3,
     &   LTNNPS, LTNNP3, FACNPS, FACNP3,
     &   IXNP, NRNP)
C=======================================================================

C   --*** NEWNPS *** (GEN3D) Calculate 3D node sets
C   --   Written by Amy Gilkey - revised 05/05/86
C   --
C   --NEWNPS calculates the node set information for the 3D database.
C   --The front and back node sets are not calculated (since they
C   --are easily derived), and they are not included in the length tally.
C   --
C   --Parameters:
C   --   IDFRO - IN - ids for front surface node sets; (0) = length
C   --   IDBCK - IN - ids for back surface node sets; (0) = length
C   --   IDNPS - IN - the 2D node sets ids
C   --   NNNPS - IN - the number of nodes for each 2D set
C   --   NNNP3 - OUT - the number of nodes for each 3D set
C   --   IXNNPS - IN - the index of the first node for each 2D set
C   --   IXNNP3 - OUT - the index of the first node for each 3D set
C   --   LTNNPS - IN - the nodes for all 2D sets
C   --   LTNNP3 - OUT - the nodes for all 3D sets
C   --   FACNPS - IN - the distribution factors for all 2D sets
C   --   FACNP3 - OUT - the distribution factors for all 3D sets
C   --   IXNP - IN - the new index for each node
C   --   NRNP - IN - the number of new nodes generated for each node
C   --
C   --Common Variables:
C   --   Uses NDBOUT of /DBASE/
C   --   Uses NUMNPS, LNPSNL of /DBNUMS/
C   --   Sets LNPSNO of /DBNUM3/
C   --   Uses NNREPL of /PARAMS/

      INCLUDE 'g3_dbtitl.blk'
      INCLUDE 'g3_dbnums.blk'
      INCLUDE 'g3_dbnum3.blk'
      INCLUDE 'g3_params.blk'

      INTEGER IDFRO(0:*)
      INTEGER IDBCK(0:*)
      INTEGER IDNPS(*)
      INTEGER NNNPS(*), NNNP3(*)
      INTEGER IXNNPS(*), IXNNP3(*)
      INTEGER LTNNPS(*), LTNNP3(*)
      REAL FACNPS(*), FACNP3(*)
      INTEGER IXNP(*), NRNP(*)

C   --Node set ids - unchanged

      CONTINUE

C   --Node set nodes and distribution factors - add on nodes for each
C   --plate/slice

      N = 0
      DO 30 NPS = 1, NUMNPS
C      --Index of each node set nodes
         IXNNP3(NPS) = N + 1

         IP0 = IXNNPS(NPS) - 1
         DO 20 I = 1, NNNPS(NPS)
            INP = LTNNPS(IP0+I)
            JNP = IXNP(INP)
            DO 10 NR = 1, NRNP(INP)
               N = N + 1
               LTNNP3(N) = JNP + NR-1
               FACNP3(N) = FACNPS(IP0+I)
   10       CONTINUE
   20    CONTINUE

C      --Number of nodes in each set
         NNNP3(NPS) = N - IXNNP3(NPS) + 1
   30 CONTINUE

C   --Number of nodes in all sets
      LNPSNO = N

      RETURN
      END
