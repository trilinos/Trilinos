C Copyright(C) 1999-2021 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE INIGEN (A, FIRST,
     &   KXN, KYN, KZN, KMAPEL,
     &   KIDELB, KNELB, KNLNK, KNATR, KLINK, KATRIB,
     &   KIDNS, KNNNS, KIXNNS, KLTNNS, KFACNS,
     &   KIDSS, KNESS, KNDSS, KIXESS, KIXDSS, KLTESS, kltsss,
     &   kltsnc, KFACSS, KNMLB, KNMBK, KNMNS, KNMSS)
C=======================================================================

C   --*** INIGEN *** (GJOIN) Initialize the memory for GENESIS database
C   --   Written by Amy Gilkey - revised 10/14/87
C   --
C   --INIGEN initializes the arrays for the GENESIS database.  If this is
C   --the first time through, the memory is reserved (with a length of 0),
C   --otherwise it is just set to a length of 0.
C   --
C   --Parameters:
C   --   A - IN/OUT - the dynamic memory base array
C   --   FIRST - IN - true iff memory is to be reserved
C   --   KXN, KYN, KZN - OUT - index of XN, YN, ZN; nodal coordinates
C   --   KMAPEL - OUT - index of MAPEL; the element order map
C   --   KIDELB - OUT - index of IDELB; the element block IDs for each block
C   --   KNELB - OUT - index of NUMELB; the number of elements in each block
C   --   KNLNK - OUT - index of NUMLNK; the number of nodes per element
C   --      in each block
C   --   KNATR - OUT - index of NUMATR; the number of attributes in each block
C   --   KLINK - OUT - index of LINK; the connectivity for each block
C   --   KATRIB - OUT - index of ATRIB; the attributes for each block
C   --   KIDNS - OUT - index of IDNPS; the nodal point set ID for each set
C   --   KNNNS - OUT - index of NNNPS; the number of nodes for each set
C   --   KIXNNS - OUT - index of IXNNPS; the index of the first node
C   --      for each set
C   --   KLTNNS - OUT - index of LTNNPS; the nodes for all sets
C   --   KFACNS - OUT - index of FACNPS; the distribution factors for all sets
C   --   KIDSS - OUT - index of IDESS; the element side set ID for each set
C   --   KNESS - OUT - index of NEESS; the number of elements for each set
C   --   KNDSS - OUT - index of NDESS; the number of dist-fact for each set
C   --   KIXESS - OUT - index of IXEESS; the index of the first element
C   --      for each set
C   --   KIXDSS - OUT - index of IXDESS; the index of the first dist-fact for each set
C   --   KLTESS - OUT - index of LTEESS; the elements for all sets
C   --   kltsss - OUT - index of LTSESS; the sides for all sets
C   --   kltsnc - OUT - index of LTSSNC; the index of the first dist factor for each face.
C   --   KFACSS - OUT - index of FACESS; the distribution factors for all sets

      include 'gj_namlen.blk'

      DIMENSION A(*)
      LOGICAL FIRST

      if (FIRST) THEN
         namlen = 0
      end if

      IF (FIRST) THEN
         CALL MDRSRV ('XN', KXN, 0)
         CALL MDRSRV ('YN', KYN, 0)
         CALL MDRSRV ('ZN', KZN, 0)
      ELSE
         CALL MDLONG ('XN', KXN, 0)
         CALL MDLONG ('YN', KYN, 0)
         CALL MDLONG ('ZN', KZN, 0)
      END IF

      IF (FIRST) THEN
         CALL MDRSRV ('MAPEL', KMAPEL, 0)
      ELSE
         CALL MDLONG ('MAPEL', KMAPEL, 0)
      END IF

      IF (FIRST) THEN
         CALL MDRSRV ('IDELB', KIDELB, 0)
         CALL MDRSRV ('NUMELB', KNELB, 0)
         CALL MDRSRV ('NUMLNK', KNLNK, 0)
         CALL MDRSRV ('NUMATR', KNATR, 0)
         CALL MDRSRV ('LINK', KLINK, 0)
         CALL MDRSRV ('ATRIB', KATRIB, 0)
      ELSE
         CALL MDLONG ('IDELB', KIDELB, 0)
         CALL MDLONG ('NUMELB', KNELB, 0)
         CALL MDLONG ('NUMLNK', KNLNK, 0)
         CALL MDLONG ('NUMATR', KNATR, 0)
         CALL MDLONG ('LINK', KLINK, 0)
         CALL MDLONG ('ATRIB', KATRIB, 0)
      END IF

      IF (FIRST) THEN
         CALL MDRSRV ('IDNPS', KIDNS, 0)
         CALL MDRSRV ('NNNPS', KNNNS, 0)
         call mdrsrv ('NDNPS', kndns, 0)   ! Node set df count array
         CALL MDRSRV ('IXNNPS', KIXNNS, 0)
         call mdrsrv ('IXDNPS', kixdns, 0) ! Node set df index array
         CALL MDRSRV ('LTNNPS', KLTNNS, 0)
         CALL MDRSRV ('FACNPS', KFACNS, 0) ! Expanded df list array
         call mdrsrv ('CFACNP', kcfacn, 0) ! Exo II df list array
      ELSE
         CALL MDLONG ('IDNPS', KIDNS, 0)
         CALL MDLONG ('NNNPS', KNNNS, 0)
         call mdlong ('NDNPS', kndns, 0)   ! Node set df count array
         CALL MDLONG ('IXNNPS', KIXNNS, 0)
         call mdlong ('IXDNPS', kixdns, 0) ! Node set df index array
         CALL MDLONG ('LTNNPS', KLTNNS, 0)
         CALL MDLONG ('FACNPS', KFACNS, 0) ! Expanded df list array
         call mdlong ('CFACNP', kcfacn, 0) ! Exo II df list array
      END IF

      IF (FIRST) THEN
         CALL MDRSRV ('IDESS', KIDSS, 0)
         CALL MDRSRV ('NEESS', KNESS, 0)
         call mdrsrv ('NDESS', kndss, 0)   ! number of dist factors array
         CALL MDRSRV ('IXEESS', KIXESS, 0)
         call mdrsrv ('IXDESS', kixdss, 0) ! index into dist factors array
         CALL MDRSRV ('LTEESS', KLTESS, 0)
         call mdrsrv ('LTSESS', kltsss, 0) ! side list
         call mdrsrv ('LTSSNC', kltsnc, 0) ! dist-fact index list
         CALL MDRSRV ('FACESS', KFACSS, 0) ! Expanded df list array
      ELSE
         CALL MDLONG ('IDESS', KIDSS, 0)
         CALL MDLONG ('NEESS', KNESS, 0)
         call mdlong ('NDESS', kndss, 0)   ! number of dist factors array
         CALL MDLONG ('IXEESS', KIXESS, 0)
         call mdlong ('IXDESS', kixdss, 0) ! index into dist factors array
         CALL MDLONG ('LTEESS', KLTESS, 0)
         call mdlong ('LTSESS', kltsss, 0) ! side list
         call mdlong ('LTSSNC', kltsnc, 0) ! dist-face index list
         CALL MDLONG ('FACESS', KFACSS, 0) ! Expanded df list array
      END IF

      IF (FIRST) THEN
         CALL MCRSRV ('NAMELB', KNMLB, 0)
         CALL MCRSRV ('NAMBK', KNMBK, 0)
         CALL MCRSRV ('NAMNS', KNMNS, 0)
         CALL MCRSRV ('NAMSS', KNMSS, 0)
      ELSE
         CALL MCLONG ('NAMELB', KNMLB, 0)
         CALL MCLONG ('NAMBK', KNMBK, 0)
         CALL MCLONG ('NAMNS', KNMNS, 0)
         CALL MCLONG ('NAMSS', KNMSS, 0)
      END IF

      RETURN
      END
