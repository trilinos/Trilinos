C Copyright (c) 2008 Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
C certain rights in this software
C 
C Redistribution and use in source and binary forms, with or without
C modification, are permitted provided that the following conditions are
C met:
C 
C     * Redistributions of source code must retain the above copyright
C       notice, this list of conditions and the following disclaimer.
C 
C     * Redistributions in binary form must reproduce the above
C       copyright notice, this list of conditions and the following
C       disclaimer in the documentation and/or other materials provided
C       with the distribution.
C 
C     * Neither the name of Sandia Corporation nor the names of its
C       contributors may be used to endorse or promote products derived
C       from this software without specific prior written permission.
C 
C THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
C "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
C LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
C A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
C OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
C SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
C LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
C DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
C THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
C (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
C OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
C 

C=======================================================================
      SUBROUTINE INIGEN (A, FIRST,
     &   KXN, KYN, KZN, KMAPEL,
     &   KIDELB, KNELB, KNLNK, KNATR, KLINK, KATRIB,
     &   KIDNS, KNNNS, KIXNNS, KLTNNS, KFACNS, 
     &   KIDSS, KNESS, KNDSS, KIXESS, KIXDSS, KLTESS, kltsss,
     &   kltsnc, KFACSS, KNMLB)
C=======================================================================
C $Id: inigen.f,v 1.3 2001/06/26 17:38:54 gdsjaar Exp $

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
C   --   KIXDSS - OUT - index of IXDESS; the indes of the first dist-fact for each set      
C   --   KLTESS - OUT - index of LTEESS; the elements for all sets
C   --   kltsss - OUT - index of LTSESS; the sides for all sets
C   --   kltsnc - OUT - index of LTSSNC; the index of the first dist factor for each face.
C   --   KFACSS - OUT - index of FACESS; the distribution factors for all sets

      DIMENSION A(*)
      LOGICAL FIRST

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
      ELSE
         CALL MDLONG ('NAMELB', KNMLB, 0)
      END IF

      RETURN
      END
