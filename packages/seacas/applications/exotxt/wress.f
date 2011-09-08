C Copyright (c) 2007 Sandia Corporation. Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Governement
C retains certain rights in this software.
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
      SUBROUTINE WRESS (NTXT, NUMESS, LESSEL, LESSNL, LESSDF,
     &                  IDESS, NEESS, NDESS, IXEESS, IXDESS, LTEESS,
     &                  LTSESS, FACESS, NNESS, IXNESS, LTNESS, LTNNN )
C=======================================================================

C   --*** WRESS *** (EXOTXT) Write database nodal point sets
C   --   Written by Amy Gilkey - revised 02/27/86
C   --   Modified for ExodusIIv2 database format 10/12/95
C   --
C   --WRESS writes the element side set information from the database.
C   --
C   --Parameters:
C   --   NTXT   - IN - the text file
C   --   NUMESS - IN - the number of side sets
C   --   LESSEL - IN - the length of the element side sets element list
C   --   LESSNL - IN - the length of the element side sets node list
C   --   LESSDF - IN - the length of the side set dist factors list
C   --   IDESS  - IN - array containing side set IDS
C   --   NEESS  - IN - array containing the number of sides for each sets
C   --   NDESS  - IN - array containing the number of dist 
C                      factors for each set
C   --   IXEESS - IN - array containing the indices into the 
C                      LTEESS array which are the locations of the 1st
C                      element of each set
C   --   IXDESS - IN - array containing the indices into the
C                      FACESS array which are the locations of the 1st
C                      distribution factor for each set.
C   --   LTEESS - IN - array containing the elements for all side
C                      sets. Internal element IDS are used in this list
C   --   LTSESS - IN - array containing the sides for all side sets
C   --   FACESS - IN - aray containing dist factors for all side sets
C   --   NNESS  - IN - the number of nodes for each side set
C   --   IXNESS - IN - index into LTNESS - the 1st node for each side set
C   --   LTNESS - IN - array of nodes for all side sets
C   --   LTNNN  - IN - array of number of nodes for each side in a side set

      INTEGER NTXT
      INTEGER NUMESS
      INTEGER LESSEL
      INTEGER LESSNL
      INTEGER LESSDF
      INTEGER IDESS(*)
      INTEGER NEESS(*)
      INTEGER NDESS(*)
      INTEGER IXEESS(*)
      INTEGER IXDESS(*)
      INTEGER LTEESS(*)
      INTEGER LTSESS(*)
      REAL    FACESS(*)
      INTEGER NNESS(*)
      INTEGER IXNESS(*)
      INTEGER LTNESS(*)
      INTEGER LTNNN(*)

      WRITE (NTXT,'(A, I7)') '! Number of Side Sets', NUMESS
      WRITE (NTXT, 10020) LESSEL, LESSNL, LESSDF,
     &      '! Element list, Node list, Dist Fact lengths'

      DO 100 IESS = 1, NUMESS
        WRITE (NTXT, 10025) IDESS(IESS),
     &    NEESS(IESS), NNESS(IESS), NDESS(IESS),
     &    '! ID, #sides, #nodes, #dist fact'
        
        WRITE (NTXT, '(A)') '! Elements and sides for side sets'
C     Element in side set
      IS = IXEESS(IESS)
      IE = IXEESS(IESS) + NEESS(IESS) - 1
      if (is .le. ie)
     &   WRITE (NTXT, 10000) (LTEESS(NL), LTSESS(NL), NL=IS,IE)

      WRITE (NTXT, '(A)') 
     &   '! Nodes and Distribution Factors for side sets'
C     Nodes and Distribution Factors for each side set
      IS = IXDESS(IESS)
      IE = IXDESS(IESS) + NDESS(IESS) - 1
      if (is .le. ie)
     &   WRITE (NTXT, 10010) (LTNESS(NL), FACESS(NL), NL=IS,IE)
 100  CONTINUE

      RETURN
10000  FORMAT (6(I10,1x,i1,1x))
10010  FORMAT (I10, 1pE16.7)
10015  FORMAT (I10, 22X, A)
10020  FORMAT (3I10, 4X, A)
10025  FORMAT (4I10, 4X, A)
      END
