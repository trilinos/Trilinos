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
C   --   FACNPS - IN - array containing the distribtion factors
C                      for all node sets

      INTEGER NTXT, NUMNPS, LNPSNL, LNPSDF
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
         if (ins .le. ine)
     &      WRITE (NTXT, 10020) (LTNNPS(NL), FACNPS(NL), NL=INS,INE)

  100 CONTINUE

      RETURN
10000 FORMAT (3I10, 6X, A)
10005 FORMAT (A, I10)
10020 FORMAT (I10, 2X, 1pE16.7)
      END
