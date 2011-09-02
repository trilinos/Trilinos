C    Copyright(C) 2008 Sandia Corporation.  Under the terms of Contract
C    DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
C    certain rights in this software
C    
C    Redistribution and use in source and binary forms, with or without
C    modification, are permitted provided that the following conditions are
C    met:
C    
C    * Redistributions of source code must retain the above copyright
C       notice, this list of conditions and the following disclaimer.
C              
C    * Redistributions in binary form must reproduce the above
C      copyright notice, this list of conditions and the following
C      disclaimer in the documentation and/or other materials provided
C      with the distribution.
C                            
C    * Neither the name of Sandia Corporation nor the names of its
C      contributors may be used to endorse or promote products derived
C      from this software without specific prior written permission.
C                                                    
C    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
C    "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
C    LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
C    A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
C    OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
C    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
C    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
C    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
C    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
C    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
C    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
C    

C=======================================================================
      SUBROUTINE ZMNPS (NUMNP, NUMNPO, IXNODE, NUMNPS, LNPSNL,
     &   IDNPS, NNNPS, IXNNPS, NDNPS, IXDNPS, LTNNPS, FACNPS,
     &   NEWIX, IXNPS)
C=======================================================================

C   --*** ZMNPS *** (ALGEBRA) Compress database node sets
C   --   Written by Amy Gilkey - revised 07/13/89
C   --
C   --ZMNPS compresses the node set information by renumbering the
C   --nodes and removing deleted nodes.
C   --
C   --Parameters:
C   --   NUMNP      - IN - the number of nodes
C   --   NUMNPO     - IN - the number of nodes
C   --   IXNODE     - IN - the indices of the output nodes(iff NUMNPO <> NUMNP)
C   --   NUMNPS - IN/OUT - the number of node sets
C   --   LNPSNL - IN/OUT - the length of the node sets node list
C   --   IDNPS  - IN/OUT - the node set ID for each set
C   --   NNNPS  - IN/OUT - the number of nodes for each set
C   --   IXNNPS - IN/OUT - the index of the first node for each set
C   --   NDNPS  - IN/OUT - number of distribution factors for each node set
C   --   IXDNPS - IN/OUT - indices into FACNPS; location of 1st dist fact/n set
C   --   LTNNPS - IN/OUT - the nodes for all sets
C   --   FACNPS - IN/OUT - the distribution factors for all sets
C   --   NEWIX - SCRATCH - size = NUMNP
C   --   IXNPS - SCRATCH - size = LNPSNL

      INTEGER NUMNP
      INTEGER NUMNPO
      INTEGER IXNODE(*)
      INTEGER NUMNPS
      INTEGER LNPSNL
      INTEGER IDNPS(*)
      INTEGER NNNPS(*)
      INTEGER IXNNPS(*)
      INTEGER NDNPS(*)
      INTEGER IXDNPS(*)
      INTEGER LTNNPS(*)
      REAL    FACNPS(*)
      INTEGER NEWIX(*)
      INTEGER IXNPS(*)

      IF (NUMNP .EQ. NUMNPO) RETURN

      DO 100 INP = 1, NUMNP
         NEWIX(INP) = 0
  100 CONTINUE
      DO 110 IX = 1, NUMNPO
         NEWIX(IXNODE(IX)) = IX
  110 CONTINUE

      NLO = 0
      DO 120 NL = 1, LNPSNL
         IF (NEWIX(LTNNPS(NL)) .GT. 0) THEN
            NLO = NLO + 1
            LTNNPS(NLO) = NEWIX(LTNNPS(NL))
            FACNPS(NLO) = FACNPS(NL)
            IXNPS(NL) = NLO
         ELSE
            IXNPS(NL) = 0
         END IF
  120 CONTINUE
      LNPSNL = NLO

      NNPSO = 0
      DO 140 INPS = 1, NUMNPS
         NN = 0
         IX0 = 0
         DO 130 IX = IXNNPS(INPS), IXNNPS(INPS)+NNNPS(INPS)-1
            IF (IXNPS(IX) .GT. 0) THEN
               NN = NN + 1
               IF (IX0 .EQ. 0) IX0 = IXNPS(IX)
            END IF
  130    CONTINUE
         IF (NN .GT. 0) THEN
            NNPSO = NNPSO + 1
            IDNPS(NNPSO)  = IDNPS(INPS)
            NNNPS(NNPSO)  = NN
            IXNNPS(NNPSO) = IX0
            NDNPS(NNPSO)  = NN
            IXDNPS(NNPSO) = IX0
         END IF
  140 CONTINUE
      NUMNPS = NNPSO

      RETURN
      END
