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
      SUBROUTINE ZMESS (NUMEL, NUMELO, IXELEM, NUMNP, NUMNPO, IXNODE,
     &           NUMESS, LESSEL, LESSDF, IDESS, NEESS, IXEESS,
     &           LTEESS, LTSESS, NDESS, LTNNSSF, FACESS, 
     &           NEWIX, NEWSD, IXESS)
C=======================================================================

C   --*** ZMESS *** (ALGEBRA) Compress database side sets
C   --   Written by Amy Gilkey - revised 07/13/89
C   --   Modified for ExodusII version 2 database format 10/2/95
C   --
C   --ZMESS compresses the side set information by renumbering the
C   --elements and nodes and removing deleted elements (and their nodes).
C   --
C   --Parameters:
C   --   NUMEL      - IN - the number of elements
C   --   NUMELO     - IN - the number of output elements
C   --   IXELEM     - IN - the indices of the output elements
C   --                     (iff NUMELO <> NUMEL)
C   --   NUMNP      - IN - the number of nodes
C   --   NUMNPO     - IN - the number of nodes
C   --   IXNODE     - IN - the indices of the output nodes
C   --                     (iff NUMNPO <> NUMNP)
C   --   NUMESS - IN/OUT - the number of side sets
C   --   LESSEL - IN/OUT - the length of the side sets element list
C   --   IDESS  - IN/OUT - the side set ID for each set
C   --   NEESS  - IN/OUT - the number of elements for each set
C   --   IXEESS - IN/OUT - the index of the first element for each set
C   --   LTEESS - IN/OUT - the elements for all sets
C   --   LTSESS - IN/OUT - array of sides for all side sets
C   --   FACESS - IN/OUT - the distribution factors for all sets
C   --   NDESS  - IN/OUT - the number of distribution factors in each side set
C   --   LTNNSSF - IN/OUT - the number of nodes for each side/face in the side set
C   --   NEWIX - SCRATCH - size = MAX (NUMEL, NUMNP)
C   --   NEWSD - SCRATCH - size = LESSEL
C   --   IXESS - SCRATCH - size = LESSEL

      INTEGER NUMEL
      INTEGER NUMELO
      INTEGER IXELEM(*)
      INTEGER NUMNP
      INTEGER NUMNPO
      INTEGER IXNODE(*)
      INTEGER NUMESS
      INTEGER LESSEL
      INTEGER IDESS(*)
      INTEGER NEESS(*)
      INTEGER IXEESS(*)
      INTEGER LTEESS(*)
      INTEGER LTSESS(*)
      REAL    FACESS(*)
      INTEGER NDESS(*)
      INTEGER LTNNSSF(*)
      INTEGER NEWIX(*)
      INTEGER NEWSD(*)
      INTEGER IXESS(*)

      IF (NUMEL .EQ. NUMELO) RETURN

C   --Get new element number for each existing element number

      DO 100 IEL = 1, NUMEL
         NEWIX(IEL) = 0
  100 CONTINUE
      DO 110 IX = 1, NUMELO
         NEWIX(IXELEM(IX)) = IX
  110 CONTINUE
      DO 115 IX = 1, LESSEL
         NEWSD(IX) = LTSESS(IX)
 115  CONTINUE

C   --Zero out deleted elements in side sets and renumber elements while
C   --packing element list

      NLO = 0
      NDFO = 0
      NDFI = 0
      DO 120 NL = 1, LESSEL
         IF (NEWIX(LTEESS(NL)) .GT. 0) THEN
            NLO = NLO + 1
            LTEESS(NLO) = NEWIX(LTEESS(NL))
            LTSESS(NLO) = NEWSD(NL)
            IXESS(NL) = NLO
            LTNNSSF(NLO) = LTNNSSF(NL)
            if (lessdf .gt. 0) then
              DO IDF = 1, LTNNSSF(NLO)
                FACESS(NDFO + IDF) = FACESS(NDFI + IDF)
              END DO
              NDFO = NDFO + LTNNSSF(NLO)
            end if
         ELSE
            IXESS(NL) = 0
          END IF
          NDFI = NDFI + LTNNSSF(NL)
  120 CONTINUE
      LESSEL = NLO
      LESSDF = NDFO

C   --Adjust side set pointers, etc

      NESSO = 0
      NLO = 0
      DO 210 IESS = 1, NUMESS
         NE = 0
         ND = 0
         IX0 = 0
         DO 200 IX = IXEESS(IESS), IXEESS(IESS)+NEESS(IESS)-1
            IF (IXESS(IX) .GT. 0) THEN
               NE = NE + 1
               NLO = NLO + 1
               ND = ND + LTNNSSF(NLO)
               IF (IX0 .EQ. 0) IX0 = IXESS(IX)
            END IF
  200    CONTINUE
         if (lessdf .eq. 0) nd = 0
         IF (NE .GT. 0) THEN
            NESSO = NESSO + 1
            IDESS(NESSO) = IDESS(IESS)
            NEESS(NESSO) = NE
            NDESS(NESSO) = ND
            IXEESS(NESSO) = IX0
         END IF
  210 CONTINUE
      NUMESS = NESSO

      RETURN
      END
