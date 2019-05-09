C Copyright(C) 2011-2017 National Technology & Engineering Solutions of
C Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C Redistribution and use in source and binary forms, with or without
C modification, are permitted provided that the following conditions are
C met:
C
C * Redistributions of source code must retain the above copyright
C    notice, this list of conditions and the following disclaimer.
C
C * Redistributions in binary form must reproduce the above
C   copyright notice, this list of conditions and the following
C   disclaimer in the documentation and/or other materials provided
C   with the distribution.
C
C * Neither the name of NTESS nor the names of its
C   contributors may be used to endorse or promote products derived
C   from this software without specific prior written permission.
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

C=======================================================================
      SUBROUTINE ZNPS (IDNPS, NNNPS, IXNNPS, LTNNPS, FACNPS, NUMNPS,
     &   LNPSNL, MAXVAL)
C=======================================================================
C   --   NUMNPS - IN - the number of node sets
C   --   LNPSNL - IN - the length of the node sets node list
C   --   IDNPS - OUT - the node set ID for each set
C   --   NNNPS - OUT - the number of nodes for each set
C   --   IXNNPS - OUT - the index of the first node for each set
C   --   LTNNPS - OUT - the nodes for all sets
C   --   FACNPS - OUT - the distribution factors for all sets

      INTEGER NUMNPS, LNPSNL
      INTEGER IDNPS(*)
      INTEGER NNNPS(*)
      INTEGER IXNNPS(*)
      INTEGER LTNNPS(*)
      REAL FACNPS(*)

      IOFF = 0
      ISET = 0
      DO 20 I=1, NUMNPS
         IBEG = IXNNPS(I)
         IEND = IBEG + NNNPS(I) - 1
         IF (IDNPS(I) .NE. 0) THEN
C ... This set may have been partially or totally deleted
C ... Count number of deleted nodes (node number = MAXVAL)
            NUMDEL = NUMEQI (MAXVAL, NNNPS(I), LTNNPS(IBEG))
            IF (NUMDEL .EQ. NNNPS(I)) THEN
C ... This set has been totally deleted
               IDNPS(I) = 0
            ELSE
C ... This set has been partially deleted, NNNPS(I)-NUMDEL nodes left
               ISET = ISET + 1
               IDNPS(ISET)  = IDNPS(I)
               NNNPS(ISET)  = NNNPS(I) - NUMDEL
               IXNNPS(ISET) = IOFF + 1
               DO 10 INOD = IBEG, IEND
                  IF (LTNNPS(INOD) .LT. MAXVAL) THEN
                     IOFF = IOFF + 1
                     LTNNPS(IOFF) = LTNNPS(INOD)
                     FACNPS(IOFF) = FACNPS(INOD)
                  END IF
   10          CONTINUE
            END IF
         END IF
         IF (IDNPS(I) .EQ. 0) THEN
C ... This set has been totally deleted
         END IF
   20 CONTINUE

      NUMNPS = ISET
      LNPSNL = IOFF

      RETURN
      END
