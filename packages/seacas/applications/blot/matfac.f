C Copyright(C) 2009-2017 National Technology & Engineering Solutions of
C Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
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
C     * Neither the name of NTESS nor the names of its
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

C=======================================================================
      INTEGER FUNCTION MATFAC (LINKF1, MAXNPF, NPFS, iel, numnp, IERR)
C=======================================================================

C   --*** MATFAC *** (MESH) Match face with existing faces
C   --   Written by Amy Gilkey - revised 02/19/88
C   --   Revised by John Glick - 10/20/88
C   --              Sam Key, 06/01/85
C   --
C   --MATFAC searches for the given face in the list of previously
C   --identified faces.  If a match is found, a pointer to the
C   --matching face is returned.
C   --
C   --The search is done by looking at the list of the faces which
C   --contain each node (the NPFS array).  The lists of two diagonally
C   --opposite nodes in the given face are searched for a reference
C   --to a common face, which is a "match".  The references to the
C   --common face are deleted for all nodes in the face.
C   --
C   --Parameters:
C   --   LINKF1 - IN - the nodes of the face
C   --   MAXNPF - IN - the maximum length of the NPFS entry
C   --   NPFS - IN/OUT - the list of unmatched faces containing a node;
C   --      (0,i) = the length of the list
C   --   IEL - IN - the element containing the face (for error message)
C   --   IERR - OUT - = 0 if routine executes with no errors
C   --                = 1 if an error was detected.

      include 'minmax.blk'
      INTEGER LINKF1(4)
      INTEGER NPFS(NUMNP, 0:MAXNPF)

      MATFAC = 0
      IERR = 0

      IF ((NPFS(LINKF1(1),0) .GT. 0)
     &   .AND. (NPFS(LINKF1(2),0) .GT. 0)
     &   .AND. (NPFS(LINKF1(3),0) .GT. 0)
     &   .AND. (NPFS(LINKF1(4),0) .GT. 0)) THEN

C      --Check all prior faces using this node, by looking for the
C      --diagonally opposite node of this face

         INF1 = LINKF1(1)
         INF3 = LINKF1(3)
         DO 150 I1 = 1, NPFS(INF1,0)
            IOLDF = NPFS(INF1,i1)
            DO 140 I3 = 1, NPFS(INF3,0)
               IF (IOLDF .EQ. NPFS(INF3,i3)) THEN

C               --Diagonal matches so check other nodes

                  INF2 = LINKF1(2)
                  DO 100 I2 = 1, NPFS(INF2,0)
                     IF (IOLDF .EQ. NPFS(INF2,i2)) GOTO 110
  100             CONTINUE
                  WRITE (*, 10000) IEL, LINKF1
                  IERR = 1
                  GOTO 140
  110             CONTINUE

                  INF4 = LINKF1(4)
                  DO 120 I4 = 1, NPFS(INF4,0)
                     IF (IOLDF .EQ. NPFS(INF4,i4)) GOTO 130
  120             CONTINUE
                  WRITE (*, 10000) IEL, LINKF1
                  IERR = 1
                  GOTO 140
  130             CONTINUE

C               --Matching faces so point to face and zero out
C               --nodal references to face

                  L = NPFS(INF1,0)
                  IF (L .GT. I1) NPFS(INF1,i1) = NPFS(INF1,l)
                  NPFS(INF1,0) = L - 1
                  L = NPFS(INF3,0)
                  IF (L .GT. I3) NPFS(INF3,i3) = NPFS(INF3,l)
                  NPFS(INF3,0) = L - 1
                  L = NPFS(INF2,0)
                  IF (L .GT. I2) NPFS(INF2,i2) = NPFS(INF2,l)
                  NPFS(INF2,0) = L - 1
                  L = NPFS(INF4,0)
                  IF (L .GT. I4) NPFS(INF4,i4) = NPFS(INF4,l)
                  NPFS(INF4,0) = L - 1

                  MATFAC = IOLDF
                  maxnod = max(maxnod, inf1, inf2, inf3, inf4)
                  minnod = min(minnod, inf1, inf2, inf3, inf4)
                  GO TO 160
               END IF
  140       CONTINUE
  150    CONTINUE
      END IF

  160 CONTINUE
      RETURN
10000  FORMAT (' Poss. Contiguity Prob. at Hex ', I7,
     *  ', Nodes ', 4(I7,1X), 10(I7,1X))
      END
