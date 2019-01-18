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
      SUBROUTINE ZMFIXD (NELBLK, NUMELB, NUMLNK, LINK, NUMNP, IXNP)
C=======================================================================

C   --*** ZMFIXD *** (GJOIN) Find nodes in element blocks
C   --   Written by Amy Gilkey - revised 01/20/88
C   --
C   --ZMFIXD finds the nodes that are within the existing element blocks.
C   --The new indices of the nodes are stored.
C   --
C   --Parameters:
C   --   NELBLK - IN - the number of element blocks
C   --   NUMELB - IN - the number of elements in each element block
C   --   NUMLNK - IN - the number of nodes per element in each element block
C   --   LINK - IN - the connectivity for all elements
C   --   NUMNP - IN/OUT - the number of nodes; returned compressed
C   --   IXNP - OUT - the index of the compressed node; <= 0 if deleted

      INTEGER NUMELB(*)
      INTEGER NUMLNK(*)
      INTEGER LINK(*)
      INTEGER IXNP(*)

      DO 100 INP = 1, NUMNP
         IXNP(INP) = 0
  100 CONTINUE

      IXL0 = 0
      DO 130 IELB = 1, NELBLK
         DO 120 NE = 1, NUMELB(IELB)
            DO 110 K = 1, NUMLNK(IELB)
               IXNP(LINK(IXL0+K)) = 1
  110       CONTINUE
            IXL0 = IXL0 + NUMLNK(IELB)
  120    CONTINUE
  130 CONTINUE

C   --Index nodes of selected element block within the zoom mesh

      JNP = 0
      DO 140 INP = 1, NUMNP
         IF (IXNP(INP) .GT. 0) THEN
            JNP = JNP + 1
            IXNP(INP) = JNP
         END IF
  140 CONTINUE

      NUMNP = JNP

      RETURN
      END
