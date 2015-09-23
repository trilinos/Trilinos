C Copyright(C) 2009 Sandia Corporation. Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
C certain rights in this software.
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

C=======================================================================
      SUBROUTINE FNODES (IFACE, LINKE, NODEF)
C=======================================================================

C   --*** FNODES *** (MESH) Get nodes for element's face
C   --   Written by Amy Gilkey - revised 05/23/86
C   --              Sam Key, 03/01/85
C   --
C   --FNODES returns the nodes for a given face of a hexahedral element.
C   --The 4-tuple sequence defining the face is counter-clockwise looking
C   --into the element.
C   --
C   --Parameters:
C   --   IFACE - IN - the face number
C   --   LINKE - IN - the hexahedral element connectivity
C   --   NODEF - OUT - the nodes in the extracted face

      INTEGER LINKE(8), NODEF(4)

      INTEGER KFACE(4,6)
      SAVE KFACE
C      --KFACE(*,i) - the indices of the 4 nodes in face i

      DATA KFACE / 1, 4, 3, 2,  5, 6, 7, 8,
     &   3, 4, 8, 7,  2, 3, 7, 6,  1, 2, 6, 5,  1, 5, 8, 4 /

      DO 100 K = 1, 4
         NODEF(K) = LINKE(KFACE(K,IFACE))
  100 CONTINUE

      RETURN
      END
