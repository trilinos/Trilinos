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
      SUBROUTINE MATEDG (LENF, IELBST, IEDSET, NEDGES)
C=======================================================================

C   --*** MATEDG *** (MESH) Delete edges of selected element block
C   --   Written by Amy Gilkey - revised 07/06/87
C   --
C   --MATEDG deletes all edges that are of a selected element block.
C   --
C   --Parameters:
C   --   LENF - IN - the cumulative face counts by element block
C   --   IELBST - IN - the element block status (>0 if selected)
C   --   IEDSET - IN/OUT - the edge line set;
C   --      (0) = face defining edge; 0 to delete edge
C   --   NEDGES - IN/OUT - the number of lines in the edge set
C   --
C   --Common Variables:
C   --   Uses NELBLK of /DBNUMS/

      include 'debug.blk'
      include 'dbnums.blk'

      INTEGER LENF(0:NELBLK)
      INTEGER IELBST(NELBLK)
      INTEGER IEDSET(0:2,*)

      nhid = 0
      DO 120 IEDG = 1, NEDGES

         IFAC = IEDSET(0,IEDG)
         IF (IFAC .EQ. 0) GOTO 120

C      --Find the face element block
         DO 100 IELB = 1, NELBLK
            IF (IFAC .LE. LENF(IELB)) GOTO 110
  100    CONTINUE
  110    CONTINUE

C      --Delete edge if face is of a selected element block

         IF (IELBST(IELB) .GT. 0) THEN
            IEDSET(0,IEDG) = 0
            nhid = nhid + 1
         END IF
  120 CONTINUE
      if ((cdebug .eq. 'HIDDEN') .and. (idebug .ge. 1))
     &   write (*, '(1x,a,i5)') 'edges in selected block =', nhid

      RETURN
      END
