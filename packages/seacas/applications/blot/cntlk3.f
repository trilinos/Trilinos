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
      SUBROUTINE CNTLK3 (NELBLK, LENE, NLNKE, LENLNK, NFACES, NAMELB,
     $     NLNKSC)
C=======================================================================

C   --*** CNTLK3 *** (MESH) Return link array length and number faces (3D)
C   --   Written by Amy Gilkey - revised 03/04/88
C   --
C   --CNTLK3 returns the maximum length of the connectivity array for
C   --the faces and the maximum number of faces that may be defined.
C   --
C   --Parameters:
C   --   NELBLK - IN - the number of element blocks
C   --   LENE - IN - the cumulative element counts by element block
C   --   NLNKE - IN - the number of nodes per element
C   --   LENLNK - OUT - the length of the connectivity array
C   --   NELEMS - OUT - the number of elements

      INTEGER LENE(0:*)
      INTEGER NLNKE(*)
      CHARACTER*(*) NAMELB(*)
      INTEGER NLNKSC(*)

      NFACES = 0
      LENLNK = 0
      DO 100 IELB = 1, NELBLK
         NUME = LENE(IELB) - LENE(IELB-1)
         IF (NAMELB(IELB)(:3) .EQ. 'HEX' .and. NLNKE(IELB) .GE. 8) THEN
            NF = 6
            NL = 4
         ELSE IF (NAMELB(IELB)(:3) .EQ. 'TET') THEN
            NF = 4

C ... Even though there are only 3 nodes per face, we fake that there
C     are 4 nodes for consistency with the rest of blot...
            NL = 4
         ELSE IF (NAMELB(IELB)(:3) .EQ. 'WED') THEN
            NF = 5
            NL = 4
         ELSE IF (NAMELB(IELB)(:3) .EQ. 'PYR') THEN
            NF = 5
            NL = 4
         ELSE
            NF = 1
            NL = NLNKE(IELB)
         END IF
         NLNKSC(IELB) = NL
         NFACES = NFACES + NUME * NF
         LENLNK = LENLNK + NUME * NF * NL
  100 CONTINUE

      RETURN
      END
