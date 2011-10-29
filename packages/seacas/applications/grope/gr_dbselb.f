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
      SUBROUTINE DBSELB (NELBLK, NUMEL, LENE, INELB, NLISEL, LISEL)
C=======================================================================

C   --*** DBSBEL *** (BLOT) Select elements list of element blocks
C   --
C   --DBSBEL creates the element block selection array and the element
C   --selection array (by block) given a list of selected element blocks.
C   --
C   --Parameters:
C   --   NELBLK - IN - the number of element blocks
C   --   NUMEL - IN - the number of elements
C   --   LENE - IN - the cumulative element counts by element block
C   --   INELB - IN - the indices of the selected element blocks
C   --   NLISEL - IN/OUT - the number of selected elements for each block
C   --   LISEL - IN/OUT - the indices of the selected elements (by block)

      INTEGER LENE(0:NELBLK)
      INTEGER INELB(0:*)
      INTEGER NLISEL(0:NELBLK)
      INTEGER LISEL(0:*)

      NLISEL(0) = 0
      CALL INIINT (NELBLK, 0, NLISEL(1))
      DO 100 IX = 1, INELB(0)
         IELB = INELB(IX)
         NLISEL(IELB) = LENE(IELB) - LENE(IELB-1)
  100 CONTINUE
      NLISEL(0) = INELB(0)

      LISEL(0) = 0
      CALL INIINT (NUMEL, 0, LISEL(1))

      DO 120 IELB = 1, NELBLK
         IF (NLISEL(IELB) .GT. 0) THEN
            LISEL(0) = LISEL(0) + NLISEL(IELB)
            DO 110 IEL = LENE(IELB-1)+1, LENE(IELB)
               LISEL(IEL) = IEL
  110       CONTINUE
         END IF
  120 CONTINUE

      END
