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
      SUBROUTINE SELMAP (NOUT, TYPE, IDGLO, NUMEL, MAPEL)
C=======================================================================

C   --SELMAP locates the local node corresponding to the global id
C   --
C   --Parameters:
C   --   NOUT - IN - the output file, <=0 for standard
C   --   TYPE - IN - 'node' or 'element'
C   --   NUMEL - IN - the number of node/elements
C   --   MAPEL - IN - the node/element number map

C ... TYPE = Node or Element
      CHARACTER*(*) TYPE
      INTEGER MAPEL(*)

      CHARACTER*128 STRA

      LOGICAL FOUND

      found = .FALSE.
      do 80 i = 1, numel
        if (idglo .eq. mapel(i)) then
          found = .TRUE.
          IF (NOUT .GT. 0) THEN
            write (nout, 10030) TYPE, IDGLO, TYPE, I
          ELSE
            write (*, 10030) TYPE, IDGLO, TYPE, I
          END IF
          go to 90
        end if
 80   continue
 90   continue
      if (.not. found) then
        write (stra, 10000) type, idglo
        call sqzstr(stra, lstra)
        call prterr('WARNING', stra(:lstra))
      end if
      RETURN

10000  FORMAT (' No local ',A,' has global id equal to ',I10)
10030  format (1x, 3x, 'Global ',A,I10,' is local ',A,I10)
      END


