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
      SUBROUTINE RWNMAP (NDBIN, NDBOUT, NUMNP, NUMNPO, IXNODE,
     &                   MAPND, NEWIX)
C=======================================================================
C $Id: rwnmap.f,v 1.1 2009/04/24 22:28:49 gdsjaar Exp $
C   --*** RWNMAP *** (ALGEBRA) Read and write database node number map
C   --
C   --
C   --Parameters:
C   --   NDBIN, NDBOUT - IN - the input and output database file
C   --   NUMNP  - IN - the number of nodes
C   --   NUMNPO - IN - the number of output nodes
C   --   IXNODE - IN - the indices of the output nodes (iff NUMNPO <> NUMNP)
C   --   IOERR  - OUT - input/output error flag
C   --   MAPND - SCRATCH - the node number map
C   --   NEWIX - SCRATCH - size = NUMNP (iff NUMNPO <> NUMNP)

      INTEGER IXNODE(*)
      INTEGER MAPND(*)
      INTEGER NEWIX(*)

      call exgnnm(ndbin, mapnd, ierr)

      IF ((NUMNPO .GT. 0) .AND. (NUMNP .NE. NUMNPO)) THEN
        do ix = 1, numnpo
          newix(ix) = mapnd(ixnode(ix))
        end do
        do ix = 1, numnpo
          mapnd(ix) = newix(ix)
        end do
        
      END IF

      call expnnm(ndbout, mapnd, ierr)

      RETURN
      END
