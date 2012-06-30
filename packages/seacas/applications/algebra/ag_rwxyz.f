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
      SUBROUTINE RWXYZ (NDBIN, NDBOUT, NDIM, NUMNP, NUMNPO,
     &                  IXNODE, CORD, CRDSCR)
C=======================================================================

C   --*** RWXYZ *** (ALGEBRA) Read and write database coordinates
C   --   Written by Amy Gilkey - revised 11/30/87
C   --   Modified for EXODUSIIV2 format 8/29/95
C   --
C   --RWXYZ reads and writes the coordinate array to the database.
C   --Deleted nodes are removed.
C   --
C   --Parameters:
C   --   NDBIN, NDBOUT - IN - the input and output database file
C   --   NDIM   - IN - the number of coordinates per node
C   --   NUMNP  - IN - the number of nodes
C   --   NUMNPO - IN - the number of nodes
C   --   IXNODE - IN - the indices of the output nodes (iff NUMNPO <> NUMNP)
C   --   CORD   - SCRATCH - coordinate I/O
C   --   CRDSCR - SCRATCH - coordinate I/O

      INTEGER NDBIN, NDBOUT
      INTEGER NDIM
      INTEGER NUMNP
      INTEGER IXNODE(*)
      REAL    CORD(NUMNP,NDIM)
      REAL    CRDSCR(NUMNPO,NDIM)


      if (ndim .eq. 2) then
         CALL EXGCOR(ndbin, cord(1,1), cord(1,2), rdum, ierr)
      else if (ndim .eq. 3) then
         CALL EXGCOR(ndbin, cord(1,1), cord(1,2),
     &               cord(1,3), ierr)
      else
         call prterr('FATAL', 'Illegal model dimension')
         RETURN
      end if

      IF ((NUMNPO .GT. 0) .AND. (NDIM .GT. 0)) THEN
         IF (NUMNP .EQ. NUMNPO) THEN
           if (ndim .eq. 2) then
              CALL EXPCOR(ndbout, cord(1,1), cord(1,2), rdum, ierr)
           else if (ndim .eq. 3) then
              CALL EXPCOR(ndbout, cord(1,1), cord(1,2),
     &                    cord(1,3), ioerr)
           else
              call prterr('FATAL', 'Illegal model dimension')
              RETURN
           end if
         ELSE
           do 20 idim=1, ndim
              do 10 ix=1, numnpo
                 crdscr(ix,idim) = cord(ixnode(ix),idim)
 10           continue
 20        continue
           if (ndim .eq. 2) then
              CALL EXPCOR(ndbout, crdscr(1,1), crdscr(1,2),
     &                    rdum, ioerr)
           else if (ndim .eq. 3) then
              CALL EXPCOR(ndbout, crdscr(1,1), crdscr(1,2),
     &                    crdscr(1,3), ioerr)
           else
              call prterr('FATAL', 'Illegal model dimension')
           end if
         END IF
      ELSE
         continue
      END IF

      RETURN
      END
