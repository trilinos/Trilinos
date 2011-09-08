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
      SUBROUTINE RWMAP (NDBIN, NDBOUT, NUMEL, NUMELO, IXELEM,
     &                  MAPEL, NEWIX)
C=======================================================================
C $Id: rwmap.f,v 1.6 2009/04/24 22:26:53 gdsjaar Exp $
C   --*** RWMAP *** (ALGEBRA) Read and write database element order map
C   --   Written by Amy Gilkey - revised 04/28/88
C   --   Modified for EXODUSIIV2 format 8/29/95
C   --
C   --RWMAP reads and writes the element order map to the database.
C   --Deleted elements are removed.
C   --
C   --Parameters:
C   --   NDBIN, NDBOUT - IN - the input and output database file
C   --   NUMEL  - IN - the number of elements
C   --   NUMELO - IN - the number of output elements
C   --   IXELEM - IN - the indices of the output elements (iff NUMELO <> NUMEL)
C   --   IOERR  - OUT - input/output error flag
C   --   MAPEL - SCRATCH - the element order map
C   --   NEWIX - SCRATCH - size = NUMEL (iff NUMELO <> NUMEL)
C   --
C   --Database must be positioned at start of map upon entry;
C   --upon exit at end of map.

      INTEGER IXELEM(*)
      INTEGER MAPEL(*)
      INTEGER NEWIX(*)

      call exgmap (ndbin, mapel, ierr)

      IF ((NUMELO .GT. 0) .AND. (NUMEL .NE. NUMELO)) THEN
        do ix = 1, numelo
          newix(ix) = mapel(ixelem(ix))
        end do
        do ix = 1, numelo
          mapel(ix) = newix(ix)
        end do
      END IF

      call expmap(ndbout, mapel, ierr)

      call exgenm(ndbin, mapel, ierr)

      IF ((NUMELO .GT. 0) .AND. (NUMEL .NE. NUMELO)) THEN
        do ix = 1, numelo
          newix(ix) = mapel(ixelem(ix))
        end do
        do ix = 1, numelo
          mapel(ix) = newix(ix)
        end do
      END IF

      call expenm(ndbout, mapel, ierr)

      RETURN
      END
