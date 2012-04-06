C Copyright(C) 2011 Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
C certain rights in this software
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
C * Neither the name of Sandia Corporation nor the names of its
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
C 

C=======================================================================
      SUBROUTINE DBOELB (A, NDB,
     &   IDELB, NUMELB, NUMLNK, NUMATR, LINK, NAMELB, ATRIB, ATRIBNW)
C=======================================================================

C   --*** DBOELB *** (EXOLIB) Write database element blocks
C   --   Written by Amy Gilkey - revised 10/12/87
C   --
C   --DBOELB writes the element block information to the database.
C   --Some dynamic dimensioning is done.
C   --
C   --Parameters:
C   --   NDB - IN - the database file
C   --   NELBS, NELBE - IN - the number of first and last element blocks
C   --      to write
C   --   IDELB - IN - the element block IDs for each block
C   --   NUMELB - IN - the number of elements in each block
C   --   NUMLNK - IN - the number of nodes per element in each block
C   --   NUMATR - IN - the number of attributes in each block
C   --   LINK - IN - the connectivity for each block
C   --   ATRIB - IN - the attributes for each block
C   --   ATRIBNW - IN - the new attributes if block is a 3D beam
      include 'exodusII.inc'

      REAL A(*)
      INTEGER NDB
      INTEGER IDELB
      INTEGER NUMELB
      INTEGER NUMLNK
      INTEGER NUMATR
      INTEGER LINK(*)
      REAL ATRIB(*)
      REAL ATRIBNW(7)
      CHARACTER*(mxstln) NAMELB
      
      IELNK = 0
      IEATR = 0

      call expelb(ndb, IDELB, NAMELB, NUMELB, NUMLNK, NUMATR, IERR)
      call expelc(ndb, idelb, link, ierr)

      if (numatr .gt. 0) then
        if (namelb(:4) .eq. 'BEAM') then
C ... A 3D beam needs special treatment since it has 7 attributes
C     and the input 2D beam will only have 1 or 3 attributes.
C     The attributes have been specified in the 'ATRIBNW' array
C     by the user (or the defaults are used).  Need to expand the
C     single value per block values from the ATRIBNW array into 
C     7 values per element.
          call mdlong('ATRIBNW', KATRIB, 7*numelb)
          CALL MDSTAT (NERR, MEM)
          IF (NERR .GT. 0) THEN
            CALL MEMERR
            STOP 'memory error'
          END IF
          call geneat(ndb, idelb, a(katrib), atribnw, numelb)
        else
          call expeat(ndb, idelb, atrib, ierr)
        end if
      end if

      RETURN
      END

      subroutine geneat(ndb, idelb, atrib, atribnw, numelb)
      real atrib(*)
      real atribnw(7)
      
      i = 0
      do 20 ie = 1, numelb
          do 10 ia=1, 7
            i = i + 1
            atrib(i) = atribnw(ia)
 10     continue
 20   continue
      call expeat(ndb, idelb, atrib, ierr)

      return
      end
      
