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
      SUBROUTINE DBIEBI (NDB, OPTION, IELB, NUMELB, NUMLNK, NUMATR,
     &                   LINK, ATRIB, NATRDM, NLNKDM, *)
C=======================================================================

C   --*** DBIEB1 *** (EXOLIB) Read database element block misc.
C   --
C   --DBIEB1 reads the element block connectivity and attribute information
C   --from the database.  An error message is displayed if the end of file
C   --is read.
C   --
C   --Parameters:
C   --   NDB    - IN  - the database file
C   --   OPTION - IN  - ' ' to not store, '*' to store all, else store options:
C   --                  'C' to store connectivity
C   --                  'A' to store attributes
C   --   IELB   - IN  - the element block number
C   --   NUMELB - IN  - the number of elements in the block
C   --   NUMLNK - IN  - the number of nodes per element;
C   --                  negate to not store connectivity
C   --   NUMATR - IN  - the number of attributes;
C   --                  negate to not store attributes
C   --   LINK   - OUT - the element connectivity for this block
C   --   ATRIB  - OUT - the attributes for this block
C   --   NATRDM - IN  - dimension of atrib array
C   --   NLNKDM - IN  - dimension of link array
C   --   *      - OUT - return statement if end of file or read error

      include 'exodusII.inc'
      INTEGER NDB
      CHARACTER*(*) OPTION
      INTEGER NUMELB, NUMLNK, NUMATR
      INTEGER LINK(NLNKDM, *)
      REAL ATRIB(NATRDM,*)

      IF ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'C') .GT. 0)) THEN
        if (numelb .gt. 0 .and. numlnk .gt. 0) then
          call exgelc(ndb, ielb, link, ierr)
        end if
      END IF

      IF ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'A') .GT. 0)) THEN
        if (numatr .gt. 0) then
           call exgeat(ndb, ielb, atrib, ierr)
        end if
      END IF
      
      RETURN

      END
