C Copyright (c) 2007 Sandia Corporation. Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Governement
C retains certain rights in this software.
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
C 
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
C 

C=======================================================================
      SUBROUTINE WRMAP (NTXT, OPTION, NUMNP, NUMEL,
     &                  NPMAP, ELMAP, MAPEL)
C=======================================================================

C   --*** WRMAP *** (EXOTXT) Write database node number map,
C   --                       element number map, and/or element order map
C   --   Written by Amy Gilkey - revised 02/27/86
C   --   Modified for ExodusIIv2 - 10/24/95
C   --
C   --Parameters:
C   --   NTXT   - IN - the database text file
C   --   OPTION - IN - '*' write all
C   --                 'N' write node number map
C   --                 'E' write element number map
C   --                 'O' write element order map
C   --   NUMNP  - IN - number of nodes
C   --   NUMEL  - IN - number of elements
C   --   NPMAP  - IN - node number map (if OPTION)
C   --   ELMAP  - IN - element number map (if OPTION)
C   --   MAPEL  - IN - element order map (if OPTION)

      INTEGER NTXT
      CHARACTER*(*) OPTION
      INTEGER NUMNP
      INTEGER NUMEL
      INTEGER NPMAP(*)
      INTEGER ELMAP(*)
      INTEGER MAPEL(*)
      LOGICAL ALL, NOPT, EOPT, OOPT
      LOGICAL INORDR
      
      ALL  = (OPTION .EQ. '*')
      NOPT = (INDEX (OPTION, 'N') .GT. 0)
      EOPT = (INDEX (OPTION, 'E') .GT. 0)
      OOPT = (INDEX (OPTION, 'O') .GT. 0)

C     Write node number map
      IF (ALL .OR. NOPT) THEN
         WRITE (NTXT, 1000) '! Node number map'
         if (inordr(npmap, numnp)) then
           write (ntxt, 1000) 'sequence 1..numnp'
         else
           write (ntxt, 1000) 'explicit map'
           WRITE (NTXT, 1010) (NPMAP(I), I = 1, NUMNP)
         end if
      END IF

C     Write element number map
      IF (ALL .OR. EOPT) THEN
         WRITE (NTXT, 1000) '! Element number map'
         if (inordr(elmap, numel)) then
           write (ntxt, 1000) 'sequence 1..numel'
         else
           write (ntxt, 1000) 'explicit map'
           WRITE (NTXT, 1010) (ELMAP(I), I = 1, NUMEL)
         end if
      END IF

C     Write element order map
      IF (ALL .OR. OOPT) THEN
         WRITE (NTXT, 1000) '! Element order map'
         if (inordr(mapel, numel)) then
           write (ntxt, 1000) 'sequence 1..numel'
         else
           write (ntxt, 1000) 'explicit map'
           WRITE (NTXT, 1010) (MAPEL(I), I = 1, NUMEL)
         end if
      END IF

 1000 FORMAT (A)
 1010 FORMAT (8I10)

      RETURN
      END

C=======================================================================
      logical function inordr(MAP, ISIZE)
C=======================================================================

C ... Determine if the passed in map is simply a sequence from 1..isize
      
      integer map(isize)
      
      inordr = .FALSE.
      do 10 i=1, isize
        if (map(i) .ne. i) then
          inordr = .FALSE.
          return
        end if
 10   continue
      inordr = .TRUE.
      return
      end
