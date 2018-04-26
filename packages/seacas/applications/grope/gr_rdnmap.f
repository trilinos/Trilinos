C    Copyright(C) 2008-2017 National Technology & Engineering Solutions of
C    Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
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
C    * Neither the name of NTESS nor the names of its
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
      SUBROUTINE RDNMAP (NDB, NUMNP, MAPNO, ISEOF)
C=======================================================================

C   --*** RDNMAP *** (GROPE) Read database node order map
C   --
C   --RDMAP reads the node order map from the database.  An error
C   --message is displayed if the end of file is read.
C   --
C   --Parameters:
C   --   NDB - IN - the database file
C   --   NUMNP - IN - the number of nodes
C   --   MAPNO - OUT - the node order map
C   --   ISEOF - IN/OUT - set true if end of file read

      include 'exodusII.inc'

      INTEGER MAPNO(*)
      LOGICAL ISEOF

      CHARACTER*80 ERRMSG
      CHARACTER*1  cdum

      CALL INIINT (NUMNP, 0, MAPNO)

C ... Don't warn about no map stored in file
      call exopts (0, ierr1)
      call exgnnm(ndb, mapno, ierr)
      call exopts (EXVRBS, ierr1)
      if (ierr .lt. 0) then
        go to 100
      else if (ierr .ne. 0) then
        do 10 i=1, numnp
          mapno(i) = i
 10     continue
      end if

      RETURN

  100 CONTINUE
      WRITE (ERRMSG, 10000) 'NODE ORDER MAP'
      GOTO 110
  110 CONTINUE
      CALL WDBERR (IERR, ERRMSG)
      ISEOF = .TRUE.

      RETURN

10000  FORMAT (A)
      END
