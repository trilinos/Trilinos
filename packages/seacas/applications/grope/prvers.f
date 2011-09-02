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
      SUBROUTINE PRVERS (NDB, NOUT)
C=======================================================================
      include 'exodusII.inc'
      include 'progqa.blk'

      character*1 cdum
      real libver
      
      call exinq(ndb, EXVERS, idum, apiver, cdum, ierr)
      call exinq(ndb, EXDBVR, idum, dbver,  cdum, ierr)
      call exinq(ndb, EXLBVR, idum, libver, cdum, ierr)

      IF (NOUT .GT. 0) WRITE (NOUT, 10000)

      CALL BANNER (NOUT, QAINFO,' ', ' ', ' ')
      IF (NOUT .GT. 0) THEN
         WRITE (NOUT, 10010) apiver, dbver, libver
      ELSE
         WRITE (*, 10010) apiver, dbver, libver
      END IF

10000 format(/, 1x, 'VERSION INFORMATION')
10010 format(/,
     &  1x, 'Database written with Exodus API version: ', F6.3,/,
     &  1x, 'Exodus database version:                  ', F6.3,/,
     &  1x, 'Exodus API Library version (grope using): ', F6.3)
      return
      end
