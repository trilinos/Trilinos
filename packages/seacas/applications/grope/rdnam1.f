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
      SUBROUTINE RDNAM1 (NDB, TYPE, NBLK, NVAR, ISEVOK)
C=======================================================================

C   --*** RDNAM1 *** (EXOLIB) Internal to RDNAME
C   --
C   --RDNAM1 reads the element block variable truth table.
C   --An error message is displayed if the end of file is read.
C   --
C   --Parameters:
C   --   NDB - IN - the database number
C   --   NBLK - IN - the number of element blocks
C   --   NVAR - IN - the number of element variables
C   --   ISEVOK - OUT - the element block variable truth table;
C   --      variable i of block j exists iff ISEVOK(i,j) is NOT 0

      include 'exodusII.inc'
      
      CHARACTER*1 TYPE
      INTEGER ISEVOK(NVAR,*)

      CHARACTER*80 ERRMSG

      CALL INIINT (NBLK*NVAR, 999, ISEVOK)

      if (nblk .gt. 0 .and. nvar .gt. 0) then
         if (type .eq. 'E') then
            call exgvtt (ndb, nblk, nvar, isevok, ierr)
         else if (type .eq. 'M') then
            call exgnstt (ndb, nblk, nvar, isevok, ierr)
         else if (type .eq. 'S') then
            call exgsstt (ndb, nblk, nvar, isevok, ierr)
         end if
        if (ierr .ne. 0) go to 100
      end if

      RETURN

 100  CONTINUE
      WRITE (ERRMSG, 10000) 'VARIABLE TRUTH TABLE'
      GOTO 110
 110  CONTINUE
      CALL WDBERR (IERR, ERRMSG)
 120  CONTINUE
      RETURN
      
10000 FORMAT (A)
      END
