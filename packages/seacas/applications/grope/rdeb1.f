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
      SUBROUTINE RDEB1 (NDB, IDELB, NUMELB, NUMLNK, NUMATR,
     &   LINK, ATRIB, ATRNM, NAMLEN)
C=======================================================================

C   --*** RDEB1 *** (GROPE) Read database element block misc.
C   --
C   --RDEB1 reads the element block connectivity and attribute information
C   --from the database.  An error message is displayed if the end of file
C   --is read.
C   --
C   --Parameters:
C   --   NDB - IN - the database file
C   --   IDELB - IN - the element block number (for errors)
C   --   NUMELB - IN - the number of elements in the block
C   --   NUMLNK - IN - the number of nodes per element
C   --   NUMATR - IN - the number of attributes
C   --   LINK - OUT - the element connectivity for this block
C   --   ATRIB - OUT - the attributes for this block
C   --
      include 'exodusII.inc'

      INTEGER LINK(*)
      REAL ATRIB(*)

      CHARACTER*80 ERRMSG
      CHARACTER*(NAMLEN) ATRNM(*)

      IF (NUMELB .GT. 0 .AND. NUMLNK .GT. 0) THEN
        CALL EXGELC (NDB, IDELB, LINK, IERR)
        IF (IERR .NE. 0) GO TO 100
      END IF
      
      IF (NUMELB .GT. 0 .AND. NUMATR .GT. 0) THEN
        CALL EXGEAT (NDB, IDELB, ATRIB, IERR)
        IF (IERR .NE. 0) GO TO 110

        CALL EXGEAN (NDB, IDELB, NUMATR, ATRNM, IERR)
        IF (IERR .NE. 0) GO TO 115
      END IF
      
      RETURN

  100 CONTINUE
      WRITE (ERRMSG, 10000, IOSTAT=IDUM) 'CONNECTIVITY for block', IDELB
      GOTO 120
  110 CONTINUE
      WRITE (ERRMSG, 10000, IOSTAT=IDUM) 'ATTRIBUTES for block', IDELB
      GOTO 120
 115  CONTINUE
      WRITE (ERRMSG, 10000, IOSTAT=IDUM) 'ATTRIBUTE NAMES for block', IDELB
      GOTO 120
  120 CONTINUE
      CALL WDBERR (IERR, ERRMSG)
  130 CONTINUE
      RETURN

10000  FORMAT (5 (A, I10))
      END
