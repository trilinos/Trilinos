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
      SUBROUTINE RDCORD (NDB, NDIM, NUMNP, CORD, NAMECO, ISEOF, NAMLEN)
C=======================================================================

C   --*** RDCORD *** (GROPE) Read database coordinates
C   --
C   --RDCORD reads the coordinate array from the database.  An error
C   --message is displayed if the end of file is read.
C   --
C   --Parameters:
C   --   NDB - IN - the database file
C   --   NDIM - IN - the number of coordinates per node
C   --   NUMNP - IN - the number of nodes
C   --   CORD - OUT - the coordinates
C   --   ISEOF - IN/OUT - set true if end of file read
C   --
C   --Database must be positioned at start of coordinates upon entry;
C   --upon exit at end of coordinates.

      include 'exodusII.inc'

      REAL CORD(NUMNP,*)
      CHARACTER*(NAMLEN) NAMECO(*)
      LOGICAL ISEOF

      CHARACTER*80 ERRMSG

      CALL INIREA (NDIM*NUMNP, 0.0, CORD)
      CALL INISTR (NDIM, '--------------------------------', NAMECO)

      call exgcor(ndb, cord(1,1), cord(1,2), cord(1,3), ierr)
      if (ierr .ne. 0) go to 100

      call exgcon(ndb, nameco, ierr)
      if (ierr .ne. 0) go to 105

      RETURN

 100  CONTINUE
      WRITE (ERRMSG, 10000) 'COORDINATES'
      GOTO 110
 105  CONTINUE
      WRITE (ERRMSG, 10000) 'COORDINATE NAMES'
      GOTO 110
 110  CONTINUE
      CALL WDBERR (IERR, ERRMSG)
      ISEOF = .TRUE.
 120  CONTINUE
      RETURN
      
10000 FORMAT (A)
      END
      
