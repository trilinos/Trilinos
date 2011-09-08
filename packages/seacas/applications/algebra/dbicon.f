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
      SUBROUTINE DBICON (NDB, NDIM, NAMECO)
C=======================================================================

C   --*** DBICON *** Read and pack coordinate names
C   --   Modified from DBINAM for ExodusIIV2 database format 8/26/95
C   --*** DBINAM *** (EXOLIB) Read database names
C   --   Written by Amy Gilkey - revised 02/08/88
C   --
C   --DBINAM performed a number of different input file read base
C   --on the passed in option argument.  DBINAM was split up
C   --into a number of different subroutins
C   --
C   --Parameters:
C   --   NDB    - IN  - the database number
C   --   NDIM   - IN  - the number of coordinates per node
C   --   NAMECO - OUT - the names of the coordinates; max size = 6 (if OPTION)

C   --Routines Called:
C   --   EXUPCS - (SUPES) Convert to uppercase and blank non-standard
C   --   PCKSTR - (STRLIB) Remove embedded blanks

      include 'namlen.blk'
      PARAMETER (MAXDIM=6)

      INTEGER NDB
      INTEGER NDIM
      CHARACTER*(namlen) NAMECO(*)

C     Read and pack coordinate names
      IF (NDIM .GT. MAXDIM) THEN
         write(6, '(A)')'Too many coordinate names in the database'
         IOERR = 1
         RETURN
      END IF

C     Read the name of the coordinate arrays from the database
      call exgcon(ndb, nameco, ierr)

C     Make upper case & remove blanks
      DO 100 I = 1, MIN(NDIM,MAXDIM)
         CALL EXUPCS (NAMECO(I))
  100 CONTINUE
      CALL PCKSTR (MIN(NDIM,MAXDIM), NAMECO)

      RETURN
      END
