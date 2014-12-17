C Copyright (C) 2009 Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
C certain rights in this software
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
      SUBROUTINE DBIST1 (NDB, NVARNP, NUMNP, VARNP, IVAR, *)
C=======================================================================
C   --*** DBIST1 *** (EXOLIB) Internal to DBISTE, Read nodal variables 
C   --
C   --DBIST1 reads the database nodal variables for one time step.
C   --
C   --Parameters:
C   --   NDB - IN - the database number
C   --   NVARNP - IN - the number of nodal variables
C   --   NUMNP - IN - the number of nodes
C   --   VARNP - OUT - the nodal variables for the time step (if OPTION)
C   --   IVAR  - OUT - the nodal variable number if error on read.
C   --   * - return statement if error encountered, including end-of-file;
C   --      message is printed
C   --
      INTEGER NDB
      INTEGER NVARNP, NUMNP
      REAL VARNP(*)

      iov = 0
      DO 100 IVAR = 1, NVARNP
         READ (NDB, END=190, ERR=190, IOSTAT=IERR)
     &        (VARNP(iov+INP), INP=1,NUMNP)
         iov = iov+numnp
  100 CONTINUE
      RETURN
  190 CONTINUE
      RETURN 1
      END
