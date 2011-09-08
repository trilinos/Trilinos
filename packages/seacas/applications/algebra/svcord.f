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
      SUBROUTINE SVCORD (CORD, MAXNE, VARVAL)
C=======================================================================

C   --*** SVCORD *** (ALGEBRA) Save referenced coordinates in VARVAL
C   --   Written by Amy Gilkey - revised 03/15/88
C   --
C   --SVCORD saves the coordinates needed in the equation evaluation
C   --in the VARVAL array.
C   --
C   --Parameters:
C   --   CORD - IN - the coordinates
C   --   MAXNE - IN - the VARVAL dimension (max of NUMEL and NUMNP)
C   --   VARVAL - OUT - the returned coordinates needed
C   --
C   --Common Variables:
C   --   Uses IDVAR, ISTVAR of /VAR../
C   --   Uses NUMNP of /DBNUMS/
C   --   Uses ICOBEG, ICOEND of /DBXVAR/

      PARAMETER (ICURTM = 1, ILSTTM = 2, IONETM = 3)
      include 'namlen.blk'
      include 'var.blk'
      include 'dbnums.blk'
      include 'dbxvar.blk'

      REAL CORD(NUMNP,NDIM)
      REAL VARVAL(MAXNE,*)

C   --Save any coordinates needed in VARVAL

      DO 100 NVAR = ICOBEG, ICOEND
         ID = IDVAR(NVAR)
         NSTO = ISTVAR(ICURTM,NVAR)
         CALL CPYREA (NUMNP, CORD(1,ID), VARVAL(1,NSTO))
  100 CONTINUE

      RETURN
      END
