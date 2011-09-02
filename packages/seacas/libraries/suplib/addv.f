C Copyright(C) 2009 Sandia Corporation. Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
C certain rights in this software.
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

C $Id: addv.f,v 1.2 2009/03/25 12:46:01 gdsjaar Exp $
C $Log: addv.f,v $
C Revision 1.2  2009/03/25 12:46:01  gdsjaar
C Add copyright and license notice to all files.
C
C Revision 1.1  1993/11/18 21:32:06  gdsjaar
C Added scilib routines saxpy.f scopy.f sdot.f snrm2.f subv.f
C Added utility routines addv.f subv.f
C
C-----------------------------------------------------------------------
      SUBROUTINE ADDV( N,A,B,C )
C
C***********************************************************************
C
C     DESCRIPTION: This routine adds two vectors
C
C     FORMAL PARAMETERS:
C        N        INTEGER   Number of entries in A, B
C        A        REAL      First vector
C        B        REAL      Vector to be added
C        C        REAL      Vector with the result
C
C***********************************************************************
C
      DIMENSION A(N),B(N),C(N)
C
      DO 100 I = 1,N
        C(I) = A(I) + B(I)
  100 CONTINUE
C
      RETURN
      END
