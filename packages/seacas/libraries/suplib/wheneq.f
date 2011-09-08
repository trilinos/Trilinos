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

C $Id: wheneq.f,v 1.2 2009/03/25 12:46:02 gdsjaar Exp $
C $Log: wheneq.f,v $
C Revision 1.2  2009/03/25 12:46:02  gdsjaar
C Add copyright and license notice to all files.
C
C Revision 1.1  1993/07/06 22:01:27  gdsjaar
C Added wheneq function for non-cray systems.
C
c
      subroutine wheneq(n,array,inc,target,index,nval)
C see WHENEQ (3SCI) Math and Scientific Library, Cray Research, Inc 
      DIMENSION ARRAY(*), INDEX(*)
      INA = 1
      NVAL = 0
      IF(INC .LT. 0) INA = (-INC) * (N-1) + 1
      DO 100 I=1,N
          IF( ARRAY(INA) .EQ. TARGET)THEN
            NVAL = NVAL + 1
            INDEX(NVAL) = I
          ENDIF
          INA = INA + INC
 100  CONTINUE
      RETURN
      END

