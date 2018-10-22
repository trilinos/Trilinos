C Copyright (C) 2009-2017 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
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
C     * Neither the name of NTESS nor the names of its
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

C $Id: ldrota.f,v 1.1 1993/07/16 16:46:34 gdsjaar Exp $ 
C $Log: ldrota.f,v $
C Revision 1.1  1993/07/16 16:46:34  gdsjaar
C Changed plt to library rather than single source file.
C 
C=======================================================================
      SUBROUTINE LDROTA(AXIS,COSANG,SINANG,MAT)
      REAL MAT(4,4)
      CHARACTER*(*) AXIS
      CHARACTER*1 TAXIS

      CALL MXIDEN(4,MAT)
      TAXIS = AXIS
      CALL CHRUP(TAXIS,TAXIS)
      IF (TAXIS.EQ.'X') THEN
         MAT(2,2) = COSANG
         MAT(2,3) = SINANG
         MAT(3,2) = -SINANG
         MAT(3,3) = COSANG

      ELSE IF (TAXIS.EQ.'Y') THEN
         MAT(1,1) = COSANG
         MAT(1,3) = -SINANG
         MAT(3,1) = SINANG
         MAT(3,3) = COSANG

      ELSE IF (TAXIS.EQ.'Z') THEN
         MAT(1,1) = COSANG
         MAT(1,2) = SINANG
         MAT(2,1) = -SINANG
         MAT(2,2) = COSANG
      END IF

      RETURN

      END
