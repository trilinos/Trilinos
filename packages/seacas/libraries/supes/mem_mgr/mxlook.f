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
      SUBROUTINE MXLOOK (MNGET, VOID, LVOID, NVOIDS, VROW, LASTER)
C
      IMPLICIT INTEGER (A-Z)
      INCLUDE 'params.inc'
C
C     This routine looks for space in the void table.
C
C***********************************************************************
C
C     MNGET    Amount of space requested.
C     VOID     Void table.
C     LVOID    Dimension of void table.
C     NVOIDS   Number of voids.
               DIMENSION VOID(LVOID,2)
C     VROW     Row number that contains void to satisfy space request.
C     LASTER   Error return.
C
C***********************************************************************
C
C     CHECK TO SEE IF A VOID WILL CONTAIN THE MEMORY REQUEST.
C
      VROW = 0
      VLEN = 0
      DO 100 I = 1, NVOIDS
         IF (VOID(I,2) .GE. MNGET) THEN
C
C           THIS VOID HAS ENOUGH ROOM - FIND THE SMALLEST VOID THAT
C           IS LARGE ENOUGH.
C
            IF (VLEN .EQ. 0 .OR. VOID(I,2) .LT. VLEN) THEN
               VROW = I
               VLEN = VOID(I,2)
            END IF
         END IF
  100 CONTINUE
      IF (VROW .NE. 0) THEN
         LASTER = SUCESS
      ELSE
         LASTER = NOGET
      END IF
      RETURN
      END
