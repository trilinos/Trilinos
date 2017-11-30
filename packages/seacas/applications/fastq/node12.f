C    Copyright(C) 2014-2017 National Technology & Engineering Solutions of
C    Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
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
C    * Neither the name of NTESS nor the names of its
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

C $Id: node12.f,v 1.1 1990/11/30 11:12:44 gdsjaar Exp $
C $Log: node12.f,v $
C Revision 1.1  1990/11/30 11:12:44  gdsjaar
C Initial revision
C
C
CC* FILE: [.PAVING]NODE12.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE NODE12 (MXND, MLN, LNODES, I1, I2, NLOOP1, NLOOP2,
     &   NODE1, NODE2, NODE, ERR)
C***********************************************************************
C
C  SUBROUTINE NODE12 = FINDS THE CURRENT NODE IN BOTH NEW LOOPS, AND
C                      KEEPS IT A CONSTANT
C
C***********************************************************************
C
      DIMENSION LNODES (MLN, MXND)
C
      LOGICAL ERR
C
      ERR = .FALSE.
C
      KOUNT = 0
      NTEST = I1
  100 CONTINUE
      KOUNT = KOUNT + 1
      IF (KOUNT .GT. NLOOP2) GOTO 110
      IF (NTEST .EQ. NODE) THEN
         NODE1 = NODE
         NODE2 = I2
         GOTO 130
      ENDIF
      NTEST = LNODES (3, NTEST)
      GOTO 100
C
  110 CONTINUE
C
      KOUNT = 0
      NTEST = I2
  120 CONTINUE
      KOUNT = KOUNT + 1
      IF (KOUNT .GT. NLOOP1) THEN
         CALL MESAGE ('** PROBLEMS IN NODE12 FINDING NODE **')
         ERR = .TRUE.
         GOTO 130
      ENDIF
      IF (NTEST .EQ. NODE) THEN
         NODE1 = I1
         NODE2 = NODE
         GOTO 130
      ENDIF
      NTEST = LNODES (3, NTEST)
      GOTO 120
C
  130 CONTINUE
C
      RETURN
C
      END
