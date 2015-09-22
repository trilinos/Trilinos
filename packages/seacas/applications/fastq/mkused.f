C    Copyright (c) 2014, Sandia Corporation.
C    Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
C    the U.S. Government retains certain rights in this software.
C    
C    Redistribution and use in source and binary forms, with or without
C    modification, are permitted provided that the following conditions are
C    met:
C    
C        * Redistributions of source code must retain the above copyright
C          notice, this list of conditions and the following disclaimer.
C    
C        * Redistributions in binary form must reproduce the above
C          copyright notice, this list of conditions and the following
C          disclaimer in the documentation and/or other materials provided
C          with the distribution.
C    
C        * Neither the name of Sandia Corporation nor the names of its
C          contributors may be used to endorse or promote products derived
C          from this software without specific prior written permission.
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

C $Id: mkused.f,v 1.1 1990/11/30 11:12:17 gdsjaar Exp $
C $Log: mkused.f,v $
C Revision 1.1  1990/11/30 11:12:17  gdsjaar
C Initial revision
C
C
CC* FILE: [.QMESH]MKUSED.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE MKUSED (MXNL, MP, ML, LISTL, IPOINT, NINT, LINKP,
     &   LINKL, LCON, NL)
C***********************************************************************
C
C  SUBROUTINE MKUSED = MARKS ALL LINES AND POINTS USED IN THE PERIMETER
C
C***********************************************************************
C
C  NOTE:
C     THE MESH TABLES ARE EFFECTIVELY DESTROYED BY THIS ROUTINE
C
C***********************************************************************
C
      DIMENSION LISTL (MXNL), IPOINT (MP), NINT (ML)
      DIMENSION LINKP (2, MP), LINKL (2, ML)
      DIMENSION LCON (3, ML)
C
      LOGICAL ADDLNK
C
      ADDLNK = .FALSE.
C
C  FLAG ALL LINES AND POINTS AS BEING USED
C
      DO 100 I = 1, NL
         CALL LTSORT (ML, LINKL, LISTL (I), IL, ADDLNK)
         NINT (IL) =  - IABS (NINT (IL))
         CALL LTSORT (MP, LINKP, LCON (1, IL), IP1, ADDLNK)
         CALL LTSORT (MP, LINKP, LCON (2, IL), IP2, ADDLNK)
         IPOINT (IP1) =  - IABS (IPOINT (IP1))
         IPOINT (IP2) =  - IABS (IPOINT (IP2))
  100 CONTINUE
C
      RETURN
C
      END
