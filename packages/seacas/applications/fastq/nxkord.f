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

C $Id: nxkord.f,v 1.1 1990/11/30 11:13:00 gdsjaar Exp $
C $Log: nxkord.f,v $
C Revision 1.1  1990/11/30 11:13:00  gdsjaar
C Initial revision
C
C
CC* FILE: [.QMESH]NXKORD.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE NXKORD (NODES, N1)
C***********************************************************************
C
C  SUBROUTINE NXKORD = ROTATES THE LIST OF FOUR NODES SO N1 APPEARS
C                      FIRST IF IT IS IN THE LIST
C
C***********************************************************************
C
      DIMENSION NODES (4)
C
      DO 100 I = 1, 4
         IF (NODES (I) .EQ. N1) THEN
            IF (I .EQ. 1) THEN
               RETURN
            ELSEIF (I .EQ. 2) THEN
               NSAVE = NODES (1)
               NODES (1) = NODES (2)
               NODES (2) = NODES (3)
               NODES (3) = NODES (4)
               NODES (4) = NSAVE
            ELSEIF (I .EQ. 3) THEN
               NSAVE = NODES (1)
               NSAVE2 = NODES (2)
               NODES (1) = NODES (3)
               NODES (2) = NODES (4)
               NODES (3) = NSAVE
               NODES (4) = NSAVE2
            ELSE
               NSAVE = NODES (4)
               NODES (4) = NODES (3)
               NODES (3) = NODES (2)
               NODES (2) = NODES (1)
               NODES (1) = NSAVE
            ENDIF
            RETURN
         ENDIF
  100 CONTINUE
C
      RETURN
C
      END
