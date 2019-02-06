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

C $Id: amxbm.f,v 1.2 1998/07/14 18:18:20 gdsjaar Exp $
C $Log: amxbm.f,v $
C Revision 1.2  1998/07/14 18:18:20  gdsjaar
C Removed unused variables, cleaned up a little.
C
C Changed BLUE labels to GREEN to help visibility on black background
C (indirectly requested by a couple users)
C
C Revision 1.1.1.1  1990/11/30 11:03:33  gdsjaar
C FASTQ Version 2.0X
C
c Revision 1.1  90/11/30  11:03:31  gdsjaar
c Initial revision
c
C
CC* FILE: [.MAIN]AMXBM.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE AMXBM (NPNODE, NPELEM, NXK, AMESUR, BMESUR, KNODE)
C***********************************************************************
C
C  SUBROUTINE AMXBM = ROUTINE TO TRANSFER ELEMENT VARIABLES TO NODES
C
C***********************************************************************
C
      DIMENSION NXK (9, NPELEM)
      DIMENSION AMESUR (NPELEM), BMESUR (NPNODE)
      DIMENSION KNODE (NPNODE)
C
      DO 100 I = 1, NPNODE
         BMESUR(I) = 0.
         KNODE (I) = 0
  100 CONTINUE
C
C  GATHER ALL THE VARIABLES TO THE NODES AND COUNT HOW MANY AT EACH NODE
C
      DO 120 I = 1, NPELEM
         DO 110 J = 1, 4
            NODE = NXK (J, I)
            BMESUR (NODE) = BMESUR (NODE) + AMESUR (I)
            KNODE (NODE) = KNODE (NODE) + 1
  110    CONTINUE
  120 CONTINUE
C
C  GET THE AVERAGE VALUE AT EACH NODE
C
      DO 130 NODE = 1, NPNODE
         BMESUR (NODE) = BMESUR(NODE) / FLOAT (KNODE (NODE))
  130 CONTINUE
C
      RETURN
C
      END
