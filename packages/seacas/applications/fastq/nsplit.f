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

C $Id: nsplit.f,v 1.2 1991/05/10 17:42:09 gdsjaar Exp $
C $Log: nsplit.f,v $
C Revision 1.2  1991/05/10 17:42:09  gdsjaar
C Changed VMS JNINT intrinsic to ANSI NINT
C
c Revision 1.1.1.1  1990/11/30  11:12:56  gdsjaar
c FASTQ Version 2.0X
c
c Revision 1.1  90/11/30  11:12:54  gdsjaar
c Initial revision
c 
C
CC* FILE: [.PAVING]NSPLIT.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE NSPLIT (MXND, MLN, LNODES, ANGLE, NSTART, KANG,
     &   INODE, NNODE, NWANT, MAXSIZ)
C***********************************************************************
C
C  SUBROUTINE NSPLIT = SPLITS UP THE KANG CONSECUTIVE NODES STARTING
C                      AT NSTART INTO NWANT INTERVALS (OR AS CLOSE
C                      AS AN BE EXPECTED).  THE MAXIMUM NWANT SHOULD
C                      BE IS 4.
C
C***********************************************************************
C
      DIMENSION LNODES (MLN, MXND), ANGLE (MXND), INODE (4)
C
      LOGICAL MAXSIZ
C
      NNODE = 0
C
      IF (KANG .LE. NWANT) THEN
         NNOW = NSTART
         DO 100 I = 1, KANG
            INODE (I) = NNOW
            NNOW = LNODES (3, NNOW)
  100    CONTINUE
         NNODE = KANG
C
      ELSEIF (NWANT .EQ. 1) THEN
         NNODE = 1
         IF (KANG .EQ. 2) THEN
            IF (MAXSIZ) THEN
               IF (ANGLE (NSTART) .GT. ANGLE (LNODES (3, NSTART)) ) THEN
                  INODE (1) = NSTART
               ELSE
                  INODE (1) = LNODES (3, NSTART)
               ENDIF
            ELSE
               IF (ANGLE (NSTART) .GT. ANGLE (LNODES (3, NSTART)) ) THEN
                  INODE (1) = LNODES (3, NSTART)
               ELSE
                  INODE (1) = NSTART
               ENDIF
            ENDIF
         ELSE
            INODE (1) = JUMPLP (MXND, MLN, LNODES, NSTART, KANG / 2)
         ENDIF
C
      ELSEIF (NWANT .EQ. 2) THEN
         NNODE = 2
         NJUMP = NINT (FLOAT (KANG + 1) / 4.)
         INODE (1) = JUMPLP (MXND, MLN, LNODES, NSTART,
     &      NJUMP - 1)
         INODE (2) = JUMPLP (MXND, MLN, LNODES, NSTART,
     &      KANG - NJUMP)
C
      ELSEIF (NWANT .EQ. 3) THEN
         NNODE = 3
         NJUMP1 = NINT (FLOAT (KANG + 1) / 6.)
         NJUMP2 = NINT (FLOAT (KANG + 1) / 2.)
         INODE (1) = JUMPLP (MXND, MLN, LNODES, NSTART,
     &      NJUMP1 - 1)
         INODE (2) = JUMPLP (MXND, MLN, LNODES, NSTART,
     &      NJUMP2 - 1)
         INODE (3) = JUMPLP (MXND, MLN, LNODES, NSTART,
     &      KANG - NJUMP1)
C
      ELSEIF (NWANT .EQ. 4) THEN
         NNODE = 4
         XKANG = KANG + 1
         NJUMP1 = NINT (XKANG / 8.) - 1
         NJUMP2 = NINT (XKANG / 2.) - NINT (XKANG / 8.) -1
         NJUMP3 = NINT (XKANG / 2.) + NINT (XKANG / 8.) -1
         INODE (1) = JUMPLP (MXND, MLN, LNODES, NSTART,
     &      NJUMP1)
         INODE (2) = JUMPLP (MXND, MLN, LNODES, NSTART,
     &      NJUMP2)
         INODE (3) = JUMPLP (MXND, MLN, LNODES, NSTART,
     &      NJUMP3)
         INODE (4) = JUMPLP (MXND, MLN, LNODES, NSTART,
     &      KANG - NJUMP1 - 1)
      ENDIF
C
      RETURN
C
      END
