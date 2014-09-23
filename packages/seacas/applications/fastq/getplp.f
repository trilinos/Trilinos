C    Copyright (c) 2014, Sandia Corporation.
C    Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
C    the U.S. Governement retains certain rights in this software.
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

C $Id: getplp.f,v 1.1 1990/11/30 11:08:34 gdsjaar Exp $
C $Log: getplp.f,v $
C Revision 1.1  1990/11/30 11:08:34  gdsjaar
C Initial revision
C
C
CC* FILE: [.RENUM]GETPLP.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE GETPLP (NPNODE, NPELEM, MAXKXN, NNXK, MXLIST, KXN,
     &   NXK, NUID, IP1, LINE, IP2, LIST, NLIST, NNN, LASTN, NOROOM,
     &   ERR)
C***********************************************************************
C
C  SUBROUTINE GETPLP = PRODUCES THE LIST OF NODES FROM POINT IP1
C                      THRU LINE TO POINT IP2
C
C***********************************************************************
C
C  NOTE:
C     THIS LIST WILL BE (LIST (I), I=1,NLIST) AND THESE WILL BE INDICES
C     INTO THE NODE TABLE
C
C***********************************************************************
C
      DIMENSION NXNLST (20)
      DIMENSION KXN (NNXK, MAXKXN), NXK (NNXK, NPELEM), NUID (NPNODE)
      DIMENSION LIST (MXLIST)
C
      LOGICAL ERR, ALL, NOROOM
C
      ERR = .FALSE.
      NOROOM = .FALSE.
C
C  FIND FIRST POINT
C
      IF (NLIST .EQ. 0) THEN
         N = INDX (NNN, NUID, IP1)
         IF (N .LE. 0) THEN
            WRITE ( * , 10000)IP1
            ERR = .TRUE.
            RETURN
         ENDIF
         NLIST = NLIST + 1
         IF (NLIST .GT. MXLIST) THEN
            NOROOM = .TRUE.
            RETURN
         ENDIF
         LIST (NLIST) = N
      ELSE
         NLIST = 0
         N = LASTN
      ENDIF
C
C  FOLLOW THE LINE
C
      IF (LINE .LE. 0) RETURN
      NPREV = 0
  100 CONTINUE
      ALL = .FALSE.
      CALL GETNXN (NPNODE, NPELEM, MAXKXN, NNXK, KXN, NXK, NUID, N,
     &   NXNLST, NUMN, ALL, ERR)
      IF (ERR) RETURN
      DO 110 I = 1, NUMN
         NEW = NXNLST (I)
         IF ( (NEW .NE. NPREV) .AND. (NUID (NEW) .GE. 1000000000)) THEN
            L =  (NUID (NEW) - 1000000000) / 100000
            IF (L .EQ. LINE) THEN
               NLIST = NLIST + 1
               IF (NLIST .GT. MXLIST) THEN
                  NOROOM = .TRUE.
                  RETURN
               ENDIF
               LIST (NLIST) = NEW
               NPREV = N
               N = NEW
               GOTO 100
            ENDIF
         ENDIF
  110 CONTINUE
C
C  LINE FINISHED  -  FIND IP2
C
      IF (IP2 .LE. 0) RETURN
      DO 120 I = 1, NUMN
         NEW = NXNLST (I)
         IF (NUID (NEW) .EQ. IP2) THEN
            NLIST = NLIST + 1
            IF (NLIST .GT. MXLIST) THEN
               NOROOM = .TRUE.
               RETURN
            ENDIF
            LIST (NLIST) = NEW
            LASTN = NEW
            RETURN
         ENDIF
  120 CONTINUE
C
C  LINE DID NOT MATCH UP RIGHT
C
      WRITE ( * , 10010)IP1, LINE, IP2
      ERR = .TRUE.
      RETURN
C
10000 FORMAT (' POINT', I5, ' IS NOT IN THE MESH')
10010 FORMAT (' P-L-P SEQUENCE OF', I5, ' -', I5, ' -', I5,
     &   'IS AN ILLEGAL SEQUENCE')
      END
