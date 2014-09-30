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

C $Id: savreg.f,v 1.1 1990/11/30 11:15:18 gdsjaar Exp $
C $Log: savreg.f,v $
C Revision 1.1  1990/11/30 11:15:18  gdsjaar
C Initial revision
C
C
CC* FILE: [.QMESH]SAVREG.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE SAVREG (MXND, MAXNBC, MAXSBC, XN, YN, NUID, LXK, NXL,
     &   LXN, LSTNBC, LSTSBC, KNBC, KSBC, NNN, KKK, NUMREG, IUNIT, BAR,
     &   M1)
C***********************************************************************
C
C  SUBROUTINE SAVREG = SAVES THE NODE AND ELEMENT DESCRIPTIONS AS WELL
C                      AS THE BOUNDARY CONDITIONS
C
C***********************************************************************
C
C  NOTE:
C     THE MESH TABLES ARE EFFECTIVELY DESTROYED BY THIS ROUTINE
C
C***********************************************************************
C
      DIMENSION XN (MXND), YN (MXND), NUID (MXND)
      DIMENSION LXK (4, MXND), NXL (2, MXND*3), LXN (4, MXND)
      DIMENSION LSTNBC (MAXNBC), LSTSBC (MAXSBC), NODES (4)
C
      LOGICAL CCW, BAR
C
      CCW = .TRUE.
      IF (.NOT.BAR) THEN
C
C  DEFINE NUID-S FOR INTERIOR NODES.
C  SKIP DELETED NODES AND CONTINUATIONS.
C
         K = 0
         DO 100 I = 1, NNN
            IF ((NUID (I) .EQ. 0) .AND. (LXN (1, I) .GT. 0)) THEN
               K = K+1
               NUID (I) = NUMREG*100000+K
            ENDIF
  100    CONTINUE
C
C  GET COUNTER-CLOCKWISE NODE LISTS.
C   (THESE LISTS WILL OVERWRITE THE LXK ARRAY.)
C  DELETED ELEMENTS WILL BE SKIPPED.
C
         J = 0
         IDEL = 0
         DO 130 K = 1, KKK
            IF (LXK (1, K) .LE. 0) THEN
               DO 110 JJ = 2, KSBC, 3
                  IF (LSTSBC (JJ) .GE. (K - IDEL)) THEN
                     LSTSBC (JJ) = LSTSBC (JJ) - 1
                  ENDIF
  110          CONTINUE
               IDEL = IDEL + 1
            ELSE
               CALL GNXKA (MXND, XN, YN, K, NODES, AREA, LXK, NXL, CCW)
               J = J+1
               DO 120 I = 1, 4
                  N = NODES (I)
                  LXK (I, J) = IABS (NUID (N))
  120          CONTINUE
            ENDIF
  130    CONTINUE
         KKK = J
      ELSE
         DO 140 I = 1, KKK
            LXK (1, I) = IABS (NUID (LXK (1, I)))
            LXK (2, I) = IABS (NUID (LXK (2, I)))
  140    CONTINUE
      ENDIF
C
C  COLLAPSE THE NODE ARRAYS TO ELIMINATE DELETED NODES,
C  CONTINUATIONS,  AND NODES ALREADY WRITTEN OUT.
C
      K = 0
      DO 150 I = 1, NNN
         IF ( ( (LXN (1, I) .GT. 0) .OR. (BAR))
     &      .AND. (NUID (I) .GT. 0) ) THEN
            K = K+1
            XN (K) = XN (I)
            YN (K) = YN (I)
            NUID (K) = NUID (I)
         ENDIF
  150 CONTINUE
      NNN = K
C
C  WRITE HEADER,  NODE LIST,  ELEMENT LIST,  AND BOUNDARY LISTS
C
      WRITE (IUNIT)KKK, NNN, KNBC, KSBC, NUMREG, BAR, M1
      IF (NNN .GE. 1) WRITE (IUNIT) (NUID (I), XN (I), YN (I),
     &   I = 1, NNN)
      WRITE (IUNIT) ((LXK (I, J), I = 1, 4), J = 1, KKK)
      IF (KNBC .GT. 0)WRITE (IUNIT) (LSTNBC (I), I = 1, KNBC)
      IF (KSBC .GT. 0)WRITE (IUNIT) (LSTSBC (I), I = 1, KSBC)
C
      RETURN
C
      END
