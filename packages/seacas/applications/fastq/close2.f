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

C $Id: close2.f,v 1.2 1998/07/14 18:18:29 gdsjaar Exp $
C $Log: close2.f,v $
C Revision 1.2  1998/07/14 18:18:29  gdsjaar
C Removed unused variables, cleaned up a little.
C
C Changed BLUE labels to GREEN to help visibility on black background
C (indirectly requested by a couple users)
C
C Revision 1.1.1.1  1990/11/30 11:04:47  gdsjaar
C FASTQ Version 2.0X
C
c Revision 1.1  90/11/30  11:04:46  gdsjaar
c Initial revision
c 
C
CC* FILE: [.PAVING]CLOSE2.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE CLOSE2 (MXND, MLN, NUID, XN, YN, ZN, LXK, KXL, NXL,
     &   LXN, LNODES, IAVAIL, NAVAIL, NNN, LLL, N1, XMIN, XMAX, YMIN,
     &   YMAX, ZMIN, ZMAX, PGRAPH, VIDEO, DEV1, KREG, NOROOM, ERR)
C***********************************************************************
C
C  SUBROUTINE CLOSE2 = SEALS OFF THE LAST 2 OPEN LINES WHILE CHECKING
C                      FOR FORMING A 2-LINE NODE ON THE INTERIOR
C                      (A 2-LINE NODE GENERATES 2 DEGENERATE QUADS)
C
C***********************************************************************
C
      DIMENSION NUID (MXND), XN (MXND), YN (MXND), ZN (MXND)
      DIMENSION LXK (4, MXND), KXL (2, 3*MXND)
      DIMENSION NXL (2, 3*MXND), LXN (4, MXND)
      DIMENSION LNODES (MLN, MXND)
C
      LOGICAL ERR, NOROOM, FOUND, PGRAPH, DONE, CHECK, VIDEO
C
      CHARACTER*3 DEV1
C
      ERR = .FALSE.
      CHECK = .FALSE.
C
      N0 = LNODES (2, N1)
      LINE1 = LNODES (5, N0)
      LINE2 = LNODES (5, N1)
C
C  CHECK TO MAKE SURE THAT AT LEAST ONE OF THE LINES
C  IS NOT A BOUNDARY LINE AND GET THE NODE TO BE DELETED
C
  100 CONTINUE
      IF ((KXL (1, LINE1) .GT. 0) .OR.
     &   (KXL (1, LINE2) .GT. 0)) THEN
C
         FOUND = .TRUE.
C
         IF (KXL (1, LINE1) .GT. 0) THEN
            LNEW = LINE2
            LOLD = LINE1
         ELSE
            LNEW = LINE1
            LOLD = LINE2
         ENDIF
         KOLD = KXL (1, LOLD)
         KNEW = KXL (1, LNEW)
C
C  CHECK FOR ONE OF THE NODES BEING A TWO LINE NODE
C
         IF (KOLD. EQ. KNEW) THEN
            IF (LXN (3, N0) .EQ. 0) THEN
               NGONE = N0
               NTHERE = N1
            ELSEIF (LXN (3, N1) .EQ. 0) THEN
               NGONE = N1
               NTHERE = N0
            ELSE
               CALL MESAGE ('** PROBLEMS WITH NO TWO LINE NODE'//
     &            ' ATTACHED IN CLOSE2 **')
               ERR = .TRUE.
               GOTO 150
            ENDIF
C
C  DELETE THE TWO-LINE NODE, THE TWO LINES, AND THE ELEMENT
C
            KXL (1, LOLD) = 0
            KXL (2, LOLD) = 0
            NXL (1, LOLD) = 0
            NXL (2, LOLD) = 0
            KXL (1, LNEW) = 0
            KXL (2, LNEW) = 0
            NXL (1, LNEW) = 0
            NXL (2, LNEW) = 0
C
C  UNHOOK BOTH LINES FROM NTHERE
C
            CALL DELLXN (MXND, LXN, NUID, NAVAIL, IAVAIL, NTHERE,
     &         LOLD, NNN, ERR, NOROOM)
            IF (ERR) THEN
               CALL MESAGE ('** PROBLEMS IN CLOSE2 DELETING LOLD'//
     &            ' FROM NTHERE **')
               GOTO 150
            ENDIF
            CALL DELLXN (MXND, LXN, NUID, NAVAIL, IAVAIL, NTHERE,
     &         LNEW, NNN, ERR, NOROOM)
            IF (ERR) THEN
               CALL MESAGE ('** PROBLEMS IN CLOSE2 DELETING LNEW'//
     &            ' FROM NTHERE **')
               GOTO 150
            ENDIF
C
C  NOW DELETE NGONE AND THE ELEMENT
C
            DO 110 I = 1, 4
               LXN (I, NGONE) = 0
               IF ( (LXK (I, KOLD) .EQ. LNEW) .OR.
     &            (LXK (I, KOLD) .EQ. LOLD) ) LXK (I, KOLD) = 0
  110       CONTINUE
            LOLD = 0
            LNEW = 0
            DO 120 I = 1, 4
               IF (LXK (I, KOLD) .NE. 0) THEN
                  IF (LOLD .EQ. 0) THEN
                     LOLD = LXK (I, KOLD)
                  ELSE
                     LNEW = LXK (I, KOLD)
                  ENDIF
                  LXK (I, KOLD) = 0
               ENDIF
  120       CONTINUE
            KXL (1, LNEW) = KXL (1, LNEW) + KXL (2, LNEW) - KOLD
            KXL (2, LNEW) = 0
            KXL (1, LOLD) = KXL (1, LOLD) + KXL (2, LOLD) - KOLD
            KXL (2, LOLD) = 0
C
C  NOW RESET THE NECESSARY VARIABLES
C
            N1 = NXL (1, LNEW)
            N0 = NXL (2, LNEW)
            LINE1 = LOLD
            LINE2 = LNEW
            GOTO 100
         ENDIF
C
C  DELETE THE OLD LINE AND REDO LINK ARRAYS
C
         IF (KNEW .EQ. 0) THEN
            KXL (1, LNEW) = KOLD
            KXL (2, LNEW) = 0
         ELSE
            KXL (1, LNEW) = KNEW
            KXL (2, LNEW) = KOLD
         ENDIF
         KXL (1, LOLD) = 0
         KXL (2, LOLD) = 0
         NXL (1, LOLD) = 0
         NXL (2, LOLD) = 0
C
C  FIX THE LINES PER ELEMENT ARRAY FOR THE ONE ELEMENT CHANGING
C
         DO 130 II = 1, 4
            IF (LXK (II, KOLD) .EQ. LOLD) THEN
               LXK (II, KOLD) = LNEW
               GOTO 140
            ENDIF
  130    CONTINUE
         CALL MESAGE ('** PROBLEMS IN CLOSE2 WITH CHANGING ELEMENT **')
         ERR = .TRUE.
         GOTO 150
  140    CONTINUE
C
C  FIX LXN ARRAY
C  UNHOOK LOLD FROM N0 AND FROM N1
C
         CALL DELLXN (MXND, LXN, NUID, NAVAIL, IAVAIL, N0,
     &      LOLD, NNN, ERR, NOROOM)
         IF (ERR) THEN
            CALL MESAGE ('** PROBLEMS IN CLOSE2 DELETING NNN LINES **')
            GOTO 150
         ENDIF
         CALL DELLXN (MXND, LXN, NUID, NAVAIL, IAVAIL, N1,
     &      LOLD, NNN, ERR, NOROOM)
         IF (ERR) THEN
            CALL MESAGE ('** PROBLEMS IN CLOSE2 DELETING N1 LINES **')
            GOTO 150
         ENDIF
C
C NOW FIX THE LNODES ARRAY
C
         LNODES (4, N1) = - 2
         LNODES (4, N0) = - 2
C
      ELSE
         CALL MESAGE ('** PINCHED TOO FAR IN CLOSE2 **')
         GOTO 150
      ENDIF
C
C  NOW SEE IF THE CLOSURE HAS PRODUCED A 2-LINE NODE AND
C  THUS REQUIRES THAT ONE OF THE ELEMENTS MUST BE SQUASHED
C
      IF ((LXN (3, N0) .EQ. 0) .AND. (LXN (2, N0) .GT. 0)) THEN
         CALL DELEM (MXND, XN, YN, NUID, LXK, KXL, NXL, LXN,
     &      NNN, NAVAIL, IAVAIL, N0, KXL (1, LNEW), IDUM1, IDUM2,
     &      DONE, CHECK, NOROOM, ERR)
         IF ((NOROOM) .OR. (ERR)) GOTO 150
      ELSEIF ((LXN (3, N1) .EQ. 0) .AND. (LXN (2, N1) .GT. 0)) THEN
         CALL DELEM (MXND, XN, YN, NUID, LXK, KXL, NXL, LXN,
     &      NNN, NAVAIL, IAVAIL, N1, KXL (1, LNEW), IDUM1, IDUM2,
     &      DONE, CHECK, NOROOM, ERR)
         IF ((NOROOM) .OR. (ERR)) GOTO 150
      ENDIF
C
      IF ( (FOUND) .AND. ((PGRAPH) .OR. (VIDEO)) ) THEN
         CALL RPLOTL (MXND, XN, YN, ZN, NXL, XMIN, XMAX, YMIN, YMAX,
     &      ZMIN, ZMAX, LLL, DEV1, KREG)
         IF (VIDEO) CALL SNAPIT (1)
      ENDIF
C
  150 CONTINUE
C
      RETURN
C
      END
