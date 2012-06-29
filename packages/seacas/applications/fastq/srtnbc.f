C $Id: srtnbc.f,v 1.1 1990/11/30 11:16:33 gdsjaar Exp $
C $Log: srtnbc.f,v $
C Revision 1.1  1990/11/30 11:16:33  gdsjaar
C Initial revision
C
C
CC* FILE: [.RENUM]SRTNBC.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE SRTNBC (MXNFLG, NPNBC, NNN, NNFLG, NNLEN, NNPTR,
     &   NODES, LSTNBC, IHERE, NNNBC, NBCNOD, NNLIST)
C***********************************************************************
C
C  SUBROUTINE SRTNBC = SORTS THE LIST OF NODAL BOUNDARY FLAGS
C
C***********************************************************************
C
C  VARIABLES USED:
C     IHERE  = AN ATTENDANCE ARRAY TO SEE IF A NODE HAS BEEN FLAGGED
C     NNFLG  = THE ARRAY OF FLAG VALUES
C     NNLEN  = NUMBER OF NODES IN THE NODE LIST ASSOCIATED WITH EACH
C              FLAG
C     NNPTR  = POINTER TO THE FIRST NODE IN LIST FOR EACH FLAG
C     NODES  = THE NODE LIST
C     NNN    = THE NUMBER OF NODES IN THE MESH
C     MXNFLG = THE NUMBER OF ENTRIES IN THE BOUNDARY LIST
C     ENTER  = .TRUE. IF THE FOLLOWING NODES ARE TO BE CHECKED "HERE"
C     FOUND  = .TRUE. IF A NEW UNIQUE FLAG HAS BEEN FOUND
C
C***********************************************************************
C
      DIMENSION NNFLG (MXNFLG), NNLEN (MXNFLG), NNPTR (MXNFLG)
      DIMENSION NODES (NPNBC), LSTNBC (NPNBC), IHERE (NNN)
C
      LOGICAL ENTER, FOUND
C
      NNLIST = 0
      IHOLD = 1
      NBCNOD = 0
C
  100 CONTINUE
      ISTART = IHOLD
      IHOLD = NNNBC
      ENTER = .FALSE.
      FOUND = .FALSE.
C
      DO 110 I = 1, NNN
         IHERE (I) = 0
  110 CONTINUE
C
      DO 120 I = ISTART, NNNBC
         IF (LSTNBC (I) .LT. 0) THEN
            IF (FOUND) THEN
               IF (ENTER)IHOLD = MIN0 (IHOLD, I - 1)
               ITEST = ABS (LSTNBC (I))
               IF (ITEST .EQ. NNFLG (NBCNOD)) THEN
                  ENTER = .TRUE.
                  LSTNBC (I) = 0
               ELSE
                  ENTER = .FALSE.
               ENDIF
            ELSE
               FOUND = .TRUE.
               ENTER = .TRUE.
               NBCNOD = NBCNOD + 1
               NNFLG (NBCNOD) = ABS (LSTNBC (I))
               NNLEN (NBCNOD) = 0
               NNPTR (NBCNOD) = NNLIST + 1
               LSTNBC (I) = 0
            ENDIF
         ELSEIF (LSTNBC (I) .GT. 0) THEN
            IF (ENTER) THEN
               IHERE (LSTNBC (I)) = 1
               LSTNBC (I) = 0
            ENDIF
         ENDIF
  120 CONTINUE
C
      IF (FOUND) THEN
         DO 130 I = 1, NNN
            IF (IHERE (I) .EQ. 1) THEN
               NNLIST = NNLIST + 1
               NNLEN (NBCNOD) = NNLEN (NBCNOD) + 1
               NODES (NNLIST) = I
            ENDIF
  130    CONTINUE
         GOTO 100
      ELSE
         RETURN
      ENDIF
C
      END
