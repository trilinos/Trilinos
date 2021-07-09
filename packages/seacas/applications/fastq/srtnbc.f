C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE SRTNBC (MXNFLG, NPNBC, NNN, NNFLG, NNLEN, NNPTR,
     &   NODES, LSTNBC, IHERE, NNNBC, NBCNOD, NNLIST)
C***********************************************************************

C  SUBROUTINE SRTNBC = SORTS THE LIST OF NODAL BOUNDARY FLAGS

C***********************************************************************

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

C***********************************************************************

      DIMENSION NNFLG (MXNFLG), NNLEN (MXNFLG), NNPTR (MXNFLG)
      DIMENSION NODES (NPNBC), LSTNBC (NPNBC), IHERE (NNN)

      LOGICAL ENTER, FOUND

      NNLIST = 0
      IHOLD = 1
      NBCNOD = 0

  100 CONTINUE
      ISTART = IHOLD
      IHOLD = NNNBC
      ENTER = .FALSE.
      FOUND = .FALSE.

      DO 110 I = 1, NNN
         IHERE (I) = 0
  110 CONTINUE

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

      END
