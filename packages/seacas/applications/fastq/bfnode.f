C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE BFNODE (NLIST, NNXK, NPNODE, NPELEM, MAXKXN, NNUID,
     &   NODE, NEWNOD, LIST, KXN, NXK, NUID, JLOC, LINE1, ERR)
C***********************************************************************

C  SUBROUTINE BFNODE = FINDS ANY NODES IN A PORTION OF THE NODAL FLAG
C                      LIST WHICH IS ATTACHED TO THE GIVEN NODE.

C***********************************************************************

C  SUBROUTINE CALLED BY:
C     ADDWT = ADDS THE WEIGHTING FACTORS TO ANY NODES WITH
C             FLAGS CONTAINING WEIGHTS

C***********************************************************************

      DIMENSION JLIST (20), NUID (NNUID)
      DIMENSION LIST (NLIST), NXK (NNXK, NPELEM), KXN (NNXK, MAXKXN)

      LOGICAL ERR, FOUND, ALL, CLOSED, NOTIFY

      ERR = .TRUE.
      FOUND = .FALSE.
      CLOSED = .FALSE.
      NOTIFY = .TRUE.
      NEWNOD = 0
      JLOC = 0

      ALL = .FALSE.
      CALL GETNXN (NPNODE, NPELEM, MAXKXN, NNXK, KXN, NXK, NUID, NODE,
     &   JLIST, IFOUND, ALL, ERR)
      IF (ERR) THEN
         WRITE ( * , 10000)NODE
         RETURN
      ENDIF

C  SEE IF ANY OF THE FOUND NODES ARE IN THE FLAGGED NODE LIST
C  ONLY 1 SHOULD BE  (REPEATS OF THE SAME NODE ARE OK FOR SIDEBC)

      NEWNOD = 0
      NODOLD = 0
      DO 110 I = 1, IFOUND
         DO 100 J = 1, NLIST
            IF (LIST (J) .EQ. JLIST (I)) THEN
               IF ( (FOUND) .AND. (LIST (J) .NE. NEWNOD)) THEN
                  IF ( (CLOSED) .AND. (NODOLD .NE. LIST (J))) THEN
                     WRITE ( * , 10010)NODE
                     RETURN
                  ELSEIF (.NOT.CLOSED) THEN

C  ASSUME IN A CLOSED LOOP THAT THE FIRST LINE IN THE SIDEBC DEFINITION
C   (LINE1) INDICATES THE APPROPRIATE DIRECTION

                     LINET1 =  (NUID (LIST (J)) - 1000000000) / 100000
                     LINET2 =  (NUID (NEWNOD) - 1000000000) / 100000
                     IF ( (LINET1 .EQ. LINET2) .AND.
     &                  (LINET1 .EQ. LINE1)) THEN
                        CLOSED = .TRUE.
                        IF (NOTIFY) WRITE (*, 10020)NODE, NEWNOD, LINE1
                        CALL MESAGE ('NOTE - NEITHER 2ND NODE IS ON '//
     &                     'FIRST LINE')
                     ELSEIF (LINE1 .EQ. LINET1) THEN
                        NODOLD = NEWNOD
                        NEWNOD = LIST (J)
                        JLOC = J
                        CLOSED = .TRUE.
                        IF (NOTIFY) WRITE (*, 10020)NODE, NEWNOD, LINE1
                        NOTIFY = .FALSE.
                     ELSEIF (LINE1 .EQ. LINET2) THEN
                        NODOLD = LIST (J)
                        CLOSED = .TRUE.
                        IF (NOTIFY) WRITE (*, 10020)NODE, NEWNOD, LINE1
                        NOTIFY = .FALSE.
                     ELSE
                        WRITE ( * , 10030)NODE
                        RETURN
                     ENDIF
                  ENDIF
               ELSE
                  FOUND = .TRUE.
                  NEWNOD = LIST (J)
                  JLOC = J
               ENDIF
            ENDIF
  100    CONTINUE
  110 CONTINUE
      ERR = .FALSE.

      RETURN

10000 FORMAT (' ERROR GETTING NODES ATTACHED TO NODE:', I5)
10010 FORMAT (' ERROR - THREE NODES HAVE BEEN FOUND FOR SEQUENCE FROM',
     &   ' NODE:', I5)
10020 FORMAT (' WARNING - CLOSED LOOP FOUND AT BEGINNING POINT OF',
     &   ' WEIGHTING',  /,
     &   ' BEGINNING POINT CORRESPONDS TO NODE:', I5,  /,
     &   ' NODE:', I5, ' ON LINE:', I5,
     &   ' USED AS SECOND WEIGHTING NODE')
10030 FORMAT (' ERROR - FOR CLOSED LOOP FOUND AT BEGINNING POINT OF',
     &   ' WEIGHTING, ',  /,
     &   'POSSIBLE SECOND NODES DO NOT LIE ON THE FIRST LINE:', I5,  /,
     &   'ATTACHED TO THE SIDEBC - DIRECTION IS THUS UNDETERMINABLE.')

      END
