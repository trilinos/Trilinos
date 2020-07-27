C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE WRABQS (MS, MR, NPNODE, NPELEM, MXNFLG, MXSFLG, NPREGN,
     &   NPNBC, NPSBC, IUNIT, NNN, KKK, NNXK, NODES, NELEMS, NNFLG,
     &   NNPTR, NNLEN, NSFLG, NSPTR, NSLEN, NVPTR, NVLEN, NSIDEN,
     &   MAPDXG, XN, YN, NXK, MAT, MAPGXD, MATMAP, NBCNOD, NNLIST,
     &   NBCSID, NSLIST, NVLIST, NUMMAT, LINKM, TITLE, ERR, EIGHT, NINE)
C***********************************************************************

C  SUBROUTINE WRABQS = WRITES ABAQUS DATABASE MESH OUTPUT FILE

C***********************************************************************

      DIMENSION XN (NPNODE), YN (NPNODE), NXK (NNXK, NPELEM)
      DIMENSION MAT (NPELEM)
      DIMENSION NODES (NPNBC), NELEMS (NPSBC), NSIDEN (NPSBC)
      DIMENSION NNFLG (MXNFLG), NNLEN (MXNFLG), NNPTR (MXNFLG)
      DIMENSION NSFLG (MXSFLG), NSLEN (MXSFLG), NSPTR (MXSFLG)
      DIMENSION NVLEN (MXSFLG), NVPTR (MXSFLG), LINKM (2,  (MS+MR))
      DIMENSION MAPDXG (NPNODE), MAPGXD (NPNODE), MATMAP (3, NPREGN)
      DIMENSION IHOLD (9)

      CHARACTER*72 TITLE, DUMMY, DUMMY2

      LOGICAL ERR, EIGHT, NINE, DEFTYP, FOUND

      ERR = .TRUE.

C  WRITE OUT HEADER TITLE AND INFORMATION

      WRITE (IUNIT, 10000, ERR = 200)TITLE
      WRITE (IUNIT, 10010, ERR = 200)NNN, KKK, NBCNOD, NBCSID

C  WRITE OUT NODE BLOCK

      WRITE (IUNIT, 10020, ERR = 200)
      Z = 0.
      DO 100 I = 1, NNN
         WRITE (IUNIT, 10030, ERR = 200)I, XN (I), YN (I), Z
  100 CONTINUE

C  QUERY THE USER FOR LOCAL CONTROL OF ELEMENT TYPE

      CALL INQTRU ('USE DEFAULT ELEMENT TYPES FOR ELSETS', DEFTYP)

C  WRITE OUT ELEMENT BLOCKS

      DO 130 I = 1, NUMMAT
         CALL GETDUM (MATMAP (1, I), DUMMY, LEN)
         IF (NXK (3, MATMAP (2, I)) .EQ. 0) THEN
            IF (DEFTYP) THEN
               WRITE (IUNIT, 10060, ERR = 200)DUMMY (1:LEN)
            ELSE
               WRITE (*, 10040)MATMAP (1, I)
               CALL INQSTR ('NEW 2 NODE ELEMENT TYPE:  ', DUMMY2)
               CALL STRLNG (DUMMY2, LEN2)
               WRITE (IUNIT, 10050, ERR = 200)DUMMY2 (1:LEN2),
     &            DUMMY (1:LEN)
            ENDIF
            INODE = 2
         ELSEIF (NXK (4, MATMAP (2, I)) .EQ. 0) THEN
            IF (DEFTYP) THEN
               WRITE (IUNIT, 10070, ERR = 200)DUMMY (1:LEN)
            ELSE
               WRITE (*, 10040)MATMAP (1, I)
               CALL INQSTR ('NEW 3 NODE ELEMENT TYPE:  ', DUMMY2)
               CALL STRLNG (DUMMY2, LEN2)
               WRITE (IUNIT, 10050, ERR = 200)DUMMY2 (1:LEN2),
     &            DUMMY (1:LEN)
            ENDIF
            INODE = 3
         ELSEIF (EIGHT) THEN
            IF (DEFTYP) THEN
               WRITE (IUNIT, 10080, ERR = 200)DUMMY (1:LEN)
            ELSE
               WRITE (*, 10040)MATMAP (1, I)
               CALL INQSTR ('NEW 8 NODE ELEMENT TYPE:  ', DUMMY2)
               CALL STRLNG (DUMMY2, LEN2)
               WRITE (IUNIT, 10050, ERR = 200)DUMMY2 (1:LEN2),
     &            DUMMY (1:LEN)
            ENDIF
            INODE = 8
         ELSEIF (NINE) THEN
            IF (DEFTYP) THEN
               WRITE (IUNIT, 10090, ERR = 200)DUMMY (1:LEN)
            ELSE
               WRITE (*, 10040)MATMAP (1, I)
               CALL INQSTR ('NEW 9 NODE ELEMENT TYPE:  ', DUMMY2)
               CALL STRLNG (DUMMY2, LEN2)
               WRITE (IUNIT, 10050, ERR = 200)DUMMY2 (1:LEN2),
     &            DUMMY (1:LEN)
            ENDIF
            INODE = 9
         ELSE
            IF (DEFTYP) THEN
               WRITE (IUNIT, 10100, ERR = 200)DUMMY (1:LEN)
            ELSE
               WRITE (*, 10040)MATMAP (1, I)
               CALL INQSTR ('NEW 4 NODE ELEMENT TYPE:  ', DUMMY2)
               CALL STRLNG (DUMMY2, LEN2)
               WRITE (IUNIT, 10050, ERR = 200)DUMMY2 (1:LEN2),
     &            DUMMY (1:LEN)
            ENDIF
            INODE = 4
         ENDIF
         IF (INODE .EQ. 8) THEN
            DO 110 KELEM = MATMAP (2, I), MATMAP (3, I)
               K = MAPGXD(KELEM)
               WRITE (IUNIT, 10110, ERR = 200) K, NXK (1, KELEM),
     &            NXK (3, KELEM), NXK (5, KELEM), NXK (7, KELEM),
     &            NXK (2, KELEM), NXK (4, KELEM), NXK (6, KELEM),
     &            NXK (8, KELEM)
  110       CONTINUE
         ELSE
            DO 120 KELEM = MATMAP (2, I), MATMAP (3, I)
               K = MAPGXD(KELEM)
               WRITE (IUNIT, 10120, ERR = 200) K, (NXK (J, KELEM),
     &            J = 1, INODE)
  120       CONTINUE
         ENDIF
  130 CONTINUE

C  WRITE OUT THE NODAL BOUNDARY CONDITIONS

      IF (NBCNOD.GT.0) THEN
         DO 140 I = 1, NBCNOD
            J1 = NNPTR (I)
            J2 = NNPTR (I)+NNLEN (I)-1
            CALL GETDUM (NNFLG (I), DUMMY, LEN)
            WRITE (IUNIT, 10130, ERR = 200)DUMMY (1:LEN)
            WRITE (IUNIT, 10110, ERR = 200) (NODES (J), J = J1, J2)
  140    CONTINUE
      ENDIF

C  WRITE OUT THE SIDE BOUNDARY FLAGS

      IF (NBCSID.GT.0) THEN
C         CALL MESAGE ('ELEMENT NUMBERING IS WRITTEN WITH ELEMENT' //
C     &      BOUNDARY FLAGS')
C         CALL INQTRU ('WOULD YOU LIKE TO CHANGE THIS TO NODES', IANS)
C         IF (IANS) THEN
C            DO 190 I = 1, NBCSID
C               J1 = NVPTR (I)
C               J2 = NVPTR (I) + NVLEN (I)-1
C               CALL GETDUM (NSFLG (I), DUMMY, LEN)
C               WRITE (IUNIT, 180, ERR = 200) DUMMY (1:LEN)
C               WRITE (IUNIT, 160, ERR = 200) (NSIDEN (J), J = J1, J2)
C 190        CONTINUE
C         ELSE
         DO 190 I = 1, NBCSID
            J1 = NSPTR (I)
            J2 = NSPTR (I)+NSLEN (I)-1
            CALL GETDUM (NSFLG (I), DUMMY, LEN)

C  WRITE OUT THE SIDE 1 ELEMENTS

            FOUND = .FALSE.
            JHOLD = 0
            DO 150 J = J1, J2
               JJ1 = NSIDEN ( (J * 2) - 1)
               JJ2 = NSIDEN (J * 2)
               K = NELEMS (J)
               IF (NXK (3, K) .EQ. 0) THEN
                  INODE = 2
               ELSEIF (NXK (4, K) .EQ. 0) THEN
                  INODE = 3
               ELSEIF (EIGHT .OR. NINE) THEN
                  INODE = 8
               ELSE
                  INODE = 4
               ENDIF
               IF ( ( (INODE .EQ. 4) .AND.
     &            ( ( (JJ1 .EQ. NXK (1, K)) .AND.
     &            (JJ2 .EQ. NXK (2, K)) ) .OR.
     &            ( (JJ2 .EQ. NXK (1, K)) .AND.
     &            (JJ1 .EQ. NXK (2, K)) ) ) ) .OR.
     &            ( (INODE .EQ. 8) .AND.
     &            ( ( (JJ1 .EQ. NXK (1, K)) .AND.
     &            (JJ2 .EQ. NXK (2, K)) ) .OR.
     &            ( (JJ2 .EQ. NXK (1, K)) .AND.
     &            (JJ1 .EQ. NXK (2, K)) ) ) ) ) THEN
C     &            ( (JJ1 .EQ. NXK (2, K)) .AND.
C     &            (JJ2 .EQ. NXK (3, K)) ) .OR.
C     &            ( (JJ2 .EQ. NXK (2, K)) .AND.
C     &            (JJ1 .EQ. NXK (3, K)) ) ) ) ) THEN

                  IF (.NOT. FOUND) THEN
                     WRITE (IUNIT, 10150, ERR = 200)
     &                  DUMMY (1:LEN) // '_1'
                     FOUND = .TRUE.
                  ENDIF

                  JHOLD = JHOLD + 1
                  IHOLD (JHOLD) = MAPGXD (K)
                  IF (JHOLD .EQ. 9) THEN
                     WRITE (IUNIT, 10110, ERR = 200)
     &                  (IHOLD (II), II = 1, 9)
                     JHOLD = 0
                  ENDIF
               ENDIF
  150       CONTINUE
            IF (JHOLD .GT. 0) WRITE (IUNIT, 10110, ERR = 200)
     &         (IHOLD (II), II = 1, JHOLD)

C  WRITE OUT THE SIDE 2 ELEMENTS

            FOUND = .FALSE.
            JHOLD = 0
            DO 160 J = J1, J2
               JJ1 = NSIDEN ( (J * 2) - 1)
               JJ2 = NSIDEN (J * 2)
               K = NELEMS (J)
               IF (NXK (3, K) .EQ. 0) THEN
                  INODE = 2
               ELSEIF (NXK (4, K) .EQ. 0) THEN
                  INODE = 3
               ELSEIF (EIGHT .OR. NINE) THEN
                  INODE = 8
               ELSE
                  INODE = 4
               ENDIF
               IF ( ( (INODE .EQ. 4) .AND.
     &            ( ( (JJ1 .EQ. NXK (2, K)) .AND.
     &            (JJ2 .EQ. NXK (3, K)) ) .OR.
     &            ( (JJ2 .EQ. NXK (2, K)) .AND.
     &            (JJ1 .EQ. NXK (3, K)) ) ) ) .OR.
     &            ( (INODE .EQ. 8) .AND.
     &            ( ( (JJ1 .EQ. NXK (3, K)) .AND.
     &            (JJ2 .EQ. NXK (4, K)) ) .OR.
     &            ( (JJ2 .EQ. NXK (3, K)) .AND.
     &            (JJ1 .EQ. NXK (4, K)) ) ) ) ) THEN
C     &            ( (JJ1 .EQ. NXK (4, K)) .AND.
C     &            (JJ2 .EQ. NXK (5, K)) ) .OR.
C     &            ( (JJ2 .EQ. NXK (4, K)) .AND.
C     &            (JJ1 .EQ. NXK (5, K)) ) ) ) ) THEN

                  IF (.NOT. FOUND) THEN
                     WRITE (IUNIT, 10150, ERR = 200)
     &                  DUMMY (1:LEN) // '_2'
                     FOUND = .TRUE.
                  ENDIF

                  JHOLD = JHOLD + 1
                  IHOLD (JHOLD) = MAPGXD (K)
                  IF (JHOLD .EQ. 9) THEN
                     WRITE (IUNIT, 10110, ERR = 200)
     &                  (IHOLD (II), II = 1, 9)
                     JHOLD = 0
                  ENDIF
               ENDIF
  160       CONTINUE
            IF (JHOLD .GT. 0) WRITE (IUNIT, 10110, ERR = 200)
     &         (IHOLD (II), II = 1, JHOLD)

C  WRITE OUT THE SIDE 3 ELEMENTS

            FOUND = .FALSE.
            JHOLD = 0
            DO 170 J = J1, J2
               JJ1 = NSIDEN ( (J * 2) - 1)
               JJ2 = NSIDEN (J * 2)
               K = NELEMS (J)
               IF (NXK (3, K) .EQ. 0) THEN
                  INODE = 2
               ELSEIF (NXK (4, K) .EQ. 0) THEN
                  INODE = 3
               ELSEIF (EIGHT .OR. NINE) THEN
                  INODE = 8
               ELSE
                  INODE = 4
               ENDIF
               IF ( ( (INODE .EQ. 4) .AND.
     &            ( ( (JJ1 .EQ. NXK (3, K)) .AND.
     &            (JJ2 .EQ. NXK (4, K)) ) .OR.
     &            ( (JJ2 .EQ. NXK (3, K)) .AND.
     &            (JJ1 .EQ. NXK (4, K)) ) ) ) .OR.
     &            ( (INODE .EQ. 8) .AND.
     &            ( ( (JJ1 .EQ. NXK (5, K)) .AND.
     &            (JJ2 .EQ. NXK (6, K)) ) .OR.
     &            ( (JJ2 .EQ. NXK (5, K)) .AND.
     &            (JJ1 .EQ. NXK (6, K)) ) ) ) ) THEN
C     &            ( (JJ1 .EQ. NXK (6, K)) .AND.
C     &            (JJ2 .EQ. NXK (7, K)) ) .OR.
C     &            ( (JJ2 .EQ. NXK (6, K)) .AND.
C     &            (JJ1 .EQ. NXK (7, K)) ) ) ) ) THEN

                  IF (.NOT. FOUND) THEN
                     WRITE (IUNIT, 10150, ERR = 200)
     &                  DUMMY (1:LEN) // '_3'
                     FOUND = .TRUE.
                  ENDIF

                  JHOLD = JHOLD + 1
                  IHOLD (JHOLD) = MAPGXD (K)
                  IF (JHOLD .EQ. 9) THEN
                     WRITE (IUNIT, 10110, ERR = 200)
     &                  (IHOLD (II), II = 1, 9)
                     JHOLD = 0
                  ENDIF
               ENDIF
  170       CONTINUE
            IF (JHOLD .GT. 0) WRITE (IUNIT, 10110, ERR = 200)
     &         (IHOLD (II), II = 1, JHOLD)

C  WRITE OUT THE SIDE 4 ELEMENTS

            FOUND = .FALSE.
            JHOLD = 0
            DO 180 J = J1, J2
               JJ1 = NSIDEN ( (J * 2) - 1)
               JJ2 = NSIDEN (J * 2)
               K = NELEMS (J)
               IF (NXK (3, K) .EQ. 0) THEN
                  INODE = 2
               ELSEIF (NXK (4, K) .EQ. 0) THEN
                  INODE = 3
               ELSEIF (EIGHT .OR. NINE) THEN
                  INODE = 8
               ELSE
                  INODE = 4
               ENDIF
               IF ( ( (INODE .EQ. 4) .AND.
     &            ( ( (JJ1 .EQ. NXK (4, K)) .AND.
     &            (JJ2 .EQ. NXK (1, K)) ) .OR.
     &            ( (JJ2 .EQ. NXK (4, K)) .AND.
     &            (JJ1 .EQ. NXK (1, K)) ) ) ) .OR.
     &            ( (INODE .EQ. 8) .AND.
     &            ( ( (JJ1 .EQ. NXK (7, K)) .AND.
     &            (JJ2 .EQ. NXK (8, K)) ) .OR.
     &            ( (JJ2 .EQ. NXK (7, K)) .AND.
     &            (JJ1 .EQ. NXK (8, K)) ) ) ) ) THEN
C     &            ( (JJ1 .EQ. NXK (8, K)) .AND.
C     &            (JJ2 .EQ. NXK (1, K)) ) .OR.
C     &            ( (JJ2 .EQ. NXK (8, K)) .AND.
C     &            (JJ1 .EQ. NXK (1, K)) ) ) ) ) THEN

                  IF (.NOT. FOUND) THEN
                     WRITE (IUNIT, 10150, ERR = 200)
     &                  DUMMY (1:LEN) // '_4'
                     FOUND = .TRUE.
                  ENDIF

                  JHOLD = JHOLD + 1
                  IHOLD (JHOLD) = MAPGXD (K)
                  IF (JHOLD .EQ. 9) THEN
                     WRITE (IUNIT, 10110, ERR = 200)
     &                  (IHOLD (II), II = 1, 9)
                     JHOLD = 0
                  ENDIF
               ENDIF
  180       CONTINUE
            IF (JHOLD .GT. 0) WRITE (IUNIT, 10110, ERR = 200)
     &         (IHOLD (II), II = 1, JHOLD)

  190    CONTINUE
      ENDIF
      CALL MESAGE ('ABAQUS OUTPUT FILE SUCCESSFULLY WRITTEN')
      ERR = .FALSE.
      RETURN

C  ERR DURING WRITE PROBLEMS

  200 CONTINUE
      CALL MESAGE ('ERR DURING WRITE TO ABAQUS OUTPUT FILE')
      CALL MESAGE ('         - NO FILE SAVED -            ')
      RETURN

10000 FORMAT ('*HEADING', /, A72)
10010 FORMAT ('**', /,
     &   '**     MESH GENERATED USING FASTQ        ', /,
     &   '**  NUMBER OF NODES:                     ', I5, /,
     &   '**  NUMBER OF ELEMENTS:                  ', I5, /,
     &   '**  NUMBER OF NODAL BOUNDARY CONDITIONS: ', I5, /,
     &   '**  NUMBER OF SIDE BOUNDARY CONDITIONS:  ', I5, /,
     &   '**')
10020 FORMAT ('*NODE')
10030 FORMAT (I10, ',', 2 (E14.7, ','), E14.7)
10040 FORMAT (' FOR ELEMENT BLOCK ID:', I7)
10050 FORMAT ('*ELEMENT, TYPE = ', A, ', ELSET = MAT', A)
10060 FORMAT ('*ELEMENT, TYPE = C1D2, ELSET = MAT', A)
10070 FORMAT ('*ELEMENT, TYPE = C1D3, ELSET = MAT', A)
10080 FORMAT ('*ELEMENT, TYPE = CPE8, ELSET = MAT', A)
10090 FORMAT ('*ELEMENT, TYPE = INTER9, ELSET = MAT', A)
10100 FORMAT ('*ELEMENT, TYPE = CPE4, ELSET = MAT', A)
10110 FORMAT (8 (I8, ','), I8)
10120 FORMAT (4 (I10, ','), I10)
10130 FORMAT ('*NSET, NSET = NB', A)
10150 FORMAT ('*ELSET, ELSET = EB', A)

      END
