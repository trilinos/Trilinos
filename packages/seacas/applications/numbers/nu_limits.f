C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE LIMITS (XYZMIN, XYZMAX, CRD, IX, MAT, NDIM, NEBLK,
     *   NNODES, EXODUS, TIME, ITMSEL, CORDSP, NUMNP)
      DIMENSION CRD(NUMNP, *), IX(NNODES,*), MAT(6,*),
     *   XYZMIN(NDIM,NEBLK), XYZMAX(NDIM,NEBLK),
     *   TIME(*), CORDSP(NUMNP, *)
      LOGICAL ITMSEL(*), ISABRT

      CHARACTER*16 ENGNOT, ENG1
      DIMENSION OVMIN (3), OVMAX (3)
      LOGICAL EXODUS
      include 'nu_io.blk'
      include 'nu_ptim.blk'

C ... IF NOT EXODUS, THEN CORDSP CONTAINS COORDINATES

      IF (EXODUS) THEN
         CALL GETDSP (CRD, CORDSP, NDIM, NUMNP, TIME, ITMSEL,
     *      'R', ISTAT)
         IF (ISTAT .NE. 0) GO TO 130
      END IF

   10 CONTINUE
      IF (EXODUS) THEN
         CALL GETDSP (CRD, CORDSP, NDIM, NUMNP, TIME, ITMSEL,
     *      'A', ISTAT)
         IF (ISTAT .NE. 0) GO TO 130
      END IF
      IF (ISABRT()) RETURN

      IF (NDIM .EQ. 2) THEN
         DO 40 IBLK = 1, NEBLK
            IF (MAT(5,IBLK) .NE. 1) GOTO 40
            IELBEG = MAT(3,IBLK)
            IELEND = MAT(4,IBLK)

            XYZMIN(1,IBLK) = CORDSP(IX(1,IELBEG),1)
            XYZMIN(2,IBLK) = CORDSP(IX(1,IELBEG),2)
            XYZMAX(1,IBLK) = CORDSP(IX(1,IELBEG),1)
            XYZMAX(2,IBLK) = CORDSP(IX(1,IELBEG),2)

            DO 30 IEL = IELBEG, IELEND
               DO 20 I = 1, NNODES
                  XYZMIN(1,IBLK) = MIN(XYZMIN(1,IBLK),
     *               CORDSP(IX(I,IEL),1))
                  XYZMIN(2,IBLK) = MIN(XYZMIN(2,IBLK),
     *               CORDSP(IX(I,IEL),2))
                  XYZMAX(1,IBLK) = MAX(XYZMAX(1,IBLK),
     *               CORDSP(IX(I,IEL),1))
                  XYZMAX(2,IBLK) = MAX(XYZMAX(2,IBLK),
     *               CORDSP(IX(I,IEL),2))
   20          CONTINUE
   30       CONTINUE
   40    CONTINUE

      ELSE
         DO 70 IBLK = 1, NEBLK
            IF (MAT(5,IBLK) .NE. 1) GOTO 70
            IELBEG = MAT(3,IBLK)
            IELEND = MAT(4,IBLK)

            XYZMIN(1,IBLK) = CORDSP(IX(1,IELBEG),1)
            XYZMIN(2,IBLK) = CORDSP(IX(1,IELBEG),2)
            XYZMIN(3,IBLK) = CORDSP(IX(1,IELBEG),3)
            XYZMAX(1,IBLK) = CORDSP(IX(1,IELBEG),1)
            XYZMAX(2,IBLK) = CORDSP(IX(1,IELBEG),2)
            XYZMAX(3,IBLK) = CORDSP(IX(1,IELBEG),3)

            DO 60 IEL = IELBEG, IELEND
               DO 50 I = 1, NNODES
                  XYZMIN(1,IBLK) = MIN(XYZMIN(1,IBLK),
     *               CORDSP(IX(I,IEL),1))
                  XYZMIN(2,IBLK) = MIN(XYZMIN(2,IBLK),
     *               CORDSP(IX(I,IEL),2))
                  XYZMIN(3,IBLK) = MIN(XYZMIN(3,IBLK),
     *               CORDSP(IX(I,IEL),3))
                  XYZMAX(1,IBLK) = MAX(XYZMAX(1,IBLK),
     *               CORDSP(IX(I,IEL),1))
                  XYZMAX(2,IBLK) = MAX(XYZMAX(2,IBLK),
     *               CORDSP(IX(I,IEL),2))
                  XYZMAX(3,IBLK) = MAX(XYZMAX(3,IBLK),
     *               CORDSP(IX(I,IEL),3))
   50          CONTINUE
   60       CONTINUE
   70    CONTINUE
      END IF
      DO 90 I=1,NDIM
         OVMIN(I) = XYZMIN(I,1)
         OVMAX(I) = XYZMAX(I,1)
         DO 80 J=2,NEBLK
            OVMIN(I) = MIN(OVMIN(I), XYZMIN(I,J))
            OVMAX(I) = MAX(OVMAX(I), XYZMAX(I,J))
   80    CONTINUE
   90 CONTINUE

      ENG1 = ENGNOT(TREAD,2)
      DO 120 IO=IOMIN, IOMAX
         IF (EXODUS) WRITE (IO, 140) ENG1
         IF (NDIM .EQ. 2) THEN
            WRITE (IO, 150)
            DO 100 ITMP=1,NEBLK
               I = MAT(6, ITMP)
               IF (MAT(5,I) .NE. 1) GOTO 100
               WRITE (IO, 170) MAT(1,I),(XYZMIN(J,I),J=1,NDIM),
     *            (XYZMAX(J,I),J=1,NDIM),
     *            ((XYZMAX(J,I)-XYZMIN(J,I)),J=1,NDIM)
  100       CONTINUE
         ELSE
            WRITE (IO, 160)
            DO 110 ITMP=1,NEBLK
               I = MAT(6, ITMP)
               IF (MAT(5,I) .NE. 1) GOTO 110
               WRITE (IO, 180) MAT(1,I),(XYZMIN(J,I),J=1,NDIM),
     *            (XYZMAX(J,I),J=1,NDIM),
     *            ((XYZMAX(J,I)-XYZMIN(J,I)),J=1,NDIM)
  110       CONTINUE
         END IF
         IF (NEBLK .GT. 1) THEN
            IF (NDIM .EQ. 2) THEN
               WRITE (IO, 190) (OVMIN(I),I=1,NDIM),(OVMAX(I),I=1,NDIM),
     *              ((OVMAX(I)-OVMIN(I)),I=1,NDIM)
            ELSE
               WRITE (IO, 200) (OVMIN(I),I=1,NDIM),(OVMAX(I),I=1,NDIM),
     *              ((OVMAX(I)-OVMIN(I)),I=1,NDIM)
            END IF
         END IF
  120 CONTINUE
      IF (EXODUS) GO TO 10
  130 CONTINUE
      RETURN

  140 FORMAT(//5X,'Time = ',A16)
  150 FORMAT(/5X,'Material',T22,'X',T36,'Y')
  160 FORMAT(/5X,'Material',T22,'X',T36,'Y',T50,'Z')
  170 FORMAT (I10,3X,2(2X,1PE12.5),'  Minimum',/
     *   ,13X,2(2X,1PE12.5),'  Maximum',/
     *   ,13X,2(2X,1PE12.5),'  Range',/)
  180 FORMAT (I10,3X,3(2X,1PE12.5),'  Minimum',/
     *   ,13X,3(2X,1PE12.5),'  Maximum',/
     *   ,13X,3(2X,1PE12.5),'  Range',/)
  190 FORMAT (5X,'Limits for Total Body:',/
     *   ,13X,2(2X,1PE12.5),'  Minimum',/
     *   ,13X,2(2X,1PE12.5),'  Maximum',/
     *   ,13X,2(2X,1PE12.5),'  Range',/)
  200 FORMAT (5X,'Limits for Total Body:',/
     *   ,13X,3(2X,1PE12.5),'  Minimum',/
     *   ,13X,3(2X,1PE12.5),'  Maximum',/
     *   ,13X,3(2X,1PE12.5),'  Range',/)
      END

