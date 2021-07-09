C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE LOCOUT (TYPE, NDIM, NODEL, TOLER, SORT, P1, P2, BOUND)
      DIMENSION P1(NDIM), P2(NDIM), TOLER(2)
      CHARACTER*(*) NODEL, BOUND, SORT, TYPE
      CHARACTER*16  BNAME
      include 'nu_io.blk'

      IF (NODEL(:1) .EQ. 'E') THEN
         NODEL = 'Elements'
      ELSE
         NODEL = 'Nodes'
      ENDIF

      IF (BOUND(:3) .EQ. 'BOU') THEN
         BNAME = 'Bounded Search'
      ELSE
         BNAME = 'Unbounded Search'
      END IF

      DO 10 IO=IOMIN, IOMAX
         WRITE (IO, 20) NODEL(:LENSTR(NODEL)), TOLER(1), TOLER(2),
     *      TYPE(:LENSTR(TYPE))
   20 FORMAT (//' Locating all ',A,' at a distance ',1PE15.8,
     *   ' plus/minus ',1PE15.8,/' from the ',A)

      IF (TYPE .EQ. 'LINE') THEN
         IF (NDIM .EQ. 2) THEN
            WRITE (IO, 30) (P1(I),I=1,NDIM), (P2(I),I=1,NDIM), BNAME
         ELSE
            WRITE (IO, 40) (P1(I),I=1,NDIM), (P2(I),I=1,NDIM), BNAME
         END IF
      ELSE IF (TYPE .EQ. 'POINT') THEN
         IF (NDIM .EQ. 2) THEN
            WRITE (IO, 50) (P1(I), I=1, NDIM)
         ELSE
            WRITE (IO, 60) (P1(I), I=1, NDIM)
         END IF
      ELSE IF (TYPE .EQ. 'PLANE') THEN
      A = P2(1)
      B = P2(2)
      C = P2(3)
      D = A * P1(1) + B * P1(2) + C * P1(3)

            WRITE (IO, 70) A, B, C, D
      END IF
      WRITE (IO, 80) SORT(:LENSTR(SORT))
   10 CONTINUE

   30 FORMAT (' from Point (',2(1PE11.3),')',/
     *        ' to   Point (',2(1PE11.3),')',2X,A)
   40 FORMAT (' from Point (',3(1PE11.3),')',/
     *        ' to   Point (',3(1PE11.3),')',2X,A)
   50 FORMAT (' (',2(1PE11.3),')')
   60 FORMAT (' (',3(1PE11.3),')')
   70 FORMAT (' ',1PE15.8,' X + ',1PE15.8,' Y + ',1PE15.8,
     *    ' Z = ',1PE15.8)
   80 FORMAT (' Sorted on field ',A,/)

      RETURN
      END
