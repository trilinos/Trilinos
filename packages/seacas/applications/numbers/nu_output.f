C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE OUTPUT (MASS, DENS, VOLM, CG, ZI, MAT, NDIM, NBLK,
     *   VOL, VOLMN, IELM, NQUAD, LABMAT, AXI, TIME)

      include 'nu_io.blk'
      CHARACTER*16 LABMAT(*)
      CHARACTER*16 ENGNOT, ENG1
      REAL MASS(*), DENS(*), VOLM(*), CG(*), ZI(*)
      DIMENSION MAT(6,*), VOLMN(4,*), IELM(4,*)
      LOGICAL AXI, FIRST
      CHARACTER*6 LABEL(3)
      DATA LABEL/'      ',' Area ','Volume'/
      DATA FIRST /.TRUE./

      DO 160 IO=IOMIN,IOMAX
         ENG1 = ENGNOT(TIME,2)
         IF (FIRST) THEN
            WRITE (IO, 10) NQUAD, ENG1
   10       FORMAT (5X,'Mass Properties Calculation',/
     *              5X,I1,'-Point Quadrature, Time = ',A16)
         ELSE
            WRITE (IO, 20) ENG1
   20       FORMAT ('      Time = ',A16)
         END IF
         TMASS = 0.
         ILAB = 3
         IF (NDIM .EQ. 2 .AND. .NOT. AXI) ILAB = 2
         WRITE (IO, 30) LABEL(ILAB)
   30    FORMAT(/5X,'Material',7X,'Density',14X,A6,16X,'Mass',
     *      13X,'Label')
         DO 50 ITMP=1,NBLK
            I = MAT(6, ITMP)
            IF (MAT(5,I) .NE. 1) GOTO 50
            WRITE (IO, 40) MAT(1,I),DENS(I),VOLM(I),MASS(I),LABMAT(I)
            TMASS = TMASS + MASS(I)
   40       FORMAT (I10,3(5X,1PE15.8),5X,A16)
   50    CONTINUE

         WRITE (IO, 60) VOL, TMASS
   60    FORMAT (25X,2(5X,'----------')/,18X,'Total: ',2(5X,1PE15.8)/)
         IF (NDIM .EQ. 2) THEN
            WRITE (IO, 70) CG(1), CG(2), 0.0
         ELSE
            WRITE (IO, 70) (CG(I),I=1,3)
   70       FORMAT (5x,'MASS PROPERTIES--- (Inertias at centroid)'/
     *         5X,' Xc = ',1PE15.8,'   Yc = ',1PE15.8,
     *         '   Zc = ',1PE15.8)
         END IF
         WRITE (IO, 80) (ZI(I),I=1,3)
   80    FORMAT (5X,'Ixx = ',1PE15.8,'  Iyy = ',1PE15.8,
     *      '  Izz = ',1PE15.8)
         IF (NDIM .EQ. 2 .AND. .NOT. AXI) THEN
            WRITE (IO, 90) ZI(4)
   90       FORMAT (5X,'Ixy = ',1PE15.8/)
         ELSE IF (NDIM .EQ. 3) THEN
            WRITE (IO, 100) (ZI(I),I=4,6)
  100       FORMAT (5X,'Ixy = ',1PE15.8,'  Ixz = ',1PE15.8,
     *         '  Iyz = ',1PE15.8/)
         ELSE
            WRITE (IO, 110)
  110       FORMAT (//)
         END IF
         WRITE (IO, 120)LABEL(NDIM),LABEL(NDIM),LABEL(NDIM)
  120    FORMAT (5X,'Mat',T11,'Minimum',T22,'Element',T32,'Maximum',
     *     T43,'Element',T53,'Average',T63,'Number of',/
     *      5X,'Num',T12,A6,T23,'Number',T33,A6,T44,'Number',T54,A6,
     *     T63,'Elements')
         DO 140 ITMP=1,NBLK
            I = MAT(6,ITMP)
            IF (MAT(5,I) .NE. 1) GOTO 140
            WRITE (IO, 130) MAT(1,I),VOLMN(1,I),IELM(1,I),VOLMN(2,I),
     *         IELM(2,I), VOLMN(3,I)/IELM(3,I), IELM(3,I)
  130       FORMAT (I7,4(2X,1PE9.2,1X,I9))
  140    CONTINUE
  160 CONTINUE
      FIRST = .FALSE.
      RETURN
      END

