C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE OVRMX3 (LSTEL, CORD, IX, NSEG, MINMAX, NIQSLV,
     *   NIQS, TEMP, LTNESS, NUMIN, NUMFAC, NUMON,
     *   NUMEL, LFACE, NUMNP)

      INTEGER   LSTEL(*), IX(8,*), NIQSLV(*), LTNESS(4,*)
      INTEGER   LFACE(6,*)
      REAL      MINMAX(6,*), CORD(NUMNP,*), TEMP(*)

      CHARACTER*16 ENGNOT, ENG1
      DIMENSION MAP(4,6), V(6), FCORD(3,4), DIST(6), SCORD(3)
      LOGICAL   INSIDE, ONFACE, INIT
      PARAMETER (MAXFAC = 6)
      include 'nu_io.blk'

      DATA MAP /1, 2, 3, 4,   6, 7, 3, 2,  6, 5, 8, 7,
     *   5, 1, 4, 8,   4, 3, 7, 8,  1, 5, 6, 2/

      INIT = .FALSE.
      NUMIN = 0
      NUMON = 0
      NUMFAC = 0

      DO 10 I=1,NSEG
         IEL = LSTEL(I)

         MINMAX(1, I) =     MIN(CORD(IX(1,IEL),1),  CORD(IX(2,IEL),1),
     *      CORD(IX(3,IEL),1),  CORD(IX(4,IEL),1),  CORD(IX(5,IEL),1),
     *      CORD(IX(6,IEL),1),  CORD(IX(7,IEL),1),  CORD(IX(8,IEL),1))
         MINMAX(2, I) =     MAX(CORD(IX(1,IEL),1),  CORD(IX(2,IEL),1),
     *      CORD(IX(3,IEL),1),  CORD(IX(4,IEL),1),  CORD(IX(5,IEL),1),
     *      CORD(IX(6,IEL),1),  CORD(IX(7,IEL),1),  CORD(IX(8,IEL),1))

         MINMAX(3, I) =     MIN(CORD(IX(1,IEL),2),  CORD(IX(2,IEL),2),
     *      CORD(IX(3,IEL),2),  CORD(IX(4,IEL),2),  CORD(IX(5,IEL),2),
     *      CORD(IX(6,IEL),2),  CORD(IX(7,IEL),2),  CORD(IX(8,IEL),2))
         MINMAX(4, I) =     MAX(CORD(IX(1,IEL),2),  CORD(IX(2,IEL),2),
     *      CORD(IX(3,IEL),2),  CORD(IX(4,IEL),2),  CORD(IX(5,IEL),2),
     *      CORD(IX(6,IEL),2),  CORD(IX(7,IEL),2),  CORD(IX(8,IEL),2))

         MINMAX(5, I) =     MIN(CORD(IX(1,IEL),3),  CORD(IX(2,IEL),3),
     *      CORD(IX(3,IEL),3),  CORD(IX(4,IEL),3),  CORD(IX(5,IEL),3),
     *      CORD(IX(6,IEL),3),  CORD(IX(7,IEL),3),  CORD(IX(8,IEL),3))
         MINMAX(6, I) =     MAX(CORD(IX(1,IEL),3),  CORD(IX(2,IEL),3),
     *      CORD(IX(3,IEL),3),  CORD(IX(4,IEL),3),  CORD(IX(5,IEL),3),
     *      CORD(IX(6,IEL),3),  CORD(IX(7,IEL),3),  CORD(IX(8,IEL),3))

   10 CONTINUE

C  ... DETERMINE WHICH FACES HAVE SSET FLAG

      CALL INIINT (MAXFAC * NUMEL, 0, LFACE)

      DO 30 ISEG = 1, NSEG
         IEL = LSTEL(ISEG)
         IFAC1 = LTNESS(1,ISEG)
         IFAC2 = LTNESS(2,ISEG)
         IFAC3 = LTNESS(3,ISEG)
         IFAC4 = LTNESS(4,ISEG)

         DO 20 IFAC = 1, MAXFAC
            INOD1 = IX(MAP(1,IFAC),IEL)
            INOD2 = IX(MAP(2,IFAC),IEL)
            INOD3 = IX(MAP(3,IFAC),IEL)
            INOD4 = IX(MAP(4,IFAC),IEL)

            ITST1 = ISIGN(1,(INOD1-IFAC1)) + ISIGN(1,(IFAC1-INOD1)) +
     *         ISIGN(1,(INOD2-IFAC1)) + ISIGN(1,(IFAC1-INOD2)) +
     *         ISIGN(1,(INOD3-IFAC1)) + ISIGN(1,(IFAC1-INOD3)) +
     *         ISIGN(1,(INOD4-IFAC1)) + ISIGN(1,(IFAC1-INOD4))

            ITST2 = ISIGN(1,(INOD1-IFAC2)) + ISIGN(1,(IFAC2-INOD1)) +
     *         ISIGN(1,(INOD2-IFAC2)) + ISIGN(1,(IFAC2-INOD2)) +
     *         ISIGN(1,(INOD3-IFAC2)) + ISIGN(1,(IFAC2-INOD3)) +
     *         ISIGN(1,(INOD4-IFAC2)) + ISIGN(1,(IFAC2-INOD4))

            ITST3 = ISIGN(1,(INOD1-IFAC3)) + ISIGN(1,(IFAC3-INOD1)) +
     *         ISIGN(1,(INOD2-IFAC3)) + ISIGN(1,(IFAC3-INOD2)) +
     *         ISIGN(1,(INOD3-IFAC3)) + ISIGN(1,(IFAC3-INOD3)) +
     *         ISIGN(1,(INOD4-IFAC3)) + ISIGN(1,(IFAC3-INOD4))

            ITST4 = ISIGN(1,(INOD1-IFAC4)) + ISIGN(1,(IFAC4-INOD1)) +
     *         ISIGN(1,(INOD2-IFAC4)) + ISIGN(1,(IFAC4-INOD2)) +
     *         ISIGN(1,(INOD3-IFAC4)) + ISIGN(1,(IFAC4-INOD3)) +
     *         ISIGN(1,(INOD4-IFAC4)) + ISIGN(1,(IFAC4-INOD4))

C ... LFACE(IFAC,IEL) = 0 IF FACE NOT ON CONTACT SURFACE
C                     > 0 IF FACE ON CONTACT SURFACE

            LFACE(IFAC,IEL) = LFACE(IFAC,IEL) +
     *         ITST1 * ITST2 * ITST3 * ITST4
   20    CONTINUE
   30 CONTINUE

C ... DETERMINE IF NODE IS CLOSE TO ELEMENT
C     TEMP = 1.0 IF INSIDE MIN/MAX BOX

      DO 150 I=1, NSEG
         IEL = LSTEL(I)

         DO 40 ISLV = 1, NIQS
            ISN = NIQSLV(ISLV)
            TEMP(ISLV) =
     *         (0.5 + SIGN( 0.5,  CORD (ISN,1) - MINMAX(1,I) )) *
     *         (0.5 + SIGN( 0.5, -CORD (ISN,1) + MINMAX(2,I) )) *
     *         (0.5 + SIGN( 0.5,  CORD (ISN,2) - MINMAX(3,I) )) *
     *         (0.5 + SIGN( 0.5, -CORD (ISN,2) + MINMAX(4,I) )) *
     *         (0.5 + SIGN( 0.5,  CORD (ISN,3) - MINMAX(5,I) )) *
     *         (0.5 + SIGN( 0.5, -CORD (ISN,3) + MINMAX(6,I) ))
   40    CONTINUE

C ... DETERMINE IF ANY INSIDE BOX ( TEMP = 1.0 )

C ... FOR EACH NODE INSIDE BOX, DETERMINE IF ACTUALLY INSIDE ELEMENT

         DO 140 ISLV = 1, NIQS
            IF (TEMP(ISLV) .EQ. 1.0) THEN
               INOD = NIQSLV(ISLV)

               X5 = CORD(INOD,1)
               Y5 = CORD(INOD,2)
               Z5 = CORD(INOD,3)

               DO 50 IPYR = 1, MAXFAC

                  X1 = CORD(IX(MAP(1,IPYR),IEL),1)
                  Y1 = CORD(IX(MAP(1,IPYR),IEL),2)
                  Z1 = CORD(IX(MAP(1,IPYR),IEL),3)

                  X2 = CORD(IX(MAP(2,IPYR),IEL),1)
                  Y2 = CORD(IX(MAP(2,IPYR),IEL),2)
                  Z2 = CORD(IX(MAP(2,IPYR),IEL),3)

                  X3 = CORD(IX(MAP(3,IPYR),IEL),1)
                  Y3 = CORD(IX(MAP(3,IPYR),IEL),2)
                  Z3 = CORD(IX(MAP(3,IPYR),IEL),3)

                  X4 = CORD(IX(MAP(4,IPYR),IEL),1)
                  Y4 = CORD(IX(MAP(4,IPYR),IEL),2)
                  Z4 = CORD(IX(MAP(4,IPYR),IEL),3)

                  Z31 = Z3 - Z1
                  Z42 = Z4 - Z2
                  Z51 = Z5 - Z1
                  Z52 = Z5 - Z2
                  Z53 = Z5 - Z3
                  Z54 = Z5 - Z4

C ... CALCULATE PYRAMIDAL VOLUMES (SHOULD BE DIVIDED BY 12 FOR VOLUME)

                  V(IPYR) = ((2.*Y5 - Y3) * Z42 + Y2 * (Z53 + Z54) -
     *               Y4 * (Z53 + Z52) ) * X1 +
     *               ( (Y4 - 2.*Y5) * Z31 + Y3 * (Z54 + Z51) -
     *               Y1 * (Z54 + Z53) ) * X2 +
     *               ( (Y1 - 2.*Y5) * Z42 + Y4 * (Z51 + Z52) -
     *               Y2 * (Z54 + Z51) ) * X3 +
     *               ( (2.*Y5 - Y2) * Z31 + Y1 * (Z52 + Z53) -
     *               Y3 * (Z52 + Z51) ) * X4 +
     *               ( (Y2 - Y4)  * (Z3 - Z1) + (Y3 - Y1) *
     *               (Z4 - Z2) ) * 2. * X5
   50          CONTINUE
               INSIDE = .TRUE.
               ONFACE = .FALSE.
               DO 60 IPYR = 1, MAXFAC
                  IF (V(IPYR) .LT. 0.0) INSIDE = .FALSE.
                  IF (V(IPYR) .EQ. 0.0) ONFACE = .TRUE.
   60          CONTINUE

C ... FLAG NODE AND ELEMENT IF INSIDE

               IF (ONFACE .AND. INSIDE) THEN
                  INSIDE = .TRUE.
                  ONFACE = .FALSE.
                  DO 70 IFAC = 1, MAXFAC
                     IF (V(IFAC) .EQ. 0.0 .AND.
     *                  LFACE(IFAC,IEL) .NE. 0 ) THEN
                        INSIDE = .FALSE.
                        ONFACE = .TRUE.
                     END IF
   70             CONTINUE
               END IF

C ... CHECK FOR NODE ON BOTH SURFACES

               IF (INSIDE) THEN
                  DO 80 INOD = 1, 8
                     IF (IX(INOD,IEL) .EQ. NIQSLV(ISLV)) THEN
                        INSIDE = .FALSE.
                        NUMON = NUMON + 1
                     END IF
   80             CONTINUE
               END IF

               IF (INSIDE) THEN
                  DO 110 IFAC = 1, MAXFAC
                     IF (LFACE(IFAC,IEL) .NE. 0 ) THEN
                        DO 100 NOD = 1, 4
                           FCORD(1,NOD) = CORD(IX(MAP(NOD,IFAC),IEL),1)
                           FCORD(2,NOD) = CORD(IX(MAP(NOD,IFAC),IEL),2)
                           FCORD(3,NOD) = CORD(IX(MAP(NOD,IFAC),IEL),3)
  100                   CONTINUE
                        SCORD(1)   = CORD(NIQSLV(ISLV),1)
                        SCORD(2)   = CORD(NIQSLV(ISLV),2)
                        SCORD(3)   = CORD(NIQSLV(ISLV),3)
                        CALL PENDIS (SCORD, FCORD, DIST(IFAC), 3, 4)
                     END IF
  110             CONTINUE
                     DO 120 IFAC = 1, MAXFAC
                        IF (LFACE(IFAC,IEL) .NE. 0 .AND.
     *                     DIST(IFAC) .NE. 0.0) THEN
                           NUMIN = NUMIN + 1
                           IF (.NOT. INIT) THEN
                              INIT = .TRUE.
                              DO 90 IO2=IOMIN, IOMAX
                                 WRITE (IO2, 170)
   90                         CONTINUE
                           END IF
                           ENG1 = ENGNOT(DIST(IFAC),2)
                           DO 130 IO=IOMIN, IOMAX
                              WRITE (IO,160) NIQSLV(ISLV), IEL,
     *                             ENG1,
     $                             CORD(NIQSLV(ISLV),1),
     $                             CORD(NIQSLV(ISLV),2),
     $                             CORD(NIQSLV(ISLV),3)

c     *                             IFAC, (IX(MAP(J,IFAC),IEL),J=1,4)
 130                    CONTINUE
                        END IF
  120                CONTINUE
               ELSE IF (ONFACE) THEN
                  NUMFAC = NUMFAC + 1
               END IF
            END IF
  140    CONTINUE
  150 CONTINUE
 160  FORMAT (T2,I6,T11,I6,T18,A16,T36,3(1PE15.8,2X))
C  160 FORMAT (T3,I6,T11,I6,T18,A16,T37,I1,T43,4(I6,2X))
  170 FORMAT (/
     &     '  Slave    Master    Penetration        Nodal Coordinates'/,
     &     '   Node   Element      Distance',
     &     '         X           Y           Z')
      RETURN
      END
