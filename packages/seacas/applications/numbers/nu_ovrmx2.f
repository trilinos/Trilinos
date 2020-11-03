C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE OVRMX2 (LSTEL, CORD, IX, NSEG, MINMAX, NIQSLV,
     *     NIQS, TEMP, LTNESS, NUMIN, NUMFAC, NUMON,
     *     NUMEL, LFACE, NUMNP)

      INTEGER   LSTEL(*), IX(4,*), NIQSLV(*), LTNESS(2,*)
      INTEGER   LFACE(4,*)
      REAL      MINMAX(4,*), CORD(NUMNP,*), TEMP(*)

      CHARACTER*16 ENGNOT, ENG1
      DIMENSION MAP(2,4), V(4), FCORD(2,2), DIST(4), SCORD(2)
      LOGICAL   INSIDE, ONFACE, INIT
      PARAMETER (MAXFAC = 4)
      include 'nu_io.blk'

      DATA MAP /1, 2,  2, 3,  3, 4,  4, 1/

      INIT  = .FALSE.
      NUMIN = 0
      NUMON = 0
      NUMFAC = 0

      DO 10 I=1,NSEG
         IEL = LSTEL(I)

         MINMAX(1, I) = MIN( CORD(IX(1,IEL),1),  CORD(IX(2,IEL),1),
     *        CORD(IX(3,IEL),1),  CORD(IX(4,IEL),1))
         MINMAX(2, I) = MAX( CORD(IX(1,IEL),1),  CORD(IX(2,IEL),1),
     *        CORD(IX(3,IEL),1),  CORD(IX(4,IEL),1))

         MINMAX(3, I) = MIN( CORD(IX(1,IEL),2),  CORD(IX(2,IEL),2),
     *        CORD(IX(3,IEL),2),  CORD(IX(4,IEL),2))
         MINMAX(4, I) = MAX( CORD(IX(1,IEL),2),  CORD(IX(2,IEL),2),
     *        CORD(IX(3,IEL),2),  CORD(IX(4,IEL),2))

 10   CONTINUE

C     ... DETERMINE WHICH FACES HAVE SSET FLAG

      CALL INIINT (MAXFAC * NUMEL, 0, LFACE)

      DO 30 ISEG = 1, NSEG
         IEL = LSTEL(ISEG)
         IFAC1 = LTNESS(1,ISEG)
         IFAC2 = LTNESS(2,ISEG)

         DO 20 IFAC = 1, MAXFAC
            INOD1 = IX(MAP(1,IFAC),IEL)
            INOD2 = IX(MAP(2,IFAC),IEL)

            ITST1 = ISIGN(1,(INOD1-IFAC1)) + ISIGN(1,(IFAC1-INOD1)) +
     *           ISIGN(1,(INOD2-IFAC1)) + ISIGN(1,(IFAC1-INOD2))

            ITST2 = ISIGN(1,(INOD1-IFAC2)) + ISIGN(1,(IFAC2-INOD1)) +
     *           ISIGN(1,(INOD2-IFAC2)) + ISIGN(1,(IFAC2-INOD2))

            LFACE(IFAC,IEL) = LFACE(IFAC,IEL) + ITST1 * ITST2
 20      CONTINUE
 30   CONTINUE

C     ... DETERMINE IF NODE IS CLOSE TO ELEMENT
C     TEMP = 1.0 IF INSIDE MIN/MAX BOX

      DO 130 I=1, NSEG
         IEL = LSTEL(I)

         DO 40 ISLV = 1, NIQS
            ISN = NIQSLV(ISLV)
            TEMP(ISLV) =
     *           (0.5 + SIGN( 0.5,  CORD (ISN,1) - MINMAX(1,I) )) *
     *           (0.5 + SIGN( 0.5, -CORD (ISN,1) + MINMAX(2,I) )) *
     *           (0.5 + SIGN( 0.5,  CORD (ISN,2) - MINMAX(3,I) )) *
     *           (0.5 + SIGN( 0.5, -CORD (ISN,2) + MINMAX(4,I) ))
 40      CONTINUE

C     ... DETERMINE IF ANY INSIDE BOX ( TEMP = 1.0 )

C     ... FOR EACH NODE INSIDE BOX, DETERMINE IF ACTUALLY INSIDE ELEMENT

         DO 120 ISLV = 1, NIQS
            IF (TEMP(ISLV) .EQ. 1.0) THEN
               INOD = NIQSLV(ISLV)

               X3 = CORD(INOD,1)
               Y3 = CORD(INOD,2)

               INSIDE = .TRUE.
               ONFACE = .FALSE.
               DO 50 IPYR = 1, 4

                  X1 = CORD(IX(MAP(1,IPYR),IEL),1)
                  Y1 = CORD(IX(MAP(1,IPYR),IEL),2)

                  X2 = CORD(IX(MAP(2,IPYR),IEL),1)
                  Y2 = CORD(IX(MAP(2,IPYR),IEL),2)

C     ... CALCULATE TRIANGLE AREAS (SHOULD BE DIVIDED BY 2 FOR AREA)

                  V(IPYR) = X1 * (Y2 - Y3) + X2 *  (Y3 - Y1)
     *                 + X3 * (Y1 - Y2)

                  IF (V(IPYR) .LT. 0.0) INSIDE = .FALSE.
                  IF (V(IPYR) .EQ. 0.0) ONFACE = .TRUE.
 50            CONTINUE

C     ... FLAG NODE AND ELEMENT IF INSIDE

               IF (ONFACE .AND. INSIDE) THEN
                  INSIDE = .TRUE.
                  ONFACE = .FALSE.
                  DO 60 IFAC = 1, MAXFAC
                     IF (V(IFAC) .EQ. 0.0 .AND.
     *                    LFACE(IFAC,IEL) .NE. 0 ) THEN
                        INSIDE = .FALSE.
                        ONFACE = .TRUE.
                     END IF
 60               CONTINUE
               END IF

C     ... CHECK FOR NODE ON BOTH SURFACES

               IF (INSIDE) THEN
                  DO 70 INOD = 1, 4
                     IF (IX(INOD,IEL) .EQ. NIQSLV(ISLV)) THEN
                        INSIDE = .FALSE.
                        NUMON = NUMON + 1
                     END IF
 70               CONTINUE
               END IF

               IF (INSIDE) THEN
                  IF (.NOT. INIT) THEN
                     INIT = .TRUE.
                     DO 80 IO=IOMIN, IOMAX
                        WRITE (IO, 150)
 80                  CONTINUE
                  END IF
                  DO 90 IFAC = 1, MAXFAC
                     IF (LFACE(IFAC,IEL) .NE. 0 ) THEN
                        FCORD(1,1) = CORD(IX(MAP(1,IFAC),IEL),1)
                        FCORD(2,1) = CORD(IX(MAP(1,IFAC),IEL),2)
                        FCORD(1,2) = CORD(IX(MAP(2,IFAC),IEL),1)
                        FCORD(2,2) = CORD(IX(MAP(2,IFAC),IEL),2)
                        SCORD(1)   = CORD(NIQSLV(ISLV),1)
                        SCORD(2)   = CORD(NIQSLV(ISLV),2)
                        CALL PENDIS (SCORD, FCORD, DIST(IFAC), 2, 2)
                     END IF
 90               CONTINUE
                  DO 100 IFAC = 1, MAXFAC
                     IF (LFACE(IFAC,IEL) .NE. 0 .AND.
     *                    DIST(IFAC) .NE. 0.0) THEN
                        NUMIN = NUMIN + 1
                        ENG1 = ENGNOT(DIST(IFAC),2)
                        DO 110 IO=IOMIN, IOMAX
                           WRITE (IO,140) NIQSLV(ISLV), IEL,
     *                          ENG1,
     *                          IFAC, CORD(NIQSLV(ISLV),1),
     $                          CORD(NIQSLV(ISLV),2)
 110                    CONTINUE
                     END IF
 100              CONTINUE
               ELSE IF (ONFACE) THEN
                  NUMFAC = NUMFAC + 1
               END IF
            END IF
 120     CONTINUE
 130  CONTINUE
 140  FORMAT (T3,I6,T11,I6,T18,A16,T37,I1,T43,2(F15.8,2X))
 150  FORMAT ('  Slave    Master  Penetration  Element'/
     *     '   Node   Element    Distance     Face      Location')
      RETURN
      END
