C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE SUMNOD (CRD, DISP, SVAR, NDIM, NUMNP, INDX,
     *   NODSEL, NAME, TIME, ITMSEL, AXI, AVER, DOABS)

      include 'exodusII.inc'

      REAL CRD(NUMNP,*), DISP(NUMNP,*), SVAR(*)
      CHARACTER*(*) NAME
      REAL TIME(*)
      LOGICAL NODSEL(*), ITMSEL(*), AXI, AVER, DOABS
      CHARACTER*64  STRA
      CHARACTER*132 STRTMP
      CHARACTER*16  ENGNOT, ENG1, ENG2, ENG3, ENG4

      include 'nu_io.blk'
      include 'nu_ptim.blk'
      include 'nu_ndisp.blk'

      PI = ATAN2(0.0, -1.0)
      NLAST = 0

      IF (AVER) THEN
         IF (DOABS) THEN
            STRA = 'Average absolute value of'
         ELSE
            STRA = 'Average'
         END IF
         RSEL = DBLE(NUMEQL (.TRUE., NUMNP, NODSEL))
      ELSE
         IF (DOABS) THEN
            STRA = 'Absolute value total of'
         ELSE
            STRA = 'Total'
         END IF
         RSEL = 1.0
      END IF
      DO 10 IO=IOMIN, IOMAX
         IF (NDIM .EQ. 2 .AND. AXI) THEN
            WRITE (IO,20) STRA(:LENSTR(STRA)), NAME(:LENSTR(NAME))
         ELSE
            WRITE (IO,30) STRA(:LENSTR(STRA)), NAME(:LENSTR(NAME))
         END IF
         WRITE (IO, 40) NAME
   10 CONTINUE
   20 FORMAT (/4X,A,' ',A,' on selected nodes ',
     *   '(2 PI Radius Multiplier)')
   30 FORMAT (/4X,A,' ',A,' on selected nodes ')
   40 FORMAT (/,
     *   4X,'Time      ',A8,/
     *   4X,'----       ------ ')

      RMIN =  1.0E38
      RMAX = -1.0E38

   50 CONTINUE
      NLAST = NLAST + 1
      IF (NLAST .GT. LSTSEL) THEN
         NLAST = 0
         GO TO 120
      ELSE IF (ITMSEL(NLAST)) THEN

C     ... READ THE STEP AND STORE X-DISPLACEMENTS (2D AXISYMMETRIC)

         IF (NDIM .EQ. 2 .AND. AXI) THEN
           call exgnv(ndb, nlast, ndisp(1), numnp, disp(1,1), ierr)
         END IF
         call exgnv(ndb, nlast, indx, numnp, svar, ierr)
         TREAD = TIME(NLAST)

         SUM = 0.0
         IF (NDIM .EQ. 2 .AND. AXI) THEN
            IF (DOABS) THEN
               DO 60 I=1, NUMNP
                  IF (NODSEL(I))
     *               SUM = SUM + 2. * PI * (CRD(I,1) + DISP(I,1)) *
     $               ABS( SVAR(I) )
   60          CONTINUE
            ELSE
               DO 70 I=1, NUMNP
                  IF (NODSEL(I)) SUM = SUM +
     $                 2. * PI * (CRD(I,1) + DISP(I,1)) * SVAR(I)
   70          CONTINUE
            END IF
         ELSE
            IF (DOABS) THEN
               DO 80 I=1, NUMNP
                  IF (NODSEL(I)) SUM = SUM + ABS(SVAR(I))
   80          CONTINUE
            ELSE
               DO 90 I=1, NUMNP
                  IF (NODSEL(I)) SUM = SUM + SVAR(I)
   90          CONTINUE
            END IF
         END IF

         SUM = SUM / RSEL
         IF (SUM .GT. RMAX) THEN
            RMAX = SUM
            ITMX = NLAST
         END IF
         IF (SUM .LT. RMIN) THEN
            RMIN = SUM
            ITMN = NLAST
         END IF

         DO 100 IO=IOMIN,IOMAX
            WRITE (IO, 110) TREAD, SUM
  100    CONTINUE
  110    FORMAT (1X,2(1PE15.8,2X))

      END IF
      GO TO 50

  120 CONTINUE
      eng1 = engnot(rmin,2)
      eng2 = engnot(time(itmn), 2)
      eng3 = engnot(rmax,2)
      eng4 = engnot(time(itmx), 2)
      WRITE (STRTMP, 140) ENG1, ENG2, ENG3, ENG4
      CALL SQZSTR(STRTMP, LSTR)
      DO 130 IO=IOMIN, IOMAX
         WRITE (IO, 150) STRTMP(:LSTR)
  130 CONTINUE
      RETURN

  140 FORMAT ('Minimum = ',A16,' at time ',A16,', maximum = ',A16,
     *   ' at time ',A16,'.')
  150 FORMAT (/,1X,A)

      CALL PRTERR ('PROGRAM',
     *   'Internal code error, contact sponsor')
      STOP 'SUMNOD'
      END
