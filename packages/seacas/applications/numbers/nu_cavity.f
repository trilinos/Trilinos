C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE CAVITY (A, CRD, IDESS, NEESS, NNESS, IPEESS, IPNESS,
     *   LTEESS, LTNESS, FACESS, DISP, NUMNP, NDIM, NUMESS,
     *   TIME, ITMSEL, TITLE, CENT, CENTER)

      include 'nu_io.blk'
      DIMENSION A(*), CRD(NUMNP,NDIM), IDESS(*), NEESS(*),
     *   NNESS(*), IPEESS(*), IPNESS(*), LTEESS(*), LTNESS(*),
     *   FACESS(*), TIME(*), DISP(NUMNP,NDIM), CENT(3)
      LOGICAL ITMSEL(*)
      CHARACTER*80 TITLE
      include 'nu_logs.blk'
      include 'nu_ptim.blk'
      include 'nu_cav.blk'
      LOGICAL ERROR, CENTER

      CALL GETCAV (ERROR, IDESS, NUMESS)
      IF (ERROR) RETURN

      TVOL = 0.0
      DO 10 NCAV = 1, NUMCAV
         IFLG = IFND(NCAV)
         IPTR = IPNESS(IFLG)
         IF (NDIM .EQ. 3) THEN
            CALL VOL3D( CRD, LTNESS(IPTR), NEESS(IFLG), VOLUME,
     *         NDIM, NUMESS, CENT, NUMNP, CENTER)
         ELSE
            CALL VOL2D( CRD, LTNESS(IPTR), NEESS(IFLG), VOLUME,
     *         NDIM, NUMESS, AXI, CENT, NUMNP, CENTER)
         END IF

         TVOL = TVOL + VOLUME
   10 CONTINUE
      DO 20 IO=IOMIN, IOMAX
         WRITE (IO,30) (ICAV(I),I=1,NUMCAV)
         IF (NDIM .EQ. 2) THEN
            WRITE (IO, 40) CENT(1),CENT(2)
         ELSE
            WRITE (IO, 50) CENT(1),CENT(2),CENT(3)
         END IF
         WRITE (IO,60) TVOL
   20 CONTINUE
   30 FORMAT (/' Cavity Flag(s): ',8I8)
   40 FORMAT ( ' Apex at X =',1PE15.8,', Y =',1PE15.8)
   50 FORMAT ( ' Apex at X =',1PE15.8,', Y =',1PE15.8,', Z =',1PE15.8)
   60 FORMAT (/' Undeformed Volume of Cavity is ',1PE15.8)

C ... REWIND EXODUS FILE TO BEGINNING OF TIMESTEPS

      IF (EXODUS .AND. ISDIS) THEN
         TIMEL = STMIN
         CALL GETDSP (CRD, DISP, NDIM, NUMNP, TIME, ITMSEL, 'R', ISTAT)
         IF (ISTAT .NE. 0) GO TO 140
         DO 70 IO=IOMIN, IOMAX
            WRITE (IO, 80)
   70    CONTINUE
   80    FORMAT (/,
     *      4X,'                 Cavity           Total',
     *     '            Timestep         Rate of',/
     *      4X,'Time             Volume           Change',
     *     '           Change           Change',/
     *      4X,'----             ------           ------',
     *     '           --------         -------')

         DELLAS = 0.0
   90    CONTINUE
         CALL GETDSP (CRD, DISP, NDIM, NUMNP, TIME, ITMSEL, 'S', ISTAT)
         IF (ISTAT .NE. 0) GO TO 140
         DELVOL = 0.0
         DO 100 NCAV = 1, NUMCAV
            IFLG = IFND(NCAV)
            IPTR = IPNESS(IFLG)

C     NOTE: Positive delcav = shrink in cavity volume

            IF (NDIM .EQ. 3) THEN
               CALL DVOL3D(CRD, DISP, LTNESS(IPTR),
     *            NEESS(IFLG), DELCAV, NDIM, NUMNP)
            ELSE
               CALL DVOL2D(CRD, DISP, LTNESS(IPTR),
     *            NEESS(IFLG), DELCAV, NDIM, AXI, NUMNP)
            END IF

            DELVOL =  DELVOL + DELCAV
  100    CONTINUE
         DELDEL = DELVOL - DELLAS
         IF (TREAD .EQ. TIMEL) THEN
            DO 110 IO=IOMIN, IOMAX
               WRITE (IO, 130) TREAD, TVOL-DELVOL, -DELVOL,
     *            -DELDEL
  110       CONTINUE
         ELSE
            DELRAT = DELDEL / (TREAD - TIMEL)
            DO 120 IO=IOMIN, IOMAX
               WRITE (IO, 130) TREAD, TVOL-DELVOL, -DELVOL,
     *            -DELDEL, -DELRAT
  120       CONTINUE
         END IF
  130    FORMAT (1X,5(1PE15.8,2X))
         DELLAS = DELVOL
         TIMEL  = TREAD
         GO TO 90
      END IF
  140 CONTINUE
      CALL INIINT (NUMCAV, 0, ICAV)
      NUMCAV = 0
      RETURN
      END
