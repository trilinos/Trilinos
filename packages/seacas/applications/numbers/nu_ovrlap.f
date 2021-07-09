C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE OVRLAP (A, COORD, IDESS, NEESS, NNESS, IPEESS, IPNESS,
     *   LTEESS, LTNESS, FACESS, DISP, NUMNP, NDIM, NUMESS,
     *   TIME, ITMSEL, TITLE, IMAS, ISLV, NUMEL)

      DIMENSION A(*), COORD(NUMNP,*), IDESS(*), NEESS(*),
     *   NNESS(*), IPEESS(*), IPNESS(*), LTEESS(*), LTNESS(*),
     *   FACESS(*), TIME(*), DISP(NUMNP,*)
      LOGICAL ITMSEL(*)
      CHARACTER*80 TITLE, STRA
      include 'nu_mass.blk'
      include 'nu_logs.blk'
      include 'nu_ptim.blk'
      include 'nu_io.blk'
      LOGICAL ERROR

      DIMENSION CPTIME(10)

      DO 10 I=1,10
         CPTIME(I) = 0.0
   10 CONTINUE

      IFLGM = LOCINT (IMAS, NUMESS, IDESS)
      IFLGS = LOCINT (ISLV, NUMESS, IDESS)

      ERROR = .FALSE.
      IF (IFLGM .EQ. 0) THEN
         WRITE (STRA, 30) 'Master', IMAS
         CALL SQZSTR (STRA, LSTR)
         CALL PRTERR ('ERROR', STRA(:LSTR))
   30    FORMAT (1X,A,' Surface Flag ',I5,' not found. ')
         ERROR = .TRUE.
      END IF
      IF (IFLGS .EQ. 0) THEN
         WRITE (STRA, 30) 'Slave', ISLV
         CALL SQZSTR (STRA, LSTR)
         CALL PRTERR ('ERROR', STRA(:LSTR))
         ERROR = .TRUE.
      END IF
      IF (ERROR) RETURN

      WRITE (STRA, 50) IMAS, ISLV
      CALL SQZSTR (STRA, LSTR)
      DO 40 IO=IOMIN, IOMAX
         WRITE (IO, 55) STRA(:LSTR)
   40 CONTINUE
   50 FORMAT ('Checking Master Surface ',I5,' Versus Slave Surface ',I5)
   55 FORMAT (/,1X,A,/)

      NSEGM = NEESS(IFLGM)
      IPTRM = IPNESS(IFLGM)
      IEPTM = IPEESS(IFLGM)

      NSEGS = NEESS(IFLGS)
      IPTRS = IPNESS(IFLGS)

      MULT = 2 * NDIM - 2

C ... PROCESS SLAVE SET TO REMOVE DUPLICATE NODES

      CALL MDRSRV ('MAPSLV', IMPSL, MULT*NSEGS)
      CALL MDRSRV ('ITEMP',  ITMP,  MAX(NUMNP,3*NSEGM))
      CALL UNIQUE (LTNESS(IPTRS), MULT*NSEGS, A(IMPSL), A(ITMP),
     *   NIQS, NUMNP)
      CALL MDRSRV ('NIQSLV', INQS, NIQS)
      CALL TRANIQ (LTNESS(IPTRS), A(IMPSL), A(INQS), MULT*NSEGS, 1)

      CALL MDRSRV ('MINMAX', IMNMX, 2*NDIM*NSEGM)
      CALL MDRSRV ('LFACE',  ILFAC, 2*NDIM*NUMEL)
      CALL MDSTAT (NERRS, NUSED)
      IF (NERRS .GT. 0) THEN
         CALL MEMERR
         STOP
      END IF

C ... Beginning of Time Step Loop

      IF (EXODUS .AND. ISDIS) THEN
         CALL GETDSP (A(IR), DISP, NDIM, NUMNP, TIME, ITMSEL,
     *        'R', ISTAT)
         IF (ISTAT .NE. 0) GO TO 150
      END IF

 60   CONTINUE
      IF (EXODUS) THEN
         CALL GETDSP (A(IR), DISP, NDIM, NUMNP, TIME, ITMSEL,
     *        'A', ISTAT)
         IF (ISTAT .NE. 0) GO TO 150
         DO 160 IO=IOMIN, IOMAX
           WRITE (IO, 56) TREAD
 160      CONTINUE
 56      FORMAT (' At Time = ', 1PE15.8)
         IF (NDIM .EQ. 3) THEN
            CALL OVRMX3 (LTEESS(IEPTM), DISP, A(IX), NSEGM, A(IMNMX),
     *           A(INQS), NIQS, A(ITMP), LTNESS(IPTRM),
     *           NUMIN, NUMFAC, NUMON, NUMEL, A(ILFAC), NUMNP)
         ELSE
            CALL OVRMX2 (LTEESS(IEPTM), DISP, A(IX), NSEGM, A(IMNMX),
     *           A(INQS), NIQS, A(ITMP), LTNESS(IPTRM),
     *           NUMIN, NUMFAC, NUMON, NUMEL, A(ILFAC), NUMNP)
         END IF

      ELSE
C ... Not EXODUS
         IF (NDIM .EQ. 3) THEN
            CALL OVRMX3 (LTEESS(IEPTM), COORD, A(IX), NSEGM, A(IMNMX),
     *           A(INQS), NIQS, A(ITMP), LTNESS(IPTRM),
     *           NUMIN, NUMFAC, NUMON, NUMEL, A(ILFAC), NUMNP)
         ELSE
            CALL OVRMX2 (LTEESS(IEPTM), COORD, A(IX), NSEGM, A(IMNMX),
     *           A(INQS), NIQS, A(ITMP), LTNESS(IPTRM),
     *           NUMIN, NUMFAC, NUMON, NUMEL, A(ILFAC), NUMNP)
         END IF
      END IF

      IF (NUMIN .GT. 0) THEN
         DO 70 IO=IOMIN, IOMAX
            WRITE (IO, 80) NUMIN
   70    CONTINUE
   80    FORMAT (' Warning -',I5,' Nodes penetrate master surface.')
      ELSE
         DO 90 IO=IOMIN, IOMAX
            WRITE (IO, 100)
   90    CONTINUE
  100    FORMAT (' No nodes found penetrating master surface.')
      END IF
      IF (NUMFAC .GT. 0) THEN
         DO 110 IO=IOMIN, IOMAX
            WRITE (IO, 120) NUMFAC
  110    CONTINUE
  120    FORMAT (' ',I5,' Nodes on an element face.')
      END IF
      IF (NUMON .GT. 0) THEN
         DO 130 IO=IOMIN, IOMAX
            WRITE (IO, 140) NUMON
  130    CONTINUE
  140    FORMAT (' ',I5,' Nodes in both side sets.')
      END IF
      IF (EXODUS) GO TO 60
 150  CONTINUE

      CALL MDDEL ('LFACE')
      CALL MDDEL ('MINMAX')
      CALL MDDEL ('NIQSLV')
      CALL MDDEL ('ITEMP' )
      CALL MDDEL ('MAPSLV')
      RETURN
      END
