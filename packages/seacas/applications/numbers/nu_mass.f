C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE MASSPR (A, TIME, ITMSEL, DENS, MAT, DISP,
     *   NQUAD, LABEL)

      DIMENSION A(*), TIME(*), DENS(*), MAT(6,*),
     *   DISP(NUMNP,*)
      LOGICAL ITMSEL(*), ISABRT
      CHARACTER*16  LABEL(32)
      include 'nu_ptim.blk'
      include 'nu_numg.blk'
      include 'nu_mass.blk'
      include 'nu_logs.blk'

      DIMENSION XI2(2,4), XI3(3,8)
      LOGICAL FIRST, HAVDEN
      DATA FIRST / .TRUE. /
      DATA XI2/ -1.,-1.,  1.,-1.,  1.,1.,  -1.,1./
      DATA XI3/ 1.,-1.,-1.,  -1.,-1.,-1.,  -1.,-1.,1.,  1.,-1.,1.,
     *   1.,1.,-1.,   -1.,1.,-1.,   -1.,1.,1.,   1.,1.,1./

      save

      IF (FIRST) THEN
         FIRST = .FALSE.
         CALL MDRSRV ('MASS'  , IS, NELBLK)
         CALL MDRSRV ('VOLUME', IV, NELBLK)
         CALL MDRSRV ('CENTER', IC, 3)
         CALL MDRSRV ('INERTA', IZ, 6)
         NNODES = 2**NDIM
         NQMAX  = 2**NDIM
         CALL MDRSRV ('XXX'   , IXXX,  (NDIM+1)*NNODES*NQMAX)
         CALL MDRSRV ('XG'    , IXG,   NDIM*NQMAX)
         CALL MDRSRV ('XINI'  , IXINI, NDIM)
C ... 'JACOB' conflicts with jacob in command.f, renamed to jacob1
         CALL MDRSRV ('JACOB1', IAJ,   NDIM*NDIM)
         CALL MDRSRV ('VOL'   , IVM,   4*NELBLK)
         CALL MDRSRV ('IEL'   , IEM,   4*NELBLK)
         CALL MDSTAT (NERRS, NUSED)
         IF (NERRS .GT. 0) THEN
            CALL MEMERR
            STOP
         END IF
      END IF

      HAVDEN = .FALSE.
      DO 20 I=1,NELBLK
         IF (DENS(I) .NE. 0.0) HAVDEN = .TRUE.
   20 CONTINUE

      IF (.NOT. HAVDEN) CALL GETDEN (MAT, DENS, NELBLK, LABEL)

      IF (EXODUS .AND. ISDIS) THEN
         CALL GETDSP (A(IR), DISP, NDIM, NUMNP, TIME, ITMSEL,
     *      'R', ISTAT)
         IF (ISTAT .NE. 0) GO TO 40

   30    CONTINUE
         IF (ISABRT()) RETURN
         CALL GETDSP (A(IR), DISP, NDIM, NUMNP, TIME, ITMSEL,
     *      'A', ISTAT)
         IF (ISTAT .NE. 0) GO TO 40
         IF (NDIM .EQ. 2) THEN
            CALL CGCAL2 (DISP,A(IX),MAT,A(IS),VOL,A(ID),
     *         A(IV),A(IC),A(IZ),A(IXXX),A(IXG),XI2,
     *         A(IXINI),A(IAJ),NNODES,NDIM,NQUAD,
     *         A(IVM),A(IEM),NELBLK,AXI,NUMNP)
         ELSE IF (NDIM .EQ. 3) THEN
            CALL CGCAL3 (DISP,A(IX),MAT,A(IS),VOL,A(ID),
     *         A(IV),A(IC),A(IZ),A(IXXX),A(IXG),XI3,
     *         A(IXINI),A(IAJ),NNODES,NDIM,NQUAD,
     *         A(IVM),A(IEM),NELBLK,NUMNP)
         END IF

         CALL OUTPUT (A(IS), A(ID), A(IV), A(IC), A(IZ), MAT,
     *      NDIM,NELBLK, VOL, A(IVM), A(IEM),
     *      NQUAD, LABEL, AXI, TREAD)

         GO TO 30
   40    CONTINUE
      ELSE
         IF (NDIM .EQ. 2) THEN
            CALL CGCAL2 (A(IR),A(IX),MAT,A(IS),VOL,A(ID),
     *         A(IV),A(IC),A(IZ),A(IXXX),A(IXG),XI2,
     *         A(IXINI),A(IAJ),NNODES,NDIM,NQUAD,
     *         A(IVM),A(IEM),NELBLK,AXI,NUMNP)
         ELSE
            CALL CGCAL3 (A(IR),A(IX),MAT,A(IS),VOL,A(ID),
     *         A(IV),A(IC),A(IZ),A(IXXX),A(IXG),XI3,
     *         A(IXINI),A(IAJ),NNODES,NDIM,NQUAD,
     *         A(IVM),A(IEM),NELBLK,NUMNP)
         END IF
         CALL OUTPUT (A(IS), A(ID), A(IV), A(IC), A(IZ), MAT,
     *      NDIM,NELBLK, VOL, A(IVM), A(IEM), NQUAD, LABEL,
     *      AXI, TREAD)

      END IF
      RETURN
      END
