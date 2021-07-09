C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE GETDSP (CRD, CRDSP, NDIM, NUMNP, TIMES, ITMSEL,
     *   ACTION, ISTAT)

C     ISTAT  -  0 IF OK
C            - -1 IF END OF FILE
C     ACTION -  'R' - REWIND AND POSITION AT BEGINNING OF TIME RECORD
C               'A' - READ NEXT SELECTED STEP, ADD DISP AND COORD
C               'S' - READ NEXT SELECTED STEP, RETURN DISP SEPARATELY
C CALLS:
C     REPOS  - REPOSITION AT BEGINNING OF TIME RECORDS
C     ADDDSP - ADD DISPLACEMENTS TO COORDINATES
C     SKIP   - SKIP A TIMESTEP

      REAL    CRD(NUMNP,*), CRDSP(NUMNP,*), TIMES(*)
      LOGICAL ITMSEL(*)
      CHARACTER*(*) ACTION

      include 'exodusII.inc'
      include 'nu_ptim.blk'
      include 'nu_logs.blk'
      include 'nu_io.blk'
      include 'nu_ndisp.blk'

      IF (NDISP(1) .LE. 0) THEN
         CALL DBVIX ('N', 1, NDISP(1))
         CALL DBVIX ('N', 2, NDISP(2))
         IF (NDIM .GT. 2) CALL DBVIX ('N', 3, NDISP(3))
      END IF

      ISTAT = 0
      IF (.NOT. EXODUS) RETURN

      IF (ACTION(:1) .EQ. 'R') THEN
         NLAST = 0
         RETURN
      ELSE IF (ACTION(:1) .EQ. 'S' .OR. ACTION(:1) .EQ. 'A') THEN
   10    CONTINUE
         NLAST = NLAST + 1
         IF (NLAST .GT. LSTSEL) THEN
            NLAST = 0
            ISTAT = -1
            RETURN
         ELSE IF (ITMSEL(NLAST)) THEN

C ... READ THE STEP AND STORE DISPLACEMENTS

           call exgnv(ndb, nlast, ndisp(1), numnp, crdsp(1,1), ierr)
           call exgnv(ndb, nlast, ndisp(2), numnp, crdsp(1,2), ierr)
           if (ndim .eq. 3) then
             call exgnv(ndb, nlast, ndisp(3), numnp, crdsp(1,3), ierr)
           end if

            TREAD = TIMES(NLAST)

            IF (ACTION(:1) .EQ. 'A') CALL ADDDSP (CRD, CRDSP)
            RETURN
         END IF
         GO TO 10
      END IF

      CALL PRTERR ('PROGRAM',
     *   'Internal code error, contact sponsor')
      STOP 'GETDSP'
      END
