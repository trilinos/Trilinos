C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE TPREAD (A, NPTIMS, IPTIMS, TIMES, WHOTIM,
     &                   PLTTIM, PLTVAL)
C=======================================================================

C   --*** TPREAD *** (TPLOT) Read plot variables from database
C   --   Written by Amy Gilkey - revised 02/02/88
C   --
C   --TPREAD reads the database and stores the variables to be plotted.
C   --
C   --This routine manipulates dynamic memory, so check after return.
C   --
C   --Parameters:
C   --   A      - IN  - the dynamic memory base array
C   --   NPTIMS - IN  - the number of selected steps
C   --   IPTIMS - IN  - the selected time steps
C   --   TIMES  - IN  - the time step times
C   --   WHOTIM - IN  - true iff whole (versus history) time step
C   --   PLTTIM - OUT - the plot times (if TIMPLT)
C   --   PLTVAL - OUT - the plot variable values
C   --
C   --Common Variables:
C   --   Uses NDIM, NUMNP, NUMEL, NELBLK,
C   --      NVARHI, NVARGL, NVARNP, NVAREL of /DBNUMS/
C   --   Uses ITVID, TIMPLT of /TPVARS/

      include 'params.blk'
      include 'dbase.blk'
      include 'dbnums.blk'
      include 'tpvars.blk'

      DIMENSION A(*)
      INTEGER IPTIMS(NPTIMS)
      REAL PLTTIM(NPTIMS)
      REAL TIMES(*)
      LOGICAL WHOTIM(*)
      REAL PLTVAL(NPTIMS,NTPVAR)

      CHARACTER TYP

C   --Reserve memory for data record

      CALL MDRSRV ('VALS', KVALS, NSTEPS)
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 130

C   --Transfer element variables onto random file (for efficiency)

      do 140 n = 1, ntpvar
        call dbvtyp_bl (itvid(n), typ, inum)
         IF (TYP .EQ. 'H') THEN
           call prterr ('PROGRAM',
     *       'History variables should not exist')
         ELSE IF (TYP .EQ. 'G') THEN
           call exggvt(ndb, inum, 1, nsteps, a(kvals), ierr)
         ELSE IF (TYP .EQ. 'N') THEN
           call exgnvt(ndb, inum, itvne(n), 1, nsteps, a(kvals), ierr)
         ELSE IF (TYP .EQ. 'E') THEN
           call exgevt(ndb, inum, itvne(n), 1, nsteps, a(kvals), ierr)
         END IF

         DO 120 NPT = 1, NPTIMS
           ISTEP = IPTIMS(NPT)
C          --Read and store variable data to be plotted
           pltval(npt, n) = a(kvals-1+istep)
 120    CONTINUE
 140  CONTINUE

C      --Store time if it is used as a plot variable
      if (timplt) then
        do 150 npt = 1, nptims
           ISTEP = IPTIMS(NPT)
           PLTTIM(NPT) = TIMES(ISTEP)
 150     continue
       end if

      CALL MDDEL ('VALS')

  130 CONTINUE
      RETURN
      END
