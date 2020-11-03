C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details
C=======================================================================
      SUBROUTINE PRSTEP (OPTION, NOUT, TIME, NCSTEP, NSTEPS)
C=======================================================================

C   --*** PRSTEP *** (EXPLORE) Display current database time and step number
C   --
C   --PRSTEP displays the time and the step number of the current time step.
C   --
C   --Parameters:
C   --   OPTION - IN - '*' to print all, else print options
C   --   NOUT - IN - the output file, <=0 for standard
C   --   TIME - IN - the current time step time
C   --   NCSTEP - IN - the current time step number
C   --   NSTEPS - IN - the number of time steps

      CHARACTER*(*) OPTION

      CHARACTER*32 STRA, STRB
      CHARACTER*32 RSTR

      INTEGER NPREC
      INTEGER GETPRC

      NPREC=GETPRC()

      STRA = 'time step'
      LSTRA = LENSTR (STRA)
      CALL NUMSTR (1, NPREC, TIME, RSTR, LRSTR)
      WRITE (STRB, 10000, IOSTAT=IDUM) NCSTEP, NSTEPS
10000  FORMAT (I12, ' of ', I12)
      CALL SQZSTR (STRB, LSTRB)
      IF (NOUT .GT. 0) THEN
         WRITE (NOUT, *)
         WRITE (NOUT, 10010)
     &      RSTR(:LRSTR), STRA(:LSTRA), STRB(:LSTRB)
      ELSE
         WRITE (*, 10010)
     &      RSTR(:LRSTR), STRA(:LSTRA), STRB(:LSTRB)
      END IF

      RETURN
10010  FORMAT (1X, 'Time = ', A, ' (', A, ' ', A, ')')
      END
