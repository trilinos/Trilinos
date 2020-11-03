C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details
C=======================================================================
      SUBROUTINE PRMMAX (NOUT, MMSTEP, MMNAME, MMTYP, MMVAR, MMNUM,
     &   XMIN, XMAX, NUMELB, IDELB, ISEVOK, VARGL, VARNP, VAREL)
C=======================================================================

C   --*** PRMMAX *** (EXPLORE) Calculate and print the variable min/max
C   --
C   --PRMMAX calculates the minimum and maximum of a variable for either
C   --one step or all steps and displays the values.
C   --
C   --Parameters:
C   --   NOUT - IN - the output file, <=0 for standard
C   --   MMSTEP - IN - the requested step number, <=0 for all
C   --   MMNAME - IN - min/max variable name
C   --   MMTYP - IN/OUT - min/max variable type:
C   --      'H'istory, 'G'lobal, 'N'odal, 'E'lement
C   --   MMVAR - IN - min/max variable number
C   --   MMNUM - IN/OUT - number of sequential min/max requests for this
C   --      variable
C   --   XMIN, XMAX - IN/OUT - the minimum and maximum values,
C   --      input as last values if MMNUM > 1, output as current values
C   --   NUMELB - IN - the number of elements per block
C   --   ISEVOK - IN - the element block variable truth table;
C   --      variable i of block j exists iff ISEVOK(i,j) is NOT 0
C   --   VARGL - IN/OUT - the global variables for current time step
C   --   VARNP - IN/OUT - the nodal variables for current time step
C   --   VAREL - IN/OUT - the element variables for current time step
C   --
C   --Common Variables:
C   --   Uses NUMNP, NUMEL, NVARNP, NVAREL, NVARGL of /DBNUMS/

      include 'exodusII.inc'
      include 'exp_dbnums.blk'

      CHARACTER*(*) MMNAME
      CHARACTER MMTYP
      INTEGER NUMELB(*)
      INTEGER IDELB(*)
      INTEGER ISEVOK(*)
      REAL VARGL(*), VARNP(*), VAREL(*)

      CHARACTER*20 STR

      IF (MMNUM .GT. 1) THEN
         OLDMIN = XMIN
         OLDMAX = XMAX
      END IF
      IF (MMSTEP .LE. 0) THEN
         ISTART = 1
         IEND = NSTEPS
      ELSE
         ISTART = MMSTEP
         IEND = MMSTEP
      END IF

      DO 100 NSTEP = ISTART, IEND

C      --Read in the variable for the time step

         N = NSTEP
         CALL TOSTEP (N, NUMELB, IDELB, ISEVOK,
     &      TIME, VARGL, VARNP, VAREL)
         IF (N .NE. NSTEP) GOTO 110

C      --Update the minimum and maximum

         N = NSTEP
         IF (MMSTEP .GT. 0) N = -999
         IF (MMTYP .EQ. 'G') THEN
            CALL NXMMAX (N, MMNUM, MMVAR, NVARGL, 1, VARGL,
     &         OLDMIN, OLDMAX, XMIN, XMAX,
     &         MINSTE, MAXSTE, MINNE, MAXNE, NMIN, NMAX)
         ELSE IF (MMTYP .EQ. 'N') THEN
            CALL NXMMAX (N, MMNUM, MMVAR, NVARNP, NUMNP, VARNP,
     &         OLDMIN, OLDMAX, XMIN, XMAX,
     &         MINSTE, MAXSTE, MINNE, MAXNE, NMIN, NMAX)
         ELSE IF (MMTYP .EQ. 'E') THEN
            CALL NXMMAX (N, MMNUM, MMVAR, NVAREL, NUMEL, VAREL,
     &         OLDMIN, OLDMAX, XMIN, XMAX,
     &         MINSTE, MAXSTE, MINNE, MAXNE, NMIN, NMAX)
         END IF
  100 CONTINUE

C   --Print the minimum and maximum

      IF (MMSTEP .LE. 0) THEN
         STR = 'ALL time steps'
      ELSE
         STR = 'this time step only'
      END IF
      LSTR = LENSTR(STR)
      IF (MMNUM .LE. 1) THEN
         IF (NOUT .GT. 0) THEN
            WRITE (NOUT, 10000, IOSTAT=IDUM)
     &         MMNAME(:LENSTR(MMNAME)), STR(:LSTR)
         ELSE
            WRITE (*, 10000, IOSTAT=IDUM)
     &         MMNAME(:LENSTR(MMNAME)), STR(:LSTR)
         END IF
      ELSE
         IF (NOUT .GT. 0) THEN
            WRITE (NOUT, 10000, IOSTAT=IDUM)
     &         MMNAME(:LENSTR(MMNAME)), STR(:LSTR), MMNUM
         ELSE
            WRITE (*, 10000, IOSTAT=IDUM)
     &         MMNAME(:LENSTR(MMNAME)), STR(:LSTR), MMNUM
         END IF
      END IF

      IF ((MMNUM .LE. 1) .AND. (XMIN .EQ. XMAX)) THEN
         IF (NOUT .GT. 0) THEN
            WRITE (NOUT, 10010, IOSTAT=IDUM) XMIN
         ELSE
            WRITE (*, 10010, IOSTAT=IDUM) XMIN
         END IF

      ELSE IF (MMTYP .EQ. 'G') THEN
         IF (MMSTEP .LE. 0) THEN
            IF (NOUT .GT. 0) THEN
               WRITE (NOUT, 10020, IOSTAT=IDUM)
     &            NMIN, 'Minimum', XMIN, MINSTE
               WRITE (NOUT, 10020, IOSTAT=IDUM)
     &            NMAX, 'Maximum', XMAX, MAXSTE
            ELSE
               WRITE (*, 10020, IOSTAT=IDUM)
     &            NMIN, 'Minimum', XMIN, MINSTE
               WRITE (*, 10020, IOSTAT=IDUM)
     &            NMAX, 'Maximum', XMAX, MAXSTE
            END IF
         ELSE
            IF (NOUT .GT. 0) THEN
               WRITE (NOUT, 10020, IOSTAT=IDUM)
     &            NMIN, 'Minimum', XMIN
               WRITE (NOUT, 10020, IOSTAT=IDUM)
     &            NMAX, 'Maximum', XMAX
            ELSE
               WRITE (*, 10020, IOSTAT=IDUM)
     &            NMIN, 'Minimum', XMIN
               WRITE (*, 10020, IOSTAT=IDUM)
     &            NMAX, 'Maximum', XMAX
            END IF
         END IF

      ELSE IF ((MMTYP .EQ. 'N') .OR. (MMTYP .EQ. 'E')) THEN
         IF (MMTYP .EQ. 'N') THEN
            STR = 'node'
         ELSE
            STR = 'element'
         END IF
         LSTR = LENSTR(STR)

         IF (MMSTEP .LE. 0) THEN
            IF (NOUT .GT. 0) THEN
               WRITE (NOUT, 10030, IOSTAT=IDUM)
     &            NMIN, 'Minimum', XMIN, STR(:LSTR), MINNE, MINSTE
               WRITE (NOUT, 10030, IOSTAT=IDUM)
     &            NMAX, 'Maximum', XMAX, STR(:LSTR), MAXNE, MAXSTE
            ELSE
               WRITE (*, 10030, IOSTAT=IDUM)
     &            NMIN, 'Minimum', XMIN, STR(:LSTR), MINNE, MINSTE
               WRITE (*, 10030, IOSTAT=IDUM)
     &            NMAX, 'Maximum', XMAX, STR(:LSTR), MAXNE, MAXSTE
            END IF
         ELSE
            IF (NOUT .GT. 0) THEN
               WRITE (NOUT, 10030, IOSTAT=IDUM)
     &            NMIN, 'Minimum', XMIN, STR(:LSTR), MINNE
               WRITE (NOUT, 10030, IOSTAT=IDUM)
     &            NMAX, 'Maximum', XMAX, STR(:LSTR), MAXNE
            ELSE
               WRITE (*, 10030, IOSTAT=IDUM)
     &            NMIN, 'Minimum', XMIN, STR(:LSTR), MINNE
               WRITE (*, 10030, IOSTAT=IDUM)
     &            NMAX, 'Maximum', XMAX, STR(:LSTR), MAXNE
            END IF
         END IF
      END IF

C   --Update the min/max number

      IF (XMIN .NE. XMAX) THEN
         MMNUM = MMNUM + 1
      ELSE
         MMNUM = -999
      END IF

  110 CONTINUE
      RETURN
10000  FORMAT (/, 1X, 'Variable ', A, ' for ', A,
     &   :, ', number =', I12)
10010  FORMAT (4X, 'All values = ', 1PE17.10)
10020  FORMAT (1X, I12, ' ', A, '(s) = ', 1PE17.10, :
     &   ', first at step ', I12)
10030  FORMAT (1X, I12, ' ', A, '(s) = ', 1PE17.10,
     &   ', first at ', A, I12, :, ' of step ', I12)
      END
