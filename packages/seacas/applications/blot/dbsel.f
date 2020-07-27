C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE DBSEL (A, IA, INLINE,
     &   SELTYP, IFLD, INTYP, CFIELD, IFIELD, RFIELD,
     &   NAMES, TIMES, WHOTIM, NPTIMS, IPTIMS,
     &   IDELB, LENE, IDNPS, IDESS,
     &   NCSTEP, LISNP, NLISEL, LISEL, LISNPS, LISESS,
     &   LISHV, LISGV, LISNV, LISEV, MAPEL, MAPND)
C=======================================================================

C   --*** DBSEL *** (BLOT) Process SELECT commands
C   --   Written by Amy Gilkey - revised 03/07/88
C   --
C   --DBSEL inputs and executes a SELECT command.
C   --
C   --Parameters:
C   --   A - IN - the dynamic memory base array
C   --   INLINE - IN/OUT - the parsed input lines for the log file
C   --   SELTYP - IN - the selection type
C   --   IFLD, INTYP, CFIELD, IFIELD, RFIELD - IN/OUT - the free-field
C   --      reader index and fields
C   --   NAMES - IN - the variable names
C   --   TIMES - IN - the times for all time steps
C   --   WHOTIM - IN - true iff whole (versus history) time step
C   --   NPTIMS - IN - the number of selected time steps
C   --   IPTIMS - IN - the selected time steps
C   --   IDELB - IN - the element block ID for each block
C   --   LENE - IN - the cumulative element counts by element block
C   --   IDNPS - IN - the node set ID for each set
C   --   IDESS - IN - the side set ID for each set
C   --   NCSTEP - IN/OUT - the current step number for display
C   --   LISNP - IN/OUT - the indices of the selected coordinates
C   --   NLISEL - IN/OUT - the number of selected elements for each block
C   --   LISEL - IN/OUT - the indices of the selected elements (by block)
C   --   LISNPS - IN/OUT - the indices of the selected node sets
C   --   LISESS - IN/OUT - the indices of the selected side sets
C   --   LISHV - IN/OUT - the indices of the selected history variables
C   --   LISGV - IN/OUT - the indices of the selected global variables
C   --   LISNV - IN/OUT - the indices of the selected nodal variables
C   --   LISEV - IN/OUT - the indices of the selected element variables
C   --
C   --Common Variables:
C   --   Uses NDIM, NUMNP, NUMEL, NELBLK, NVARHI, NVARGL, NVARNP, NVAREL,
C   --      NSTEPS of /DBNUMS/
C   --   Uses NUMNPS, NUMESS of /DBNUMG/

      include 'params.blk'
      include 'dbnums.blk'
      include 'dbnumgq.blk'

      DIMENSION A(*)
      INTEGER IA(*)

      CHARACTER*(*) INLINE(*)
      CHARACTER*(*) SELTYP
      INTEGER     INTYP(*)
      CHARACTER*(*) CFIELD(*)
      INTEGER     IFIELD(*)
      REAL        RFIELD(*)
      CHARACTER*(*) NAMES(*)
      REAL TIMES(*)
      LOGICAL WHOTIM(*)
      INTEGER IPTIMS(*)
      INTEGER IDELB(*)
      INTEGER LENE(0:*)
      INTEGER IDNPS(*)
      INTEGER IDESS(*)
      INTEGER LISNP(0:*)
      INTEGER NLISEL(0:*), LISEL(0:*)
      INTEGER LISNPS(0:*), LISESS(0:*)
      INTEGER LISHV(0:*), LISGV(0:*), LISNV(0:*), LISEV(0:*)
      INTEGER MAPEL(*)
      INTEGER MAPND(*)

      CHARACTER*(MXSTLN) STYP

      LOGICAL FIRST
      SAVE FIRST
C      --FIRST - true iff first time through routine

      CHARACTER*(MXSTLN) CMDTBL(14)
      SAVE CMDTBL
C      --CMDTBL - the valid commands table

      DATA FIRST / .TRUE. /

C   --Command table follows.  Remember to change the dimensioned size when
C   --changing the table.
      DATA CMDTBL /
     1   'NODES   ', 'ELEMENTS', 'BLOCKS  ', 'MATERIAL',
     2   'NSETS   ', 'SSETS   ',
     3   'HVARS   ', 'GVARS   ', 'NVARS   ', 'EVARS   ',
     4   'READ    ', 'STEP    ', 'TIME    ',
     5   '        ' /

C   --Get the command STYP

      CALL ABRSTR (STYP, SELTYP, CMDTBL)
      IF (STYP .EQ. ' ') STYP = SELTYP

C *** Initialization ***

C   --Reset parameters

      IF (FIRST .OR. (STYP .EQ. 'RESET')) THEN

         LISNP(0) = NUMNP
         DO 100 I = 1, NUMNP
            LISNP(I) = I
  100    CONTINUE
         NLISEL(0) = NELBLK
         DO 110 I = 1, NELBLK
            NLISEL(I) = LENE(I) - LENE(I-1)
  110    CONTINUE
         LISEL(0) = NUMEL
         DO 120 I = 1, NUMEL
            LISEL(I) = I
  120    CONTINUE
         LISNPS(0) = NUMNPS
         DO 130 I = 1, NUMNPS
            LISNPS(I) = I
  130    CONTINUE
         LISESS(0) = NUMESS
         DO 140 I = 1, NUMESS
            LISESS(I) = I
  140    CONTINUE
         IF (EXODUS) THEN
            LISHV(0) = NVARHI
            DO 150 I = 1, NVARHI
               LISHV(I) = I
  150       CONTINUE
            LISGV(0) = NVARGL
            DO 160 I = 1, NVARGL
               LISGV(I) = I
  160       CONTINUE
            LISNV(0) = NVARNP
            DO 170 I = 1, NVARNP
               LISNV(I) = I
  170       CONTINUE
            LISEV(0) = NVAREL
            DO 180 I = 1, NVAREL
               LISEV(I) = I
  180       CONTINUE
         END IF

         NCSTEP = MIN (1, NSTEPS)

         FIRST = .FALSE.
      END IF

      IF ((STYP .EQ. 'RESET') .OR. (STYP .EQ. 'reset')) GOTO 230

C *** GENESIS Print Commands ***

      IF (STYP .EQ. 'NODES') THEN
         CALL FFADDC (STYP, INLINE(1))
         CALL CKNONE (NUMNP, .FALSE., 'nodes', *220)

         CALL RMIXINT (INLINE(1), IFLD, INTYP, CFIELD, IFIELD,
     &      'node number', NUMNP, LISNP(0), LISNP(1), MAPND, *220)

         IF (LISNP(0) .GT. 0) THEN
            CALL PRNSEL (LISNP(0), NUMNP, 'nodes')
         END IF

      ELSE IF (STYP .EQ. 'ELEMENTS') THEN
         CALL FFADDC (STYP, INLINE(1))
         CALL CKNONE (NUMEL, .FALSE., 'elements', *220)

         CALL MDRSRV ('SCRSEL', KLEL, 1+NUMEL)
         CALL MDSTAT (NERR, MEM)
         IF (NERR .GT. 0) GOTO 220

         CALL RMIXINT (INLINE(1), IFLD, INTYP, CFIELD, IFIELD,
     &      'element number', NUMEL, IA(KLEL), IA(KLEL+1), MAPEL, *190)
  190    CONTINUE

         CALL DBSBEL (NELBLK, NUMEL, LENE, A(KLEL), NLISEL, LISEL)

         CALL MDDEL ('SCRSEL')

         IF (NLISEL(0) .GT. 0) THEN
            CALL PRNSEL (NLISEL(0), NELBLK, 'element blocks')
            CALL PRNSEL (LISEL(0), NUMEL, 'elements')
         END IF

      ELSE IF ((STYP .EQ. 'BLOCKS') .OR. (STYP .EQ. 'MATERIAL')) THEN
         CALL FFADDC (STYP, INLINE(1))
         CALL CKNONE (NELBLK, .FALSE., 'element blocks', *220)

         CALL MDRSRV ('SCRSEL', KLELB, 1+NELBLK)
         CALL MDSTAT (NERR, MEM)
         IF (NERR .GT. 0) GOTO 220

         CALL RIXID (INLINE(1), IFLD, INTYP, CFIELD, IFIELD,
     &      'element block ID',
     &      NELBLK, IDELB, IA(KLELB), IA(KLELB+1), *200)
  200    CONTINUE

         CALL DBSELB (NELBLK, NUMEL, LENE, A(KLELB), NLISEL, LISEL)

         CALL MDDEL ('SCRSEL')

         IF (NLISEL(0) .GT. 0) THEN
            CALL PRNSEL (NLISEL(0), NELBLK, 'element blocks')
            CALL PRNSEL (LISEL(0), NUMEL, 'elements')
         END IF

      ELSE IF (STYP .EQ. 'NSETS') THEN
         CALL FFADDC (STYP, INLINE(1))
         CALL CKNONE (NUMNPS, .FALSE., 'node sets', *220)

         CALL RIXID (INLINE(1), IFLD, INTYP, CFIELD, IFIELD,
     &      'node set ID',
     &      NUMNPS, IDNPS, LISNPS(0), LISNPS(1), *220)

         IF (LISNPS(0) .GT. 0) THEN
            CALL PRNSEL (LISNPS(0), NUMNPS, 'node sets')
         END IF

      ELSE IF (STYP .EQ. 'SSETS') THEN
         CALL FFADDC (STYP, INLINE(1))
         CALL CKNONE (NUMESS, .FALSE., 'side sets', *220)

         CALL RIXID (INLINE(1), IFLD, INTYP, CFIELD, IFIELD,
     &      'side set ID',
     &      NUMESS, IDESS, LISESS(0), LISESS(1), *220)

         IF (LISESS(0) .GT. 0) THEN
            CALL PRNSEL (LISESS(0), NUMESS, 'side sets')
         END IF

C *** EXODUS Print Commands ***

      ELSE IF (STYP .EQ. 'HVARS') THEN
         CALL FFADDC (STYP, INLINE(1))
         CALL CKEXOD (EXODUS, *220)
         CALL CKNONE (NVARHI, .FALSE., 'history variables', *220)

         CALL DBVIX_BL ('H', 1, IXHV)

         CALL RIXWRD (INLINE(1), IFLD, INTYP, CFIELD,
     &      'history variable name', NVARHI, NAMES(IXHV),
     &      LISHV(0), LISHV(1), *220)

         IF (LISHV(0) .GT. 0) THEN
            CALL PRNSEL (LISHV(0), NVARHI, 'history variables')
         END IF

      ELSE IF (STYP .EQ. 'GVARS') THEN
         CALL FFADDC (STYP, INLINE(1))
         CALL CKEXOD (EXODUS, *220)
         CALL CKNONE (NVARGL, .FALSE., 'global variables', *220)

         CALL DBVIX_BL ('G', 1, IXGV)

         CALL RIXWRD (INLINE(1), IFLD, INTYP, CFIELD,
     &      'global variable name', NVARGL, NAMES(IXGV),
     &      LISGV(0), LISGV(1), *220)

         IF (LISGV(0) .GT. 0) THEN
            CALL PRNSEL (LISGV(0), NVARGL, 'global variables')
         END IF

      ELSE IF (STYP .EQ. 'NVARS') THEN
         CALL FFADDC (STYP, INLINE(1))
         CALL CKEXOD (EXODUS, *220)
         CALL CKNONE (NVARNP, .FALSE., 'nodal variables', *220)

         CALL DBVIX_BL ('N', 1, IXNV)

         CALL RIXWRD (INLINE(1), IFLD, INTYP, CFIELD,
     &      'nodal variable name', NVARNP, NAMES(IXNV),
     &      LISNV(0), LISNV(1), *220)

         IF (LISNV(0) .GT. 0) THEN
            CALL PRNSEL (LISNV(0), NVARNP, 'nodal variables')
         END IF

      ELSE IF (STYP .EQ. 'EVARS') THEN
         CALL FFADDC (STYP, INLINE(1))
         CALL CKEXOD (EXODUS, *220)
         CALL CKNONE (NVAREL, .FALSE., 'element variables', *220)

         CALL DBVIX_BL ('E', 1, IXEV)

         CALL RIXWRD (INLINE(1), IFLD, INTYP, CFIELD,
     &      'element variable name', NVAREL, NAMES(IXEV),
     &      LISEV(0), LISEV(1), *220)

         IF (LISEV(0) .GT. 0) THEN
            CALL PRNSEL (LISEV(0), NVAREL, 'element variables')
         END IF

C *** EXODUS Movement Commands ***

      ELSE IF ((STYP .EQ. 'READ') .OR. (STYP .EQ. 'STEP')
     &   .OR. (STYP .EQ. 'TIME')) THEN
         nstep = 0
         CALL FFADDC (STYP, INLINE(1))
         CALL CKEXOD (EXODUS, *220)
         CALL CKNONE (NSTEPS, .FALSE., 'time steps', *220)

         IF (STYP .EQ. 'READ') THEN
            CALL FFINTG (IFLD, INTYP, IFIELD,
     &         'number of steps to read', 1, N, *210)
            CALL FFADDI (N, INLINE(1))
            NSTEP = MAX (1, NCSTEP + N)
         ELSE IF (STYP .EQ. 'STEP') THEN
            CALL FFINTG (IFLD, INTYP, IFIELD,
     &         'step number', NCSTEP, N, *210)
            CALL FFADDI (N, INLINE(1))
            NSTEP = N
         ELSE IF (STYP .EQ. 'TIME') THEN
            CALL FFREAL (IFLD, INTYP, RFIELD,
     &         'time step time', TIMES(NCSTEP), T, *210)
            CALL FFADDR (T, INLINE(1))
            NSTEP = LOCREA (T, NSTEPS, TIMES)
         END IF

         NCSTEP = NSTEP
         IF (NCSTEP .GT. NSTEPS) THEN
            CALL PRTERR ('CMDERR',
     &         'All time steps have been read from the database')
            NCSTEP = NSTEPS
         END IF

  210    CONTINUE
         CALL PRSTEP ('*', -1,
     &      TIMES(NCSTEP), WHOTIM(NCSTEP), NCSTEP, NSTEPS)

      ELSE IF (STYP .NE. 'RESET') THEN
         CALL SHOCMD ('SELECT Options:', CMDTBL)
         GOTO 220
      END IF

      GOTO 230

  220 CONTINUE
      INLINE(1) = ' '

  230 CONTINUE
      RETURN
      END
