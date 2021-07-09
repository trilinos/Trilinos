C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE RWSTEP (NTXT, NDB, istep, idelb,
     &  NVARGL, NVARNP, NUMNP, NVAREL, NELBLK, NUMELB,
     &  ISEVOK, TIME, VAR, *)
C=======================================================================

C   --*** RDSTEP *** (TXTEXO) Read database variables for one time step
C   --   Written by Amy Gilkey - revised 03/02/88
C   --
C   --RDSTEP reads the database history, global, nodal, and element variables
C   --for one time step.
C   --
C   --Parameters:
C   --   NTXT - IN - the text file
C   --   NVARGL - IN - the number of global variables
C   --   NVARNP - IN - the number of nodal variables
C   --   NUMNP - IN - the number of nodes
C   --   NVAREL - IN - the number of element variables
C   --   NELBLK - IN - the number of element blocks
C   --   NUMELB - IN - the number of elements per block
C   --   ISEVOK - IN - the element block variable truth table;
C   --      variable i of block j exists iff ISEVOK(j,i)
C   --   TIME - OUT - the time step time
C   --   VAR - OUT - the global variables for the time step
C   --   * - return statement if error encountered, including end-of-file;
C   --      message is printed
C   --
C   --Database must be positioned in front of time step upon entry;
C   --upon exit positioned after time step.

      INTEGER NUMELB(*)
      integer idelb(*)
      LOGICAL ISEVOK(NVAREL, *)
      REAL VAR(*)

      CHARACTER*32 STRA

C   --Read step time

      READ (NTXT, *, END=180, ERR=130)
      READ (NTXT, *, END=180, ERR=130) TIME
      call exptim(ndb, istep, time, ierr)

C      --Read global variables

         IF (NVARGL .GT. 0) THEN
            READ (NTXT, *, END=150, ERR=150)
            READ (NTXT, *, END=150, ERR=150) (VAR(I), I=1,NVARGL)
            call expgv(ndb, istep, nvargl, var, ierr)
         END IF

C      --Read nodal variables

         IF (NVARNP .GT. 0) THEN
            READ (NTXT, *, END=160, ERR=160)
            DO 100 INP = 1, NUMNP
               READ (NTXT, *, END=160, ERR=160)
     &          (VAR(inp+(i-1)*numnp), I=1,NVARNP)
  100       CONTINUE

            ioff = 1
            do 105 i=1, nvarnp
              call expnv(ndb, istep, i, numnp, var(ioff), ierr)
              ioff = ioff + numnp
 105        continue
         END IF

C      --Read element variables
C ... Values for all element variables are stored in text file
C     Read an element block at a time and then write that element block
         IF (NVAREL .GT. 0) THEN
C ... Skip comment line
            READ (NTXT, *, END=170, ERR=170)
            DO 120 IELB = 1, NELBLK
C ... Skip comment line
              READ (NTXT, *, END=170, ERR=170)
              READ (NTXT, *, END=170, ERR=170)
               DO 110 N = 1, NUMELB(IELB)
                  READ (NTXT, *, END=170, ERR=170)
     &               (VAR(N+(I-1)*numelb(ielb)), I=1,NVAREL)
  110          CONTINUE
C ... Write this blocks variables
               ioff = 1
               do 115 i=1, nvarel
                 if (isevok(i,ielb)) then
                   call expev(ndb, istep, i, idelb(ielb), numelb(ielb),
     &               var(ioff), ierr)
                 end if
                 ioff = ioff + numelb(ielb)
 115           continue
  120       CONTINUE
         END IF

      RETURN

  130 CONTINUE
      CALL PRTERR ('FATAL', 'Reading TIME STEP TIME')
      GOTO 180
  150 CONTINUE
      CALL PRTERR ('FATAL', 'Reading GLOBAL VARIABLES')
      GOTO 180
  160 CONTINUE
      CALL INTSTR (1, 0, INP, STRA, LSTRA)
      CALL PRTERR ('FATAL',
     &   'Reading NODAL VARIABLES for node ' // STRA(:LSTRA))
      GOTO 180
  170 CONTINUE
      CALL INTSTR (1, 0, N, STRA, LSTRA)
      CALL PRTERR ('FATAL',
     &   'Reading ELEMENT VARIABLES for element ' // STRA(:LSTRA))
      GOTO 180
  180 CONTINUE
      RETURN 1
      END
