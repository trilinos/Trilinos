C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE WRSTEP (NTXT, ISTEP, NELBLK, TIME,
     &           NVARGL, NVARNP, NVAREL, NUMNP, IDELB, NUMELB, ISEVOK,
     &           VARGL, VARNP, VAREL, NGLDM, NNPDM, NELDM, IOERR)
C=======================================================================

C   --*** WRSTEP *** (TXTEXO) Write database variables for one time step
C   --   Written by Amy Gilkey - revised 03/02/88
C   --   Modified for ExodusIIv2 database format
C   --
C   --WRSTEP writes the database history, global, nodal, and element variables
C   --for one time step.
C   --
C   --Parameters:
C   --   NTXT   - IN - the text file
C   --   ISTEP  - IN - the time step number
C   --   NELBLK - IN - the number of element blocks
C   --   TIME   - IN - the time step time
C   --   NVARGL - IN - the number of global variables
C   --   NVARNP - IN - the number of nodal variables
C   --   NVAREL - IN - the number of element variables
C   --   NUMNP  - IN - the number of nodes
C   --   IDELB  - IN - element block ID's
C   --   NUMELB - IN - the number of elements per block
C   --   ISEVOK - IN - the element block variable truth table;
C   --                 variable i of block j exists iff ISEVOK(j,i)
C   --   VARGL  - IN - the global variables for the time step
C   --   VARNP  - IN - the nodal variables for the time step
C   --   VAREL  - IN - the element variables for the time step
C   --   NGLDM  - IN - dimension of global variable array
C   --   NNPDM  - IN - dimension of nodal variable array
C   --   NELDM  - IN - dimension of element variable array
C   --   IOERR  - IN/OUT - I/O error flag

      INTEGER NTXT, ISTEP, NELBLK
      REAL    TIME
      INTEGER NVARGL, NVARNP, NVAREL, NUMNP
      INTEGER IDELB(*)
      INTEGER NUMELB(*)
      LOGICAL ISEVOK(NELBLK,*)
      REAL    VARGL(NGLDM)
      REAL    VARNP(NUMNP,NNPDM)
      REAL    VAREL(*)
      INTEGER IOERR

C   --Write step time
      WRITE (NTXT, '(A, I10)') '! Time step', ISTEP
      WRITE (NTXT, '(1pe16.7, 16X,A)') TIME, '! Time step time value'

      IF ((NVARGL+NVARNP+NVAREL) .GT. 0) THEN

C      --Write global variables
         IF (NVARGL .GT. 0) THEN
            WRITE (NTXT, '(A)') '! Global variables'
            WRITE (NTXT, 10010) (VARGL(I), I=1,NVARGL)
         END IF

         IF (NVARNP .GT. 0) THEN
            WRITE (NTXT, '(A)') '! Nodal variables'
            DO 100 INP = 1, NUMNP
               WRITE (NTXT, 10010) (VARNP(INP,I), I=1,NVARNP)
  100       CONTINUE
         END IF

C      --Write element variables

         IF (NVAREL .GT. 0) THEN
            WRITE (NTXT, '(A)') '! Element variables'
            IEL0 = 0
            DO 120 IELB = 1, NELBLK
               NELPB = NUMELB(IELB)
               WRITE (NTXT, '(A, I3)') '! Element Block ', IELB
               WRITE (NTXT, '(A, I9)')
     &               '! Number of elements ', NELPB
               DO 110 N = 1, NELPB
C              This is a really ugly index.I cannot seem to simplify
C              this anymore.
C              The element variables are read in DBISTE with the
C              EXGEV call. IELO is the total number of elements in the
C              previous element block. N if the index for the number
C              of elements in the current element block.
                     WRITE (NTXT, 10010)
     &               (VAREL(IEL0+N+((I-1)*NELPB)), I=1,NVAREL)
  110          CONTINUE
               IEL0 = IEL0 + NELPB*NVAREL
  120       CONTINUE
         END IF

      END IF
10010 FORMAT (5(1pE16.7))

      RETURN
      END
