C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE DBISTE (NDB, OPTION, ISTEP, NELBLK, TIME,
     &                   NVARGL, NVARNP, NVAREL, NUMNP,
     &                   IDELB, NUMELB, ISEVOK,
     &                   VARGL, VARNP, VAREL, IOERR)
C=======================================================================
C   --*** DBISTE *** (EXOLIB) Read database variables for one time step
C   --   Written by Amy Gilkey - revised 10/14/87
C   --   Modified for ExodusIIv2 database format 10/16/95
C   --
C   --DBISTE reads the database global, nodal, and element variables
C   --for one time step.
C   --
C   --Parameters:
C   --   NDB    - IN  - the database number
C   --   OPTION - IN  - ' ' to not store
C   --                  '*' to store all, else store options:
C   --                  'G' to store global variables
C   --                  'E' to store element variables
C   --                  'N' to store nodal variables
C   --   ISTEP  - IN  - the time step number
C   --   NELBLK - IN  - the number of element blocks
C   --   TIME   - OUT - the time step time
C   --   NVARGL - IN  - the number of global variables
C   --   NVARNP - IN  - the number of nodal variables
C   --   NVAREL - IN  - the number of element variables
C   --   NUMNP  - IN  - the number of nodes
C   --   IDELB  - IN  - the element block ID's
C   --   NUMELB - IN  - the number of elements per block
C   --   ISEVOK - IN  - the element block variable truth table;
C   --                  variable i of block j exists iff ISEVOK(j,i)
C   --   VARGL  - OUT - the global variables for the time step (if OPTION)
C   --   VARNP  - OUT - the nodal variables for the time step (if OPTION)
C   --   VAREL  - OUT - the element variables for the time step (if OPTION)
C   --   IOERR  - OUT - I/O error flag

      INTEGER NDB
      CHARACTER*(*) OPTION
      INTEGER ISTEP
      INTEGER NELBLK
      REAL    TIME
      INTEGER NVARGL, NVARNP, NVAREL, NUMNP
      INTEGER IDELB(*)
      INTEGER NUMELB(*)
      LOGICAL ISEVOK(*)
      REAL    VARGL(*)
      REAL    VARNP(*)
      REAL    VAREL(*)
      INTEGER IOERR
      INTEGER EBID, NEPB

C   --Read step time

      CALL EXGTIM(NDB, ISTEP, TIME, IERR)

C     --Read global variables
      IF (NVARGL .GT. 0) THEN
         CALL EXGGV (NDB, ISTEP, NVARGL, VARGL, IERR)
      END IF

C     --Read nodal variables
      IF (NVARNP .GT. 0) THEN
         INP = 1
         DO 10 I = 1, NVARNP
            CALL EXGNV (NDB, ISTEP, I, NUMNP, VARNP(INP), IERR)
            INP = INP + NUMNP
 10      CONTINUE
      END IF

C      --Read element variables

      IF (NVAREL .GT. 0) THEN
         IEL = 1
         DO 20 I = 1, NELBLK
            EBID = IDELB(I)
            NEPB = NUMELB(I)
            DO 30 J = 1, NVAREL
               CALL EXGEV(NDB,ISTEP,J,EBID,NEPB,VAREL(IEL),IERR)
               IEL = IEL + NEPB
 30         CONTINUE
 20      CONTINUE
      END IF

      RETURN
      END
