C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE DBIST1 (NDB, NVARNP, NUMNP, VARNP, IVAR, *)
C=======================================================================
C   --*** DBIST1 *** (EXOLIB) Internal to DBISTE, Read nodal variables
C   --
C   --DBIST1 reads the database nodal variables for one time step.
C   --
C   --Parameters:
C   --   NDB - IN - the database number
C   --   NVARNP - IN - the number of nodal variables
C   --   NUMNP - IN - the number of nodes
C   --   VARNP - OUT - the nodal variables for the time step (if OPTION)
C   --   IVAR  - OUT - the nodal variable number if error on read.
C   --   * - return statement if error encountered, including end-of-file;
C   --      message is printed
C   --
      INTEGER NDB
      INTEGER NVARNP, NUMNP
      REAL VARNP(*)

      iov = 0
      DO 100 IVAR = 1, NVARNP
         READ (NDB, END=190, ERR=190, IOSTAT=IERR)
     &        (VARNP(iov+INP), INP=1,NUMNP)
         iov = iov+numnp
  100 CONTINUE
      RETURN
  190 CONTINUE
      RETURN 1
      END
