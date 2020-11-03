C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE DBIST2 (NDB, NVAREL,  NELBLK, NELBDM, ISEVOK,
     $     VAREL, NUMELB, IVAR, IELB, *)
C=======================================================================
C   --*** DBIST2 *** (EXOLIB) Internal to DBISTE, Read element variables
C   --
C   --DBIST2 reads the database element variables for one time step.
C   --
C   --Parameters:
C   --   NDB - IN - the database number
C   --   NVAREL - IN - the number of element variables
C   --   NVARDM - IN - the row dimension of VAREL
C   --   NELBLK - IN - the number of element blocks
C   --   NELBDM - IN - the row dimension of ISEVOK
C   --   ISEVOK - IN - the element block variable truth table;
C   --      variable i of block j exists iff ISEVOK(j,i)
C   --   VAREL - OUT - the element variables for the time step (if OPTION)
C   --   IVAR  - OUT - the nodal variable number if error on read.
C   --   IELB  - OUT - the element block number if error on read.
C   --   * - return statement if error encountered, including end-of-file;
C   --      message is printed
C   --
      INTEGER NDB
      INTEGER NVAREL, NELBLK, NELBDM
      INTEGER NUMELB(*)
c      LOGICAL ISEVOK(nvarel,*)
      integer ISEVOK(nvarel,*)
      REAL VAREL(*)

      IELo = 0
      DO 130 IELB = 1, NELBLK
         DO 120 IVAR = 1, NVAREL
            IF (ISEVOK(IVAR,IELB) .ne. 0) THEN
               READ (NDB, END=200, ERR=200, IOSTAT=IERR)
     &              (VAREL(IELo+N), N=1,NUMELB(IELB))
               ielo=ielo+numelb(ielb)
            END IF
  120    CONTINUE
  130 CONTINUE
      RETURN
  200 CONTINUE
      RETURN 1
      END
