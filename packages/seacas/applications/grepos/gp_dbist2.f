C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE DBIST2 (NDB, TYPE, ISTEP, NVAREL, NVARDM, NELBLK,
     &  NELBDM, ISEVOK, VAREL, NUMELB, IDELB, IVAR, IELB, *)
C=======================================================================

C   --*** DBIST2 *** (EXOLIB) Internal to DBISTE, Read element variables
C   --   Written by Greg Sjaardema 8/8/90, to remove MAX from dimensions
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
      include 'exodusII.inc'
      INTEGER NDB
      CHARACTER*1 TYPE
      INTEGER NVAREL, NVARDM, NELBLK, NELBDM
      INTEGER NUMELB(*)
      INTEGER IDELB(*)
      LOGICAL ISEVOK(NELBDM,*)
      REAL VAREL(*)

      iel0 = 1
      do 120 ivar = 1, nvarel
        do 130 ielb = 1, nelblk

C ... To make it easier to handle consolidating blocks, we read in
C     or store data for all variables. If variable exists, read it,
C     if not, then zero-fill. Remember to do this on modification,
C     and on output also...
C     This lets us treat the VAREL array as a rectangular array
C     of dimension (NUMEL, NVAREL)

          id = idelb(ielb)
          if (isevok(ielb,ivar)) then
            if (type .eq. 'E') then
              call exgev (ndb, istep, ivar, idelb(ielb),
     &          numelb(ielb), varel(iel0), ierr)
            else if (type .eq. 'M') then
              call exgnsv (ndb, istep, ivar, idelb(ielb),
     &          numelb(ielb), varel(iel0), ierr)
            else if (type .eq. 'S') then
              call exgssv (ndb, istep, ivar, idelb(ielb),
     &          numelb(ielb), varel(iel0), ierr)
            end if
          else
            do 110 i = 0, numelb(ielb)-1
              varel(iel0 + i) = 0.0
 110        continue
          end if
          iel0 = iel0 + numelb(ielb)
 130    continue
 120  continue
      RETURN
      END
