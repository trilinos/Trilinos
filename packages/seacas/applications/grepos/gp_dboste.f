C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE DBOSTE (NDB, ISTEP,
     &  NVARGL, NVARNP, NUMNP, NVAREL, NELEXT,
     *  NELBLK, NUMELB, ISEVOK, IDELB,
     *  NVARNS, NUMNPS, NNNPS,  ISNSOK, IDNPS,
     *  NVARSS, NUMESS, NEESS,  ISSSOK, IDESS,
     *  TIME, VARGL, VARNP, VAREL, VARNS, VARSS, VARXT)
C=======================================================================

C   --*** DBOSTE *** (EXOLIB) Write database variables for one time step
C   --   Written by Amy Gilkey - revised 10/14/87
C   --   Modified by Greg Sjaardema - 10/2/90
C   --      Removed MAX from dimension statements, added routine
C   --      DBOST1 to do all work.
C   --
C   --DBOSTE writes the database history, global, nodal, and element variables
C   --for one time step.
C   --
C   --Parameters:
C   --   NDB - IN - the database number
C   --   ISTEP - IN - the time step number
C   --   NVARGL - IN - the number of global variables
C   --   NVARNP - IN - the number of nodal variables
C   --   NUMNP - IN - the number of nodes
C   --   NVAREL - IN - the number of element variables
C   --   NELBLK - IN - the number of element blocks
C   --   NUMELB - IN - the number of elements per block
C   --   ISEVOK - IN - the element block variable truth table;
C   --      variable i of block j exists iff ISEVOK(j,i)
C   --   TIME - IN - the time step time
C   --   VARGL - IN - the global variables for the time step
C   --   VARNP - IN - the nodal variables for the time step
C   --   VAREL - IN - the element variables for the time step
C   --
C   --Database must be positioned in front of time step upon entry;
C   --upon exit positioned after time step.

      include 'exodusII.inc'
      INTEGER NDB
      INTEGER ISTEP
      INTEGER NVARGL
      INTEGER NVARNP, NUMNP
      INTEGER NVAREL, NELBLK, NUMELB(*), IDELB(*)
      LOGICAL ISEVOK(NELBLK,*)
      INTEGER NVARNS, NUMNPS, NNNPS(*),  IDNPS(*)
      LOGICAL ISNSOK(NUMNPS,*)
      INTEGER NVARSS, NUMESS, NEESS(*),  IDESS(*)
      LOGICAL ISSSOK(NUMESS,*)
      REAL TIME
      REAL VARGL(*)
      REAL VARNP(*)
      REAL VAREL(*)
      REAL VARNS(*)
      REAL VARSS(*)
      REAL VARXT(*)

C     --Write step time
      CALL EXPTIM (NDB, ISTEP, TIME, IERR)

C     --Write global variables

      IF (NVARGL .GT. 0) THEN
         call expgv (ndb, istep, nvargl, vargl, ierr)
      END IF

C     --Write nodal variables

      INP0 = 1
      DO 100 IVAR = 1, NVARNP
        CALL EXPNV (NDB, ISTEP, IVAR, NUMNP, VARNP(INP0), IERR)
        INP0 = INP0 + NUMNP
 100  CONTINUE

C     --Write element variables
      iel0 = 1
      do 110 i = 1, nvarel
        do 120 ielb = 1, nelblk
          if (isevok(ielb,i)) then
             if (numelb(ielb) .gt. 0) then
                call expev (ndb, istep, i, idelb(ielb), numelb(ielb),
     &               varel(iel0), ierr)
             end if
          end if
C ... See dbist2.f for details. There is a value stored for each
C     element variable for each block even if the truth table
C     if false. We only output when the truth table is true.
C     We need to skip over the values that shouldn't be written.
          iel0 = iel0 + numelb(ielb)
 120    continue
 110  continue

      iel0 = 1
      do i = 1, nelext
        do ielb = 1, nelblk
           if (numelb(ielb) .gt. 0) then
              call expev (ndb, istep, nvarel+i, idelb(ielb),
     $             numelb(ielb), varxt(iel0), ierr)
           end if
           iel0 = iel0 + numelb(ielb)
        end do
      end do

C     --Write nodeset variables
      iel0 = 1
      do i = 1, nvarns
        do ielb = 1, numnps
          if (isnsok(ielb,i)) then
            if (nnnps(ielb) .gt. 0) then
              call expnsv (ndb, istep, i, idnps(ielb), nnnps(ielb),
     &          varns(iel0), ierr)
            end if
          end if
          iel0 = iel0 + nnnps(ielb)
        end do
      end do

C     --Write sideset variables
      iel0 = 1
      do i = 1, nvarss
        do ielb = 1, numess
          if (isssok(ielb,i)) then
            if (neess(ielb) .gt. 0) then
              call expssv (ndb, istep, i, idess(ielb), neess(ielb),
     &          varss(iel0), ierr)
            end if
          end if
          iel0 = iel0 + neess(ielb)
        end do
      end do

      RETURN
      END
