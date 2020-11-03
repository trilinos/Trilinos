C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details
C=======================================================================
      SUBROUTINE RDSTEP (ISTEP, TIME, NUMELB, IDELB, ISEVOK,
     &                   VISELB, MAXNE, VARVAL, MERR)
C=======================================================================

C   --*** RDSTEP *** (ALGEBRA) Read database time step variables
C   --   Written by Amy Gilkey - revised 11/30/87
C   --   Modified 8/30/95
C   --
C   --RDSTEP reads the input database variables for one time step and
C   --stores the ones used in the equations in array VARVAL.
C   --
C   --Parameters:
C   --   ISTEP  - IN  - the time step number
C   --   TIME   - IN  - the time step time
C   --   NUMELB - IN  - the number of elements per block
C   --   ISEVOK - IN  - the element block variable truth table;
C   --                  variable i of block j exists iff ISEVOK(j,i)
C   --   MAXNE  - IN  - the VARVAL dimension
C   --   VARVAL - OUT - the input data needed
C   --   MERR   - OUT - error flag
C   --
C   --Common Variables:
C   --   Uses ITIME, IGVBEG, INVBEG, IEVBEG, IGVEND, INVEND, IEVEND
C   --      of /DBXVAR/

      include 'ag_namlen.blk'
      include 'ag_var.blk'
      include 'ag_dbase.blk'
      include 'ag_dbnums.blk'
      include 'ag_dbxvar.blk'

      PARAMETER (ICURTM = 1, ILSTTM = 2, IONETM = 3)

      INTEGER NUMELB(*)
      INTEGER IDELB(*)
      LOGICAL ISEVOK(NELBLK*NVAREL)
      LOGICAL VISELB(NELBLK)
      REAL    VARVAL(MAXNE,*)
      INTEGER MERR
      MERR = 0

C   --Assign time to time/history/globals entry
      VARVAL(IDVAR(ITIME),ISTVAR(ICURTM,ITIME)) = TIME

C      --Read global variables

      IF (IGVBEG .LE. IGVEND) THEN
         CALL STORE (ISTEP, 'G', IGVBEG, IGVEND, NVARGL,
     &        NUMELB, IDELB, ISEVOK, VISELB, MAXNE, VARVAL, MERR)
         IF (MERR .EQ. 1) RETURN
      END IF

C      --Read nodal variables

      IF (INVBEG .LE. INVEND) THEN
         CALL STORE (ISTEP, 'N', INVBEG, INVEND, NUMNP,
     &        NUMELB, IDELB, ISEVOK, VISELB, MAXNE, VARVAL, MERR)
         IF (MERR .EQ. 1) RETURN
      END IF

C      --Read element variables

      IF (IEVBEG .LE. IEVEND) THEN
         CALL STORE (ISTEP, 'E', IEVBEG, IEVEND, NUMEL,
     &        NUMELB, IDELB, ISEVOK, VISELB, MAXNE, VARVAL, MERR)
         IF (MERR .EQ. 1) RETURN
      END IF

      RETURN
      END
