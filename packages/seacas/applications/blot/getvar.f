C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE GETVAR (A, IVAR, IELBLK, INSTEP, LENVAR, VAR)
C=======================================================================

C   --*** GETVAR *** (BLOT) Read variable
C   --   Written by Amy Gilkey - revised 05/17/88
C   --
C   --GETVAR returns the values for the requested variable for the
C   --requested time step.  It either reads the values from the sequential
C   --database file and writes them to a direct access scratch file or it
C   --reads the values from the direct access scratch file.
C   --
C   --This routine uses MDFIND to find the following dynamic memory arrays:
C   --   NUMELB - the number of elements per element block
C   --   ISEVOK - the element block variable truth table;
C   --      variable i of block j exists iff ISEVOK(j,i)
C   --   WHOTIM - true iff whole (versus history) time step
C   --
C   --Parameters:
C   --   A - IN - the dynamic memory base array
C   --   IVAR - IN - the variable index
C   --      (0 to initialize random file only)
C   --   IELBLK - IN - the element block number, <=0 for all
C   --      (for element blocks only)
C   --   INSTEP - IN - the time step number
C   --      = +n to read time step n
C   --      = -n to transfer time step n to random file only
C   --      =  0 to transfer all time steps to random file
C   --   LENVAR - IN - the length of VAR
C   --   VAR - OUT - the variable values (indeterminate if INSTEP <= 0)
C   --
C   --Common Variables:
C   --   Uses NDB of /DBASE/
C   --   Uses NUMNP, NUMEL, NELBLK, NVARHI, NVARGL, NVARNP, NVAREL,
C   --      NSTEPS, NSTEPW of /DBNUMS/
C   --
C   --Database is rewound upon the first entry of this routine; upon
C   --exit a flag is set to keep track of the database position; the
C   --database should not be moved between calls to this routine.
C   --
C   --A scratch random file is created and read and written in this routine.
C   --It is connected to unit 90.

      DIMENSION A(*)
      REAL VAR(*)

      IF (IVAR .LE. 0) RETURN
      CALL RNDVAR (A, A, A, IVAR, IELBLK, INSTEP, LENVAR, VAR)

      RETURN
      END
