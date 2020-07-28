C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE SPREAD (A, NPTIMS, IPTIMS, NENUM, PLTVAL)
C=======================================================================

C   --*** SPREAD *** (SPLOT) Read plot variables from database
C   --   Written by Amy Gilkey - revised 02/02/88
C   --
C   --SPREAD reads the database and stores all variables needed
C   --for all time steps.
C   --
C   --This routine manipulates dynamic memory, so check after return.
C   --
C   --Parameters:
C   --   A - IN - the dynamic memory base array
C   --   NPTIMS - IN - the number of selected times
C   --   IPTIMS - IN - the selected time steps
C   --   NENUM - IN - the selected node/element numbers
C   --   PLTVAL - OUT - the plot variable values
C   --
C   --Common Variables:
C   --   Uses NDIM, NUMNP, NUMEL, NELBLK,
C   --      NVARHI, NVARGL, NVARNP, NVAREL of /DBNUMS/
C   --   Uses NODVAR, NNENUM of /SELNE/
C   --   Uses NSPVAR of /SPVARS/

      include 'params.blk'
      include 'dbnums.blk'
      include 'selne.blk'
      include 'spvars.blk'

      DIMENSION A(*)
      INTEGER IPTIMS(NPTIMS)
      INTEGER NENUM(NNENUM)
      REAL PLTVAL(NNENUM,NSPVAR,*)

      LOGICAL NEEDNV, NEEDEV

C   --Determine which types of variables are needed

      NEEDNV = .FALSE.
      NEEDEV = .FALSE.
      LDATA = 0

      IF (NODVAR) THEN
         NEEDNV = .TRUE.
         LDATA = MAX (LDATA, NUMNP)
      ELSE
         NEEDEV = .TRUE.
         LDATA = MAX (LDATA, NUMEL)
      END IF

C   --Reserve memory for data record

      CALL MDRSRV ('DATA', KDATA, LDATA)
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 110

C   --Transfer element variables onto random file (for efficiency)

      IF (NEEDEV) THEN
         CALL SPTRND (A, IPTIMS(NPTIMS), 'E', NUMEL, NVAREL, A(KDATA))
      END IF

      DO 100 NPT = 1, NPTIMS

         ISTEP = IPTIMS(NPT)

C      --Read and store nodal data to be plotted

         IF (NEEDNV) THEN
            CALL SPSTOR (A, ISTEP, 'N', NUMNP, NVARNP, NENUM,
     &         PLTVAL(1,1,NPT), A(KDATA))
         END IF

C      --Read and store element data to be plotted

         IF (NEEDEV) THEN
            CALL SPSTOR (A, ISTEP, 'E', NUMEL, NVAREL, NENUM,
     &         PLTVAL(1,1,NPT), A(KDATA))
         END IF

  100 CONTINUE

      CALL MDDEL ('DATA')

  110 CONTINUE
      RETURN
      END
