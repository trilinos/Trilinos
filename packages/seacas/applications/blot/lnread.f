C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE LNREAD (A, NPTIMS, NPTIMW, IPTIMS, TIMES, WHOTIM,
     &   XLN, YLN, ZLN)
C=======================================================================

C   --*** LNREAD *** (PATHLN) Read pathline data from database
C   --   Written by Amy Gilkey - revised 05/27/88
C   --
C   --LNREAD reads the database and stores the pathline data.
C   --
C   --This routine manipulates dynamic memory, so check after return.
C   --
C   --Parameters:
C   --   A - IN - the dynamic memory base array
C   --   NPTIMS - IN - the number of selected steps (for history variables)
C   --   NPTIMW - IN - the number of whole selected steps
C   --   IPTIMS - IN - the selected time steps
C   --   TIMES - IN - the time step times
C   --   WHOTIM - IN - true iff whole (versus history) time step
C   --   XLN, YLN, ZLN - OUT - the pathline data
C   --
C   --Common Variables:
C   --   Uses NDIM, NUMNP, NUMEL, NELBLK,
C   --      NVARHI, NVARGL, NVARNP, NVAREL of /DBNUMS/
C   --   Uses NLNCRV, ILVID of /LNVARS/

      include 'dbnums.blk'
      include 'lnvars.blk'

      DIMENSION A(*)
      INTEGER IPTIMS(NPTIMS)
      REAL TIMES(*)
      LOGICAL WHOTIM(*)
      REAL XLN(NPTIMS,NLNCRV), YLN(NPTIMS,NLNCRV), ZLN(NPTIMS,NLNCRV)

      LOGICAL NEEDHV, NEEDGV, NEEDNV, NEEDEV
      CHARACTER TYP

C   --Determine which types of variables are needed

      NEEDHV = .FALSE.
      NEEDGV = .FALSE.
      NEEDNV = .FALSE.
      NEEDEV = .FALSE.
      LDATA = 0

      DO 100 NP = 1, NLNCRV
         CALL DBVTYP_BL (ILVID(1,NP), TYP, IDUM)
         IF (TYP .EQ. 'H') THEN
            NEEDHV = .TRUE.
            LDATA = MAX (LDATA, NVARHI)
         ELSE IF (TYP .EQ. 'G') THEN
            NEEDGV = .TRUE.
            LDATA = MAX (LDATA, NVARGL)
         ELSE IF (TYP .EQ. 'N') THEN
            NEEDNV = .TRUE.
            LDATA = MAX (LDATA, NUMNP)
         ELSE IF (TYP .EQ. 'E') THEN
            NEEDEV = .TRUE.
            LDATA = MAX (LDATA, NUMEL)
         END IF
  100 CONTINUE

C   --Reserve memory for data record

      CALL MDRSRV ('DATA', KDATA, LDATA)
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 130

C   --Transfer element variables onto random file (for efficiency)

      MXSTEP = 0
      DO 110 NPT = 1, NPTIMS
         ISTEP = IPTIMS(NPT)
         IF (WHOTIM(ISTEP)) MXSTEP = ISTEP
  110 CONTINUE

      IF (NEEDEV .AND. (MXSTEP .GT. 0)) THEN
C????         CALL LNTRND (A, MXSTEP, 'E', NUMEL, NVAREL, A(KDATA))
      END IF

      NPTIMW = 0
      DO 120 NPT = 1, NPTIMS

         ISTEP = IPTIMS(NPT)
         IF (WHOTIM(NPT)) NPTIMW = NPTIMW + 1

C      --Read and store variable data to be plotted

         IF (NEEDHV) THEN
            CALL LNSTOR (A, ISTEP, 'H', NVARHI, 1,
     &         NPT, NPTIMS, XLN, YLN, ZLN, A(KDATA))
         END IF

         IF (WHOTIM(ISTEP)) THEN
            IF (NEEDGV) THEN
               CALL LNSTOR (A, ISTEP, 'G', NVARGL, 1,
     &            NPTIMW, NPTIMS, XLN, YLN, ZLN, A(KDATA))
            END IF

            IF (NEEDNV) THEN
               CALL LNSTOR (A, ISTEP, 'N', NUMNP, NVARNP,
     &            NPTIMW, NPTIMS, XLN, YLN, ZLN, A(KDATA))
            END IF

            IF (NEEDEV) THEN
               CALL LNSTOR (A, ISTEP, 'E', NUMEL, NVAREL,
     &            NPTIMW, NPTIMS, XLN, YLN, ZLN, A(KDATA))
            END IF
         END IF

  120 CONTINUE

      CALL MDDEL ('DATA')

  130 CONTINUE
      RETURN
      END
