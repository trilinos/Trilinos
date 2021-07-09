C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE LNSTOR (A, ISTEP, TYP, NWRDS, NVAR, NPT, NPTS,
     &   XLN, YLN, ZLN, DATA)
C=======================================================================

C   --*** LNSTOR *** (PATHLN) Read and store pathline data from database
C   --   Written by Amy Gilkey - revised 05/27/88
C   --
C   --LNSTOR reads variables from the database and stores any that are
C   --pathline data in the appropriate location.  It reads only one
C   --variable type (history, global, nodal, element) per call.
C   --
C   --Parameters:
C   --   A - IN - the dynamic memory base array
C   --   ISTEP - IN - the time step number
C   --   TYP - IN - the type of variable:
C   --      'H'istory, 'G'lobal, 'N'odal, 'E'lement
C   --   NWRDS - IN - the number of words in a data record
C   --   NVAR - IN - the number of records to read
C   --   NPT - IN - the XLN, YLN, ZLN time index to fill
C   --   NPTS - IN - the maximum XLN, YLN, ZLN time index
C   --   XLN, YLN, ZLN - IN/OUT - the pathline data array
C   --   DATA - SCRATCH - size = NWRDS
C   --
C   --Common Variables:
C   --   Uses NLNCRV, ILVNE, ILVID of /LNVARS/

      include 'lnvars.blk'
      include 'dbnums.blk'

      DIMENSION A(*)
      CHARACTER TYP
      REAL XLN(NPTS,NLNCRV), YLN(NPTS,NLNCRV), ZLN(NPTS,NLNCRV)
      REAL DATA(NWRDS)

      LOGICAL NEED
      CHARACTER T

      CALL DBVIX_BL (TYP, 1, ISID)
      CALL DBVIX_BL (TYP, NVAR, IEID)

      IF ((TYP .EQ. 'H') .OR. (TYP .EQ. 'G')) THEN

C      --Read history or global variables

         CALL GETVAR (A, ISID, -999, ISTEP, NWRDS, DATA)

C      --Scan plot variable information and store if a history or global
C      --point from this data is to be plotted

         DO 100 NP = 1, NLNCRV
            CALL DBVTYP_BL (ILVID(1,NP), T, IDUM)
            IF ((T .EQ. 'H') .OR. (T .EQ. 'G')) THEN
               CALL DBVTYP_BL (ILVID(1,NP), T, NE)
               XLN(NPT,NP) = DATA(NE)
               CALL DBVTYP_BL (ILVID(2,NP), T, NE)
               YLN(NPT,NP) = DATA(NE)
               IF (NDIM .GE. 3) THEN
                  CALL DBVTYP_BL (ILVID(3,NP), T, NE)
                  ZLN(NPT,NP) = DATA(NE)
               END IF
            END IF
  100    CONTINUE

      ELSE
         DO 130 ID = ISID, IEID

C         --Determine if variable is needed

            NEED = .FALSE.
            DO 110 NP = 1, NLNCRV
               NEED = NEED .OR. (LOCINT (ID, NDIM, ILVID(1,NP)) .GT. 0)
  110       CONTINUE

            IF (NEED) THEN

C            --Read nodal/element variable

               CALL GETVAR (A, ID, -1, ISTEP, NWRDS, DATA)

C            --Scan plot variable information and store if a node/element
C            --point from this data is to be plotted

               DO 120 NP = 1, NLNCRV
                  NE = ILVNE(NP)
                  IF (ID .EQ. ILVID(1,NP)) XLN(NPT,NP) = DATA(NE)
                  IF (ID .EQ. ILVID(2,NP)) YLN(NPT,NP) = DATA(NE)
                  IF (NDIM .GE. 3) THEN
                     IF (ID .EQ. ILVID(3,NP)) ZLN(NPT,NP) = DATA(NE)
                  END IF
  120          CONTINUE
            END IF

  130    CONTINUE
      END IF

      RETURN
      END
