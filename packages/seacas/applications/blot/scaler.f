C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C======================================================================
      SUBROUTINE SCALER (A, IA, IPRINT, NAME, IVAR,
     &  USESEL, IELBST, NALVAR, VALMIN, VALMAX, MAPEL, MAPND)
C=======================================================================

C   --*** SCALER *** (BLOT) Scale variable
C   --   Written by Amy Gilkey - revised 04/01/88
C   --   D. P. Flanagan, 06/23/83
C   --
C   --SCALER finds the minimum and maximum values over the entire database
C   --for the desired variable.  The results are printed to the screen.
C   --The minimum and maximum values are stored and are returned but not
C   --printed if the routine is called twice for the same variable.
C   --
C   --For elements, the minimum and maximum are calculated for each
C   --selected element block.  All calculated minimums and maximums or stored
C   --and the cumulative values for all selected element blocks are returned.
C   --Element birth/death is also considered.
C   --
C   --This routine uses MDFIND to find the following dynamic memory arrays:
C   --   TIMES - the database time step times
C   --   WHOTIM - true iff the time step is a whole (versus history) time step
C   --   XN, YN, ZN - the nodal coordinates
C   --   XE, YE, ZE - the element coordinates
C   --   LENE - the cumulative element counts by element block
C   --   ISEVOK - the element block variable truth table;
C   --      variable i of block j exists iff ISEVOK(j,i)
C   --
C   --Parameters:
C   --   A - IN - the dynamic memory base array
C   --   IPRINT - IN - the min/max display flag:
C   --      0 = calculate only, do not display
C   --      1 = print only if not previously calculated
C   --      2 = always print
C   --   NAME - IN - the variable name
C   --   IVAR - IN - the variable index
C   --   USESEL - IN - use the element blocks selected array iff true,
C   --      else all selected
C   --   IELBST - IN - the element block status (>0 if selected)
C   --      (element variable only)
C   --   NALVAR - IN - the element birth/death variable (element variable only)
C   --   VALMIN, VALMAX - OUT - the variable minimum and maximum
C   --
C   --Common Variables:
C   --   Uses NDIM, NUMNP, NUMEL, NELBLK, NVARNP, NVAREL, NSTEPS of /DBNUMS/
C   --   Uses NAMECO of /DBNAMS/

      include 'dbnums.blk'
      include 'dbnams.blk'

      DIMENSION A(*)
      INTEGER IA(*)
      CHARACTER*(*) NAME
      LOGICAL USESEL
      INTEGER IELBST(*)
      INTEGER MAPEL(*), MAPND(*)

      CHARACTER TYP
      REAL XYZMIN(3), XYZMAX(3)
      REAL RDUM(3)

      SAVE KTIMES, KWHOLE, KXN, KYN, KZN, KXE, KYE, KZE, KLENE, KIEVOK
      SAVE KVALHN, KISTHN, KVALHX, KISTHX
      SAVE KVALGN, KISTGN, KVALGX, KISTGX
      SAVE KVALNN, KNUMNN, KXYZNN, KISTNN,
     &  KVALNX, KNUMNX, KXYZNX, KISTNX
      SAVE KVALEN, KNUMEN, KXYZEN, KISTEN,
     &  KVALEX, KNUMEX, KXYZEX, KISTEX
      SAVE IXALIV

      LOGICAL HIST1, GLOB1, NODE1, ELEM1
      SAVE HIST1, GLOB1, NODE1, ELEM1

      INTEGER NALOLD
      SAVE NALOLD

      DATA HIST1, GLOB1, NODE1, ELEM1 / .TRUE., .TRUE., .TRUE., .TRUE. /
      DATA NALOLD / -999 /

      IF (HIST1 .AND. GLOB1 .AND. NODE1 .AND. ELEM1) THEN
C      --Get times storage
        CALL MDFIND ('TIMES', KTIMES, IDUM)
        CALL MDFIND ('WHOTIM', KWHOLE, IDUM)
        CALL MDSTAT (NERR, MEM)
        IF (NERR .GT. 0) GOTO 100
      END IF

      CALL DBVTYP_BL (IVAR, TYP, ID)

      IF (TYP .EQ. 'H') THEN

        IF (HIST1) THEN

C         --Get history min/max memory

          CALL MDRSRV ('VALHMN', KVALHN, NVARHI)
          CALL MDRSRV ('ISTHMN', KISTHN, NVARHI)
          CALL MDRSRV ('VALHMX', KVALHX, NVARHI)
          CALL MDRSRV ('ISTHMX', KISTHX, NVARHI)
          CALL MDSTAT (NERR, MEM)
          IF (NERR .GT. 0) GOTO 100

C         --Initialize history min/max

          CALL SCAINI (0, NVARHI, A(KISTHN))

          HIST1 = .FALSE.
        END IF

C      --Get the offset for the variable

        IX = ID-1

C      --Determine if min/max is already calculated

        CALL SCACAL (NAME, IVAR, .FALSE., IELBST,
     &    A(KISTHN+IX), ICALC)

        IF (ICALC .GT. 0) THEN

C         --Reserve space for variable array
          CALL MDRSRV ('SCAVAR', KVAR, NVARHI)
          CALL MDSTAT (NERR, MEM)
          IF (NERR .GT. 0) GOTO 100

C         --Calculate min/max for all history variables

          CALL SCAHIS (A, A(KVAR), A(KWHOLE),
     &      A(KVALHN), A(KISTHN),
     &      A(KVALHX), A(KISTHX))

          CALL MDDEL ('SCAVAR')
        END IF

C      --Assign variable min/max

        VALMIN = A(KVALHN+IX)
        VALMAX = A(KVALHX+IX)

C      --Print min/max for history variable

        IF ((IPRINT .GT. 0) .AND.
     &    ((IPRINT .GT. 1) .OR. (ICALC .GT. 0))) THEN
          CALL SCAPRT (NAMECO, NAME, IVAR, A(KTIMES),
     &      A(KVALHN+IX), IDUM, RDUM, IA(KISTHN+IX),
     &      A(KVALHX+IX), IDUM, RDUM, IA(KISTHX+IX))
        END IF

      ELSE IF (TYP .EQ. 'G') THEN

        IF (GLOB1) THEN

C         --Get global min/max memory

          CALL MDRSRV ('VALGMN', KVALGN, NVARGL)
          CALL MDRSRV ('ISTGMN', KISTGN, NVARGL)
          CALL MDRSRV ('VALGMX', KVALGX, NVARGL)
          CALL MDRSRV ('ISTGMX', KISTGX, NVARGL)
          CALL MDSTAT (NERR, MEM)
          IF (NERR .GT. 0) GOTO 100

C         --Initialize global min/max

          CALL SCAINI (0, NVARGL, A(KISTGN))

          GLOB1 = .FALSE.
        END IF

C      --Get the offset for the variable

        IX = ID-1

C      --Determine if min/max is already calculated

        CALL SCACAL (NAME, IVAR, .FALSE., IELBST,
     &    A(KISTGN+IX), ICALC)

        IF (ICALC .GT. 0) THEN

C         --Reserve space for variable array
          CALL MDRSRV ('SCAVAR', KVAR, NVARGL)
          CALL MDSTAT (NERR, MEM)
          IF (NERR .GT. 0) GOTO 100

C         --Calculate min/max for all global variables

          CALL SCAGLO (A, A(KVAR), A(KWHOLE),
     &      A(KVALGN), A(KISTGN),
     &      A(KVALGX), A(KISTGX))

          CALL MDDEL ('SCAVAR')

        END IF

C      --Assign variable min/max

        VALMIN = A(KVALGN+IX)
        VALMAX = A(KVALGX+IX)

C      --Print min/max for global variable

        IF ((IPRINT .GT. 0) .AND.
     &    ((IPRINT .GT. 1) .OR. (ICALC .GT. 0))) THEN
          CALL SCAPRT (NAMECO, NAME, IVAR, A(KTIMES),
     &      A(KVALGN+IX), IDUM, RDUM, IA(KISTGN+IX),
     &      A(KVALGX+IX), IDUM, RDUM, IA(KISTGX+IX))
        END IF

      ELSE IF (TYP .EQ. 'N') THEN

        IF (NODE1) THEN

C         --Get nodal coordinates

          CALL MDFIND ('XN', KXN, IDUM)
          CALL MDFIND ('YN', KYN, IDUM)
          IF (NDIM .GE. 3) THEN
            CALL MDFIND ('ZN', KZN, IDUM)
          ELSE
            KZN = 1
          END IF
          CALL MDSTAT (NERR, MEM)
          IF (NERR .GT. 0) GOTO 100

C         --Get nodal min/max memory

          CALL MDRSRV ('VALNMN', KVALNN, NVARNP)
          CALL MDRSRV ('NUMNMN', KNUMNN, NVARNP)
          CALL MDRSRV ('XYZNMN', KXYZNN, 3*NVARNP)
          CALL MDRSRV ('ISTNMN', KISTNN, NVARNP)
          CALL MDRSRV ('VALNMX', KVALNX, NVARNP)
          CALL MDRSRV ('NUMNMX', KNUMNX, NVARNP)
          CALL MDRSRV ('XYZNMX', KXYZNX, 3*NVARNP)
          CALL MDRSRV ('ISTNMX', KISTNX, NVARNP)
          CALL MDSTAT (NERR, MEM)
          IF (NERR .GT. 0) GOTO 100

C         --Initialize nodal min/max

          CALL SCAINI (0, NVARNP, A(KISTNN))

          NODE1 = .FALSE.
        END IF

C      --Get the offset for the variable

        IX = ID-1
        IXX = 3*IX

C      --Determine if min/max is already calculated

        CALL SCACAL (NAME, IVAR, .FALSE., IELBST,
     &    A(KISTNN+IX), ICALC)

        IF (ICALC .GT. 0) THEN

C         --Reserve space for variable array
          CALL MDRSRV ('SCAVAR', KVAR, NUMNP)
          CALL MDSTAT (NERR, MEM)
          IF (NERR .GT. 0) GOTO 100

C         --Calculate min/max for nodal variable

          CALL SCANOD (A, IVAR, A(KVAR), A(KWHOLE),
     &      A(KXN), A(KYN), A(KZN),
     &      A(KVALNN+IX), A(KNUMNN+IX), A(KXYZNN+IXX),
     &      A(KISTNN+IX),
     &      A(KVALNX+IX), A(KNUMNX+IX), A(KXYZNX+IXX),
     &      A(KISTNX+IX))

          CALL MDDEL ('SCAVAR')

        END IF

C      --Assign variable min/max

        VALMIN = A(KVALNN+IX)
        VALMAX = A(KVALNX+IX)

C      --Print min/max for nodal variable

        IF ((IPRINT .GT. 0) .AND.
     &    ((IPRINT .GT. 1) .OR. (ICALC .GT. 0))) THEN
          CALL SCAPRT (NAMECO, NAME, IVAR, A(KTIMES),
     &      A(KVALNN+IX), MAPND(IA(KNUMNN+IX)), A(KXYZNN+IXX),
     &      IA(KISTNN+IX),
     &      A(KVALNX+IX), MAPND(IA(KNUMNX+IX)), A(KXYZNX+IXX),
     *      IA(KISTNX+IX))
        end if

      ELSE IF (TYP .EQ. 'E') THEN

        IF (ELEM1) THEN

C         --Get element coordinates and element blocks

          CALL MDFIND ('XE', KXE, IDUM)
          CALL MDFIND ('YE', KYE, IDUM)
          IF (NDIM .GE. 3) THEN
            CALL MDFIND ('ZE', KZE, IDUM)
          ELSE
            KZE = 1
          END IF

          CALL MDFIND ('LENE', KLENE, IDUM)
          CALL MDFIND ('ISEVOK', KIEVOK, IDUM)
          CALL MDSTAT (NERR, MEM)
          IF (NERR .GT. 0) GOTO 100

C         --Get element min/max memory

          L = (1+NELBLK) * NVAREL
          IXALIV = L
          CALL MDRSRV ('VALEMN', KVALEN, L+L)
          CALL MDRSRV ('NUMEMN', KNUMEN, L+L)
          CALL MDRSRV ('XYZEMN', KXYZEN, 3*(L+L))
          CALL MDRSRV ('ISTEMN', KISTEN, L+L)
          CALL MDRSRV ('VALEMX', KVALEX, L+L)
          CALL MDRSRV ('NUMEMX', KNUMEX, L+L)
          CALL MDRSRV ('XYZEMX', KXYZEX, 3*(L+L))
          CALL MDRSRV ('ISTEMX', KISTEX, L+L)
          CALL MDSTAT (NERR, MEM)
          IF (NERR .GT. 0) GOTO 100

C         --Initialize element min/max (birth/death done later)

          CALL SCAINI (NELBLK, NVAREL, A(KISTEN+0))

          ELEM1 = .FALSE.
        END IF

C      --Initialize alive element min/max, if needed

        IF (NALVAR .GT. 0) THEN
          IF (NALOLD .NE. NALVAR) THEN
            IX = 0
            IF (NALVAR .GT. 0) IX = IXALIV

            CALL SCAINI (NELBLK, NVAREL, A(KISTEN+IX))

            NALOLD = NALVAR
          END IF
        END IF

C      --Get the offset for the variable

        IX = (1+NELBLK) * (ID-1)
        IF (NALVAR .GT. 0) IX = IXALIV + IX
        IXX = 3*IX

C      --Determine if min/max is already calculated

        CALL SCACAL (NAME, IVAR, USESEL, IELBST,
     &    A(KISTEN+IX), ICALC)

        IF (ICALC .GT. 1) THEN

C         --Reserve space for variable and ALIVE array
          CALL MDRSRV ('SCAVAR', KVAR, NUMEL)
          CALL MDRSRV ('SCAALV', KALIVE, NUMEL)
          CALL MDSTAT (NERR, MEM)
          IF (NERR .GT. 0) GOTO 100

C         --Calculate min/max for element variable by element block

          CALL SCAELE (A, IVAR, A(KLENE), A(KIEVOK),
     &      NALVAR, A(KALIVE), A(KVAR), A(KWHOLE),
     &      A(KXE), A(KYE), A(KZE),
     &      A(KVALEN+IX), A(KNUMEN+IX), A(KXYZEN+IXX),
     &      A(KISTEN+IX),
     &      A(KVALEX+IX), A(KNUMEX+IX), A(KXYZEX+IXX),
     &      A(KISTEX+IX))

          CALL MDDEL ('SCAALV')
          CALL MDDEL ('SCAVAR')

        END IF

C         --Calculate min/max for selected element blocks

        CALL SCAELB (A, USESEL, IELBST,
     &    A(KVALEN+IX), A(KNUMEN+IX), A(KXYZEN+IXX), A(KISTEN+IX),
     &    A(KVALEX+IX), A(KNUMEX+IX), A(KXYZEX+IXX), A(KISTEX+IX),
     &    VALMIN, NUMMIN, XYZMIN, ISTMIN,
     &    VALMAX, NUMMAX, XYZMAX, ISTMAX)

C      --Assign variable min/max

        IF (ICALC .LE. 0) THEN
          VALMIN = A(KVALEN+IX)
          VALMAX = A(KVALEX+IX)
        END IF

C      --Print min/max for element variable

        IF ((IPRINT .GT. 0) .AND.
     &    ((IPRINT .GT. 1) .OR. (ICALC .GT. 0))) THEN
          CALL SCAPRT (NAMECO, NAME, IVAR, A(KTIMES),
     &      VALMIN, MAPEL(NUMMIN), XYZMIN, ISTMIN,
     &      VALMAX, MAPEL(NUMMAX), XYZMAX, ISTMAX)
        END IF
      END IF

 100  CONTINUE
      RETURN
      END
