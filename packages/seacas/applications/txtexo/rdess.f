C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE RDESS (NTXT, NUMESS, LESSEL, LESSNL, LESSDF,
     &  IDESS, NEESS, NNESS, NDESS, IXEESS, IXDESS, LTEESS, LTSESS,
     *  FACESS, *)
C=======================================================================

C   --*** RDESS *** (TXTEXO) Read database side sets
C   --   Written by Amy Gilkey - revised 02/27/86
C   --
C   --RDESS reads the side set information from the text file.
C   --An error message is displayed if the end of file is read.
C   --
C   --Parameters:
C   --   NTXT - IN - the text file
C   --   NUMESS - IN - the number of side sets
C   --   LESSEL - IN - the length of the side sets element list
C   --   LESSNL - IN - the length of the side sets node list
C   --   IDESS - OUT - the side set ID for each set
C   --   NEESS - OUT - the number of elements for each set
C   --   NNESS - OUT - the number of nodes for each set
C   --   NDESS - OUT - the number of dist factors for each set
C   --   IXEESS - OUT - the index of the first element for each set
C   --   IXDESS - OUT - the index of the first node for each set
C   --   LTEESS - OUT - the elements for all sets
C   --   LTSESS - OUT - the element sides for all sets
C   --   FACESS - OUT - the distribution factors for all sets
C   --   * - return statement if end of file or read error
C   --
C   --Database must be positioned at start of side set information
C   --upon entry; upon exit at end of side set information.

      INTEGER IDESS(*)
      INTEGER NEESS(*)
      INTEGER NNESS(*)
      INTEGER NDESS(*)
      INTEGER IXEESS(*)
      INTEGER IXDESS(*)
      INTEGER LTEESS(*)
      INTEGER LTSESS(*)
      REAL FACESS(*)

      CHARACTER*32 STRA, STRB

C ... Skip comment records
      READ (NTXT, *, END=120, ERR=120)
      READ (NTXT, *, END=120, ERR=120)

      NE = 0
      NN = 0
      ND = 0
      DO 110 IESS = 1, NUMESS
         READ (NTXT, *, END=120, ERR=120) IDESS(IESS), NEESS(IESS),
     &    NNESS(IESS), NDESS(IESS)
C ... Skip comment record
         READ (NTXT, *, END=120, ERR=120)

         IXEESS(IESS) = NE + 1
         NE = NE + NEESS(IESS)
         if (neess(iess) .gt. 0) then
           READ (NTXT, *, END=130, ERR=130)
     &       (LTEESS(NL), LTSESS(NL), NL=IXEESS(IESS), NE)
         end if

C ... Skip comment record
         READ (NTXT, *, END=120, ERR=120)
         IXDESS(IESS) = ND + 1
         ND = ND + NDESS(IESS)
         NN = NN + NNESS(IESS)
         DO 100 NL = IXDESS(IESS), ND
            READ (NTXT, *, END=140, ERR=140) LTNESS, FACESS(NL)
  100    CONTINUE
  110 CONTINUE

      IF (NE .NE. LESSEL) THEN
         CALL INTSTR (1, 0, NE, STRA, LSTRA)
         CALL INTSTR (1, 0, LESSEL, STRB, LSTRB)
         CALL PRTERR ('WARNING',
     &      'SIDE SETS NUMBER OF ELEMENTS = ' // STRA(:LSTRA)
     &      // ' does not match TOTAL = ' // STRB(:LSTRB))
      END IF

      IF (NN .NE. LESSNL) THEN
         CALL INTSTR (1, 0, NN, STRA, LSTRA)
         CALL INTSTR (1, 0, LESSNL, STRB, LSTRB)
         CALL PRTERR ('WARNING',
     &      'SIDE SETS NUMBER OF NODES = ' // STRA(:LSTRA)
     &      // ' does not match TOTAL = ' // STRB(:LSTRB))
      END IF

      IF (ND .NE. LESSDF) THEN
         CALL INTSTR (1, 0, ND, STRA, LSTRA)
         CALL INTSTR (1, 0, LESSDF, STRB, LSTRB)
         CALL PRTERR ('WARNING',
     &      'SIDE SETS DISTRIBUTION FACTOR COUNT = ' // STRA(:LSTRA)
     &      // ' does not match TOTAL = ' // STRB(:LSTRB))
      END IF

      RETURN

  120 CONTINUE
      CALL INTSTR (1, 0, IESS, STRA, LSTRA)
      CALL PRTERR ('FATAL',
     &   'Reading HEADER DATA for SIDE SET ' // STRA(:LSTRA))
      GOTO 150
  130 CONTINUE
      CALL INTSTR (1, 0, IESS, STRA, LSTRA)
      CALL PRTERR ('FATAL',
     &   'Reading ELEMENTS for SIDE SET ' // STRA(:LSTRA))
      GOTO 150
  140 CONTINUE
      CALL INTSTR (1, 0, NL-IXDESS(IESS)+1, STRA, LSTRA)
      CALL INTSTR (1, 0, IESS, STRB, LSTRB)
      CALL PRTERR ('FATAL',
     &   'Reading NODES and FACTORS for node ' // STRA(:LSTRA)
     &   // ' for SIDE SET ' // STRB(:LSTRB))
      GOTO 150
  150 CONTINUE
      RETURN 1
      END
