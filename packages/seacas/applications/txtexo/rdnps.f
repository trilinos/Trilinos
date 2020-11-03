C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE RDNPS (NTXT, NUMNPS, LNPSNL, LNPSDF,
     &   IDNPS, NNNPS, NDNPS, IXNNPS, IXDNPS, LSTNPS, FACNPS, *)
C=======================================================================

C   --*** RDNPS *** (TXTEXO) Read database node sets
C   --   Written by Amy Gilkey - revised 02/27/86
C   --
C   --RDNPS reads the node set information from the text file.
C   --An error message is displayed if the end of file is read.
C   --
C   --Parameters:
C   --   NTXT - IN - the text file
C   --   NUMNPS - IN - the number of node sets
C   --   LNPSNL - IN - the length of the node sets node list
C   --   IDNPS - OUT - the node set ID for each set
C   --   NNNPS - OUT - the number of nodes for each set
C   --   IXNNPS - OUT - the index of the first node for each set
C   --   LSTNPS - OUT - the nodes for all sets
C   --   FACNPS - OUT - the distribution factors for all sets
C   --   * - return statement if end of file or read error
C   --
C   --Database must be positioned at start of node set information
C   --upon entry; upon exit at end of node set information.

      INTEGER IDNPS(*)
      INTEGER NNNPS(*)
      INTEGER NDNPS(*)
      INTEGER IXNNPS(*)
      INTEGER IXDNPS(*)
      INTEGER LSTNPS(*)
      REAL FACNPS(*)

      CHARACTER*32 STRA, STRB

      NN = 0
      ND = 0
      READ (NTXT, *, END=120, ERR=120)
      READ (NTXT, *, END=120, ERR=120)
      READ (NTXT, *, END=120, ERR=120)
      DO INPS = 1, NUMNPS
         READ (NTXT, *, END=120, ERR=120)
         READ (NTXT, *, END=120, ERR=120) IDNPS(INPS), NNNPS(INPS),
     &        NDNPS(INPS)

         IXNNPS(INPS) = NN + 1
         NN = NN + NNNPS(INPS)

         IXDNPS(INPS) = ND + 1
         ND = ND + NDNPS(INPS)

         if (NDNPS(INPS) .NE. 0) THEN
C ... Node set has distribution factors
            IND = IXDNPS(INPS)
            DO NL = IXNNPS(INPS), NN
               READ (NTXT, *, END=130, ERR=130) LSTNPS(NL), FACNPS(IND)
               IND = IND + 1
            end do
         else
C ... Node set doesn't have distribution factors
            DO NL = IXNNPS(INPS), NN
               READ (NTXT, *, END=130, ERR=130) LSTNPS(NL)
            end do
         endif
      end do

      IF (NN .NE. LNPSNL) THEN
         CALL INTSTR (1, 0, NN, STRA, LSTRA)
         CALL INTSTR (1, 0, LNPSNL, STRB, LSTRB)
         CALL PRTERR ('WARNING',
     &      'NODE SET NUMBER OF NODES = ' // STRA(:LSTRA)
     &      // ' does not match TOTAL = ' // STRB(:LSTRB))
      END IF

      RETURN

  120 CONTINUE
      CALL INTSTR (1, 0, INPS, STRA, LSTRA)
      CALL PRTERR ('FATAL',
     &   'Reading HEADER DATA for NODE SET ' // STRA(:LSTRA))
      GOTO 140
  130 CONTINUE
      CALL INTSTR (1, 0, NL-IXNNPS(INPS)+1, STRA, LSTRA)
      CALL INTSTR (1, 0, INPS, STRB, LSTRB)
      CALL PRTERR ('FATAL',
     &   'Reading NODES and FACTORS for node ' // STRA(:LSTRA)
     &   // ' for NODE SET ' // STRB(:LSTRB))
      GOTO 140
  140 CONTINUE
      RETURN 1
      END
