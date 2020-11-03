C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE ZESS(IDESS, NEESS, NNESS, IXEESS, IXNESS, LTEESS,
     &   LTNESS, FACESS, NUMESS, LESSEL, LESSNL, MAXVAL)
C=======================================================================

C   --   NUMESS - IN - the number of side sets
C   --   LESSEL - IN - the length of the side sets element list
C   --   LESSNL - IN - the length of the side sets node list
C   --   IDESS - OUT - the side set ID for each set
C   --   NEESS - OUT - the number of elements for each set
C   --   NNESS - OUT - the number of nodes for each set
C   --   IXEESS - OUT - the index of the first element for each set
C   --   IXNESS - OUT - the index of the first node for each set
C   --   LTEESS - OUT - the elements for all sets
C   --   LTNESS - OUT - the nodes for all sets
C   --   FACESS - OUT - the distribution factors for all sets

C ... All deleted node and element numbers have been set equal to MAXVAL,
C     all active nodes and elements have already been remapped.
C     This routine only compresses the node and element lists.

      INTEGER NUMESS, LESSEL, LESSNL, MAXVAL
      INTEGER IDESS(*)
      INTEGER NEESS(*)
      INTEGER NNESS(*)
      INTEGER IXEESS(*)
      INTEGER IXNESS(*)
      INTEGER LTEESS(*)
      INTEGER LTNESS(*)
      REAL    FACESS(*)

      IOFF = 0
      IOFFN= 0
      ISET = 0
      DO 30 I=1, NUMESS
         IBEG = IXEESS(I)
         IEND = IBEG + NEESS(I) - 1
         IF (IDESS(I) .NE. 0) THEN
C ... This set may have been partially or totally deleted
C ... Count number of deleted elements (element number = MAXVAL)
            NUMDEL = NUMEQI (MAXVAL, NEESS(I), LTEESS(IBEG))
            IF (NUMDEL .EQ. NEESS(I)) THEN
C ... This set has been totally deleted
               IDESS(I) = 0
            ELSE
C ... This set has been partially deleted, NEESS(I)-NUMDEL nodes left
               NPSET = NNESS(I) / NEESS(I)
               ISET = ISET + 1
               IDESS(ISET)  = IDESS(I)
               NEESS(ISET)  = NEESS(I) - NUMDEL
               NNESS(ISET)  = NNESS(I) - NPSET * NUMDEL
               IBEGN        = IXNESS(I)
               IXEESS(ISET) = IOFF + 1
               IXNESS(ISET) = IOFFN + 1

               DO 20 IEL  = IBEG, IEND
                  IF (LTEESS(IEL) .LT. MAXVAL) THEN
C ... This element not deleted, copy it and its nodes to new lists
                     IOFF = IOFF + 1
                     LTEESS(IOFF) = LTEESS(IEL)
                     IBNOD = IBEGN + NPSET * (IEL - IBEG)
                     IENOD = IBNOD + NPSET - 1
                     DO 10 INOD = IBNOD, IENOD
                        IOFFN = IOFFN + 1
                        LTNESS(IOFFN) = LTNESS(INOD)
                        FACESS(IOFFN) = FACESS(INOD)
   10                CONTINUE
                  END IF
   20          CONTINUE
            END IF
         END IF
         IF (IDESS(I) .EQ. 0) THEN
C ... This set has been totally deleted
         END IF
   30 CONTINUE

      NUMESS = ISET
      LESSEL = IOFF
      LESSNL = IOFFN

      RETURN
      END
