C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE MUNESS (NUMESS, ISTAT, LESSEL, LESSDL,
     &  IDESS, NEESS, NEDSS, IXEESS, IXEDSS,
     &  LTEESS, LTSSS, FACSS,
     &  LTEX, LTSX, TDX, IXESS, IXDSS, NEX, NDX, ISCR, USESDF)
C=======================================================================

C   --*** MUNESS *** (GJOIN) Compress and rearrange element side sets
C   --   Written by Amy Gilkey - revised 02/25/88
C   --
C   --MUNESS processes the element side sets according to the set status.
C   --Sets may be combined or deleted.
C   --
C   --Parameters:
C   --   NUMESS - IN/OUT - the number of element side sets
C   --   ISTAT - IN - the status of each set:
C   --      0 = same
C   --      - = delete
C   --      n = combine with set n
C   --   LESSEL - IN/OUT - the length of the element side sets element list
C   --   LESSDL - IN/OUT - the length of the element side sets dist-fact list
C   --   IDESS - IN/OUT - the element side set ID for each set
C   --   NEESS - IN/OUT - the number of elements for each set
C   --   NEDSS - IN/OUT - the number of dist-fact for each set
C   --   IXEESS - IN/OUT - the index of the first element for each set
C   --   IXEDSS - IN/OUT - the index of the first dist-fact for each set
C   --   LTEESS - IN/OUT - the elements for all sets
C   --   LTSSS - IN/OUT - the sides for all sets
C   --   FACSS - IN/OUT - the dist-fact for all sets
C   --   LTEX - SCRATCH - sized to hold the set elements
C   --   LTSX - SCRATCH - sized to hold the set sides
C   --   TDX - SCRATCH - sized to hold the set dist-fact
C   --   IXESS - SCRATCH - size = NUMESS
C   --   IXDSS - SCRATCH - size = NUMESS -- dist-fact
C   --   NEX - SCRATCH - size = NUMESS
C   --   NDX - SCRATCH - size = NUMESS -- dist-face
C   --   ISCR - SCRATCH - size = NUMESS

      INTEGER ISTAT(*)
      INTEGER IDESS(*)
      INTEGER NEESS(*), NEDSS(*)
      INTEGER IXEESS(*)
      INTEGER LTEESS(*), LTEX(*)
      INTEGER LTSSS(*), LTSX(*)
      INTEGER IXESS(*), IXEDSS(*)
      INTEGER NEX(*), NDX(*)
      INTEGER ISCR(*)
      REAL    FACSS(*), TDX(*)
      LOGICAL USESDF

      IF (NUMESS .LE. 0) RETURN

      JESS = 0
      JNE = 1
      JND = 1
      DO 110 IESS = 1, NUMESS

         IF (ISTAT(IESS) .EQ. 0) THEN
            NINSET = 1
            ISCR(NINSET) = IESS
         ELSE IF (ISTAT(IESS) .EQ. IESS) THEN
            CALL GETALL (IESS, NUMESS, ISTAT, NINSET, ISCR)
         ELSE
            NINSET = 0
         END IF

         IF (NINSET .GT. 0) THEN
            JESS = JESS + 1
            IXESS(JESS) = IESS
            NEX(JESS) = 0
            NDX(JESS) = 0
         END IF
         DO 100 ISET = 1, NINSET
            N = ISCR(ISET)
            CALL MOVINT (NEESS(N), LTEESS(IXEESS(N)), LTEX(JNE))
            CALL MOVINT (NEESS(N), LTSSS(IXEESS(N)),  LTSX(JNE))
            IF (USESDF) THEN
              CALL MOVREA (NEDSS(N), FACSS(IXEDSS(N)),  TDX(JND))
            ENDIF
            JNE = JNE + NEESS(N)
            JND = JND + NEDSS(N)
            NEX(JESS) = NEX(JESS) + NEESS(N)
            NDX(JESS) = NDX(JESS) + NEDSS(N)
  100    CONTINUE
  110 CONTINUE

      CALL ORDIX (JESS, IXESS, NUMESS, IDESS, ISCR, IDESS)
      CALL MOVINT (JESS, NEX, NEESS)
      CALL MOVINT (JESS, NDX, NEDSS)
      NUMESS = JESS
      JNE = 1
      JND = 1
      DO 120 IESS = 1, NUMESS
         IXEESS(IESS) = JNE
         IXEDSS(IESS) = JND
         JNE = JNE + NEESS(IESS)
         JND = JND + NEDSS(IESS)
  120 CONTINUE
      LESSEL = JNE - 1
      LESSDL = JND - 1
      CALL MOVINT (LESSEL, LTEX, LTEESS)
      CALL MOVINT (LESSEL, LTSX, LTSSS)
      if (usesdf) then
        CALL MOVREA (LESSDL, TDX,  FACSS)
      end if
      RETURN
      END
