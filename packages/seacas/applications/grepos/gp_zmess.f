C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE ZMESS (NUMESS, ISTAT, LESSEL, LESSDL, IDESS, NEESS,
     *     NEDSS, IXEESS, IXEDSS, LTEESS, LTSSS, LTSNC, FAC, NAMES)
C=======================================================================

C     --*** ZMESS *** (GJOIN) Compress element side sets
C     --   Written by Amy Gilkey - revised 01/20/88
C     --
C     --ZMESS compresses the element side sets by removing deleted elements
C     --and their nodes.  Assumes that the elements and nodes are already
C     --renumbered and ordered.
C     --
C     --Parameters:
C     --   NUMESS - IN/OUT - the number of element side sets
C     --   LESSEL - IN/OUT - the length of the element side sets element list
C     --   IDESS - IN/OUT - the element side set ID for each set
C     --   NEESS - IN/OUT - the number of elements for each set
C     --   NEDSS - IN/OUT - the number of dist-fac for each set
C     --   IXEESS - IN/OUT - the index of the first element for each set
C     --   IXEDSS - IN/OUT - the index of the first dist-fac for each set
C     --   LTEESS - IN/OUT - the elements for all sets
C     --   LTSSS - IN/OUT - the sides for all sets
C     --   LTSNC - IN/OUT - the face count for each element/side in the list
C     --   FACESS - IN/OUT - the distribution factors for all sets????????????

      include 'gp_namlen.blk'

      INTEGER ISTAT(*)          ! NUMESS
      INTEGER IDESS(*)          ! NUMESS
      INTEGER NEESS(*)          ! NUMESS
      INTEGER NEDSS(*)          ! NUMESS
      INTEGER IXEESS(*)         ! NUMESS
      INTEGER IXEDSS(*)         ! NUMESS
      INTEGER LTEESS(*)         ! LESSEL
      INTEGER LTSSS(*)          ! LESSEL
      INTEGER LTSNC(*)          ! LESSEL
      REAL    FAC(*)            ! LESSDL
      CHARACTER*(maxnam) NAMES(*)

      IF (NUMESS .LE. 0) RETURN

      JESS = 0
      JNE = 0
      JDF = 0
      idfe = 0
      DO 120 IESS = 1, NUMESS
         JNELST = JNE
         JDFLST = JDF
         nd1 = ixedss(iess)     ! Index of First distribution factor for this list
         IDFB = ND1
         DO 110 N = 1, NEESS(IESS)
C     ... N     is the 'local' index within the current set.
C     NE    is the 'global' index within the concatenated (LTEESS, LTSSS, LTSNC) lists
C     ND1   is the FAC index of the first df for the current list.
C     ICNT  is the number of DF for element N in the current list
            NE = N + IXEESS(IESS) - 1
            ICNT = LTSNC(NE)
C     IDFB = index of first df for local element N, global element NE
C     IDFE = index of last  df for local element N, global element NE
            IDFE = IDFB + ICNT - 1
            IF (LTEESS(NE) .GT. 0) THEN
               JNE = JNE + 1
               LTEESS(JNE) = LTEESS(NE)
               LTSSS(JNE)  = LTSSS(NE)
               LTSNC(JNE)  = LTSNC(NE)
               do 100 nd = idfb, idfe
                  JDF = JDF + 1
                  fac(JDF) = fac(ND)
 100           continue
            END IF
            IDFB = IDFE + 1
 110     CONTINUE
         N = JNE - JNELST
         IF (N .GT. 0) THEN
C     ... There is at least 1 element remaining in the list...
            JESS = JESS + 1     ! increment sideset count
            IDESS(JESS) = IDESS(IESS) ! copy the sideset id
            NAMES(JESS) = NAMES(IESS)
            NEESS(JESS) = N     ! Set the elements per list count
            IXEESS(JESS) = JNELST + 1 ! set the index
            NEDSS(JESS) = JDF - JDFLST ! set the df per list count
            IXEDSS(JESS) = JDFLST + 1
         ELSE
            ISTAT(IESS) = -IDESS(IESS)
         END IF
 120  CONTINUE
      if (idfe .ne. lessdl) stop 'ZMESS: Internal error'
      NUMESS = JESS
      LESSEL = JNE
      LESSDL = JDF
      RETURN
      END
