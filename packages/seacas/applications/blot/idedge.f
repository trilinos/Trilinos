C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE IDEDGE (IEDSET, NEDGES, HIDENP, LINSET, LENL, IXLIN,
     *  MREF, LREF)
C=======================================================================

C   --*** IDEDGE *** (MESH) Identify 3D lines on edge of visible mesh
C   --   Written by Amy Gilkey - revised 10/02/86
C   --
C   --IDEDGE finds edges from IEDSET in the line set and marks them as
C   --"must draw" with a -n in LINSET(3,x).
C   --
C   --Parameters:
C   --   IEDSET - IN - the edge line set;
C   --      (0) = face defining edge; 0 to delete edge
C   --   NEDGES - IN - the number of lines in the edge set
C   --   HIDENP - IN/OUT - node status (as in HIDDEN)
C   --   LINSET - IN/OUT - the sorted line set
C   --   LENL - IN - the length of the line set
C   --   IXLIN - SCRATCH - index line set to be searched, length = LENL

      PARAMETER (KNVIS=0, KNFOVR=10, KNHID=100)

      common /debugc/ cdebug
      common /debugn/ idebug
      character*8 cdebug

      COMMON /D3NUMS/ IS3DIM, NNPSUR, NUMNPF, LLNSET
      LOGICAL IS3DIM

      INTEGER IEDSET(0:2,*)
      INTEGER HIDENP(*)
      INTEGER LINSET(LLNSET,*)
      INTEGER IXLIN(*)
      INTEGER MREF(*), LREF(*)

C   --Get indices of line set to search and order line set nodes

      LENIX = 0
      DO 100 IL = 1, LENL
         IF (LINSET(3,IL) .GE. 1) THEN
            LENIX = LENIX + 1
            IXLIN(LENIX) = IL
            IF (LINSET(3,IL) .EQ. 1) THEN
               IF (LINSET(1,IL) .GT. LINSET(2,IL)) THEN
                  IF (HIDENP(LINSET(1,IL)) .GT. KNVIS) THEN
                     I = LINSET(1,IL)
                     LINSET(1,IL) = LINSET(2,IL)
                     LINSET(2,IL) = I
                  END IF
               END IF
            END IF
         END IF
  100 CONTINUE

C ... For each node, mref is first linset edge referencing that node
C                    lref is last        linset edge referencing that node
      call iniint(numnpf, 0, mref)
      call iniint(numnpf, 0, lref)
      DO IX=LENIX, 1, -1
        IL = IXLIN(IX)
        MREF(LINSET(1,IL)) = IX
        MREF(LINSET(2,IL)) = IX
      END DO
      DO IX=1,LENIX
        IL = IXLIN(IX)
        LREF(LINSET(1,IL)) = IX
        LREF(LINSET(2,IL)) = IX
      END DO

      DO 130 IEDG = 1, NEDGES
         IF (IEDSET(0,IEDG) .EQ. 0) GOTO 130

C      --Make sure edge set is still ordered (some nodes may have become
C      --visible)

         N1 = IEDSET(1,IEDG)
         N2 = IEDSET(2,IEDG)
         IF (N1 .GT. N2) THEN
            IF (HIDENP(N1) .GT. KNVIS) THEN
               I = N1
               N1 = N2
               N2 = I
            END IF
         END IF

C      --Find edge in ordered line set

         imax = min(lref(n1), lref(n2), lenix)
         imin = max(mref(n1), mref(n2), 1)
         DO 110 IX = imin, imax
            IL = IXLIN(IX)
            IF (LINSET(1,IL) .EQ. N1 .and.
     *          LINSET(2,IL) .EQ. N2) THEN
C      --Mark line set as an edge
                LINSET(3,IL) = - LINSET(3,IL)
                GO TO 130
              END IF
  110    CONTINUE
  130 CONTINUE

      RETURN
      END
