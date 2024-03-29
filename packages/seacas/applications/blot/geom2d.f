C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE GEOM2D (LENF, NLNKF, LINKF, IF2EL,
     &   LENL, LINSET, IEBSET, LINDEF, NREF, LREF)
C=======================================================================

C   --*** GEOM2D *** (MESH) Group lines connecting nodes
C   --   Modified by John H. Glick - 10/25/88
C   --   Written by Amy Gilkey - revised 03/29/88
C   --   D. P. Flanagan, 07/02/82
C   --
C   --GEOM2D groups the lines connecting nodes.  Each unique line in the
C   --mesh is defined by searching through the connectivity.
C   --The line set is then sorted into the following parts:
C   --   -1) Mesh boundary
C   --    0) Element block boundary
C   --    n) Interior within Element block 'n'
C   --
C   --Parameters:
C   --   LENF - IN - the cumulative face counts by element block
C   --   NLNKF - IN - the number of nodes per face
C   --   LINKF - IN - the connectivity for all faces
C   --   IF2EL - IN - the element number of each face
C   --   LENL - OUT - the cumulative line counts by element block
C   --   LINSET - OUT - the sorted line set
C   --   IEBSET - SCRATCH - the element block for LINSET;
C   --      size = maximum number of lines
C   --   LINDEF - SCRATCH - the non-boundary line set;
C   --      size = (1+LLNSET) * maximum number of lines
C   --   NREF - SCRATCH - the number of references to a node in the line set;
C   --      size = NUMNPF
C   --   LREF - SCRATCH - the last line set index of a node in the line set;
C   --      size = NUMNPF
C   --
C   --Common Variables:
C   --   Uses NELBLK of /DBNUMS/
C   --   Uses NUMNPF, LLNSET of /D3NUMS/

      common /debugc/ cdebug
      common /debugn/ idebug
      character*8 cdebug

      include 'dbnums.blk'
      include 'd3nums.blk'

      INTEGER LENF(0:NELBLK)
      INTEGER NLNKF(NELBLK)
      INTEGER LINKF(*)
      INTEGER IF2EL(*)
      INTEGER LENL(-2:NELBLK), LINSET(LLNSET,*)
      INTEGER IEBSET(*)
      INTEGER LINDEF(0:LLNSET,*)
      INTEGER NREF(NUMNPF), LREF(NUMNPF)

C   --Define line set

      DO 100 IELB = -2, NELBLK
         LENL(IELB) = 0
  100 CONTINUE

      CALL INIINT (NUMNPF, 0, NREF)
      CALL INIINT (NUMNPF, 0, LREF)

      NHEAP = 0
      NSTACK = 0

      DO 150 IELB = 1, NELBLK

         IXL0 = IDBLNK (IELB, 0, LENF, NLNKF) - 1
         NNPF = NLNKF(IELB)
         IF (NNPF .EQ. 9)NNPF = 8
         DO 140 IFAC = LENF(IELB-1)+1, LENF(IELB)

C         --Extract element lines and element block type

            NOLD = NHEAP
            N2 = LINKF(IXL0+NNPF)
            LREFN2 = LREF(N2)
            LSTREF = LREFN2

            DO 130 ILINK = 1, NNPF
               N1 = N2
               N2 = LINKF(IXL0+ILINK)
               LREFN1 = LREFN2
               IF (ILINK .LT. NNPF) THEN
                  LREFN2 = LREF(N2)
               ELSE
                  LREFN2 = LSTREF
               END IF

C            --Determine if line already exists in line set

               DO 110 IL = MIN (LREFN1, LREFN2, NOLD), 1, -1
                  IF (LINSET(1,IL) .EQ. N2) THEN
                     IF (LINSET(2,IL) .EQ. N1) GOTO 120
                  END IF
  110          CONTINUE
  120          CONTINUE

               IF (IL .GT. 0) THEN
                  IELB1 = IEBSET(IL)

C               --Delete line from boundary set

                  IEBSET(IL) = IEBSET(NHEAP)
                  LINSET(1,IL) = LINSET(1,NHEAP)
                  LINSET(2,IL) = LINSET(2,NHEAP)
                  NREF(N1) = NREF(N1) - 1
                  NREF(N2) = NREF(N2) - 1
                  IF (NREF(N1) .LE. 0) LREF(N1) = 0
                  IF (NREF(N2) .LE. 0) LREF(N2) = 0
                  NHEAP = NHEAP - 1
                  IF (NOLD .GT. NHEAP) NOLD = NHEAP

C               --Add line to element block boundary or interior set

                  IF (IELB .EQ. IELB1) THEN
C                  --Move line to interior set
                     CONTINUE
                  ELSE
C                  --Move line to element block boundary set
                     IELB1 = 0
                  END IF

                  NSTACK = NSTACK + 1
                  LINDEF(0,NSTACK) = IELB1
                  LINDEF(1,NSTACK) = N2
                  LINDEF(2,NSTACK) = N1
                  LENL(IELB1) = LENL(IELB1) + 1

               ELSE

C               --Add line to boundary set

                  NHEAP = NHEAP + 1
                  IEBSET(NHEAP) = IELB
                  LINSET(1,NHEAP) = N1
                  LINSET(2,NHEAP) = N2
                  NREF(N1) = NREF(N1) + 1
                  NREF(N2) = NREF(N2) + 1
                  LREF(N1) = NHEAP
                  LREF(N2) = NHEAP
               END IF
  130       CONTINUE

            IXL0 = IXL0 + NLNKF(IELB)
  140    CONTINUE
  150 CONTINUE

C   --Sort element block boundary and interior line sets

      LENL(-1) = NHEAP
      DO 160 IELB = 0, NELBLK
         L = LENL(IELB)
         LENL(IELB) = NHEAP
         NHEAP = NHEAP + L
  160 CONTINUE
      DO 170 I = 1, NSTACK
         IELB = LINDEF(0,I)
         NHEAP = LENL(IELB) + 1
         LENL(IELB) = NHEAP
         CALL CPYINT (2, LINDEF(1,I), LINSET(1,NHEAP))
  170 CONTINUE

      RETURN
      END
