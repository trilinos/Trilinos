C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE SORLNX (NSETS, NEWELB, IE2ELB,
     &   NLNKF, LENSC, LINKSC, IF2ESC, IF22SC,
     &   LENF, LINKF, IF2EL, IF2EL2)
C=======================================================================

C   --*** SORLNX *** (MESH) Sort faces
C   --   Written by Amy Gilkey - revised 03/04/88
C   --
C   --SORLNX resorts all the faces into new element block sets.
C   --
C   --Parameters:
C   --   NSETS - IN - the number of element block sets in LINKF
C   --   NEWELB - IN - the new element block set for each face
C   --   IE2ELB - IN - the element for each element block
C   --   NLNKF - IN - the number of nodes for each element block
C   --   LENSC - IN - the cumulative face counts by element block for LINKSC
C   --   LINKSC - IN - the unsorted connectivity for all faces
C   --   IF2ESC - IN - the element number of each face in LINKSC
C   --   IF22SC - IN - the secondary element number of each face in LINKSC
C   --   LENF - OUT - the cumulative face counts by element block
C   --   LINKF - OUT - the connectivity for all faces
C   --   IF2EL - OUT - the element number of each face
C   --   IF2EL2 - OUT - the secondary element number of each face
C   --
C   --Common Variables:
C   --   Uses NELBLK of /DBNUMS/

      common /debugc/ cdebug
      common /debugn/ idebug
      character*8 cdebug

      include 'dbnums.blk'
      include 'd3nums.blk'

      INTEGER NEWELB(*)
      INTEGER IE2ELB(*)
      INTEGER NLNKF(NELBLK)
      INTEGER LENSC(0:NSETS)
      INTEGER LINKSC(*)
      INTEGER IF2ESC(*), IF22SC(*)
      INTEGER LENF(0:NSETS)
      INTEGER LINKF(*)
      INTEGER IF2EL(*), IF2EL2(*)

C   --Move faces into appropriate element block

      IX = 0
      IXL0 = 0
      DO 120 NELB = 1, NSETS

         IXLSC0 = 0
         DO 110 IELB = 1, NSETS
            IF (IELB .LE. NELBLK) NL = NLNKF(IELB)

            DO 100 IFAC = LENSC(IELB-1)+1, LENSC(IELB)
               IF (IELB .GT. NELBLK) NL = NLNKF(IE2ELB(IF2ESC(IFAC)))

               IF (NELB .EQ. IABS (NEWELB(IFAC))) THEN

                  IF (NEWELB(IFAC) .LT. 0) THEN

C                  --IF2ESC(IFAC) holds the "live" element number
                     I = IF2ESC(IFAC)
                     IF2ESC(IFAC) = IF22SC(IFAC)
                     IF22SC(IFAC) = I

C                  --Swap nodes to simulate surface being defined
C                  --by facing element
                     IF (NL .EQ. 4) THEN
                        I = LINKSC(IXLSC0+2)
                        LINKSC(IXLSC0+2) = LINKSC(IXLSC0+4)
                        LINKSC(IXLSC0+4) = I
                     END IF
                  END IF

                  IX = IX + 1
                  CALL CPYINT (NL, LINKSC(IXLSC0+1), LINKF(IXL0+1))
                  IXL0 = IXL0 + NL
                  IF2EL(IX) = IF2ESC(IFAC)
                  IF (IS3DIM) IF2EL2(IX) = IF22SC(IFAC)
               END IF

               IXLSC0 = IXLSC0 + NL
  100       CONTINUE

  110    CONTINUE
         LENF(NELB) = IX
  120 CONTINUE
      RETURN
      END
