C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE FIXAL3 (ALIVE, LENF, IF2EL, IF2EL2, IE2ELB, NEWELB)
C=======================================================================

C   --*** FIXAL3 *** (MESH) Adjust for element birth/death (3D)
C   --   Written by Amy Gilkey - revised 03/10/88
C   --
C   --FIXAL3 adjusts the face array to reflect the new element
C   --states.  Faces that are dead are moved to LENF(NELBLK+2).
C   --
C   --An element death causes all the faces of the element to be moved
C   --as follows:
C   --   Surface  -> dead
C   --   Interior -> surface
C   --An element birth causes all the faces of the element to be moved
C   --as follows:
C   --   Dead    -> surface (if one alive element for face) or
C   --              interior (if two alive elements)
C   --   Surface -> interior
C   --
C   --Parameters:
C   --   ALIVE - IN - true iff element i is alive
C   --   LENF - IN - the cumulative face counts by element block
C   --   IF2EL - IN - the element number of each face
C   --   IF2EL2 - IN - the secondary element number of each face
C   --   IE2ELB - IN - the element block for each element
C   --   NEWELB - OUT - size = LENF(NELBLK+1)
C   --
C   --Common Variables:
C   --   Uses NELBLK of /DBNUMS/

      include 'dbnums.blk'
      include 'd3nums.blk'

      LOGICAL ALIVE(NUMEL)
      INTEGER LENF(0:NELBLK+2)
      INTEGER IF2EL(*), IF2EL2(*)
      INTEGER IE2ELB(NUMEL)
      INTEGER NEWELB(*)

C   --Check each face by its defining elements, and move face if necessary

      DO 110 IELB = 1, NELBLK+2
         DO 100 IFAC = LENF(IELB-1)+1, LENF(IELB)

C         --Determine the number of elements which are alive for this face

            NALIVE = 0
            IF (ALIVE(IF2EL(IFAC))) NALIVE = NALIVE + 1
            IEL = IF2EL2(IFAC)
            IF (IEL .GT. 0) THEN
               IF (ALIVE(IEL)) NALIVE = NALIVE + 1
            END IF

            IF (NALIVE .EQ. 0) THEN

C            --If none alive, change to a DEAD face
               NEWELB(IFAC) = NELBLK+2

            ELSE IF (NALIVE .EQ. 1) THEN

C            --If only one alive, change to a SURFACE face

               IF (ALIVE(IF2EL(IFAC))) THEN
                  NEWELB(IFAC) = IE2ELB(IF2EL(IFAC))
               ELSE
C               --Swap nodes to simulate surface being defined by facing element
                  NEWELB(IFAC) = - IE2ELB(IF2EL2(IFAC))
               END IF

            ELSE IF (NALIVE .EQ. 2) THEN

C            --If both elements are alive, change to an INTERIOR face
               NEWELB(IFAC) = NELBLK+1

            END IF

  100    CONTINUE
  110 CONTINUE

      RETURN
      END
