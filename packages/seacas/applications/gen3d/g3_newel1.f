C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE NEWEL1 (BLKTYP, NUMELB, NELB3, NRSTRT, NREND,
     &   NUMLNK, NUMLN3, NUMATR, NATRDM, NUMNP, NUMNP3,
     &   LINK, LINK3, ATRIB, ATRIB3,
     &   IXEL, INCEL, NREL, IELCOL, IXNP, NRNP)
C=======================================================================

C   --*** NEWEL1 *** (GEN3D) Write 3D element block misc.
C   --   Written by Amy Gilkey - revised 05/05/86
C   --
C   --NEWEL1 calculates and writes the element block connectivity and
C   --attributes for the 3D database.  It may write only a section of the
C   --element block (depending on NRSTRT and NREND).
C   --
C   --Parameters:
C   --   BLKTYP - IN - the element block type
C   --   NUMELB - IN - the number of elements in the block
C   --   NELB3 - IN - the number of elements in the 3D block
C   --   NRSTRT, NREND - IN - the starting and ending translation/rotation
C   --   NUMLNK - IN - the number of nodes in the element (always 4)
C                       Note: skip 3,4 for bars
C   --   NUMLN3 - IN - the number of nodes in the 3D element
C   --   NUMATR - IN - the number of attributes
C   --   NUMNP - IN - the number of nodes
C   --   NUMNP3 - IN - the number of 3D nodes
C   --   LINK - IN - the 2D connectivity array for the block
C   --   LINK3 - OUT - the 3D connectivity array (packed) for the block section
C   --   ATRIB - IN - the 2D attribute array (packed) for the block
C   --   ATRIB3 - OUT - the 3D attribute array (packed) for the block section
C   --   IXEL - IN - the new index for each element
C   --   INCEL - IN - the increment for each element, needed for blocks
C   --      that become multiple blocks
C   --   NREL - IN - the number of new elements generated for each element
C   --   IELCOL - IN - the row number for each element, 0 if not needed
C   --   IXNP - IN - the new index for each node
C   --   NRNP - IN - the number of new nodes generated for each node
C   --
C   --Common Variables:
C   --   Uses NDBOUT of /DBASE/
C   --   Uses NNREPL, NEREPL, DIM3 of /PARAMS/
C   --   Uses IX1, IX2, IX3, IX4 of /CENPAR/

      INCLUDE 'g3_dbase.blk'
      INCLUDE 'g3_params.blk'
      INCLUDE 'g3_cenpar.blk'

      CHARACTER*(*) BLKTYP
      INTEGER LINK(4,NUMELB)
      INTEGER LINK3(NUMLN3,NELB3)
      REAL ATRIB(NATRDM,NUMELB)
      REAL ATRIB3(NATRDM,NELB3)
      INTEGER IXEL(*), INCEL(*), NREL(*), IELCOL(*)
      INTEGER IXNP(*), NRNP(*)

      LOGICAL ROT360

      ROT360 = (NEREPL .EQ. NNREPL)

      N4 = NINT (DIM3 / 90)

C   --Connectivity - add on nodes for next plane/slice for bottom face
C   --(1-2-3-4-5-6-7-8), repeat for each plate/slice

C   --Find smallest element offset for block

      JELOFF = IXEL(1)
      DO IEL = 1, NUMELB
        JELOFF = MIN (JELOFF, IXEL(IEL))
      end do
      JELOFF = JELOFF - 1

      DO IEL = 1, NUMELB
         JEL = IXEL(IEL) - JELOFF

         IF (IELCOL(IEL) .EQ. 0) THEN

C         --Handle non-center elements

            DO NR = NRSTRT, NREND
              IF (NUMLNK .EQ. 2) THEN
C ... Bars to shells
                INP1 = LINK(2,IEL)
                INP2 = LINK(1,IEL)
                JNP1 = IXNP(INP1)
                JNP2 = IXNP(INP2)
                IF (NR .LT. NNREPL) THEN
                  LINK3(1,JEL) = JNP1 + NR
                  LINK3(2,JEL) = JNP1 + NR - 1
                  LINK3(3,JEL) = JNP2 + NR - 1
                  LINK3(4,JEL) = JNP2 + NR
                ELSE
                  LINK3(1,JEL) = JNP1
                  LINK3(2,JEL) = JNP1 + NR - 1
                  LINK3(3,JEL) = JNP2 + NR - 1
                  LINK3(4,JEL) = JNP2
                END IF
              ELSE
                DO J = 1, NUMLNK
                  INP = LINK(J,IEL)
                  JNP = IXNP(INP)
                  IF (NR .LT. NNREPL) THEN
                    LINK3(J,JEL) = JNP + NR
                  ELSE
                    LINK3(J,JEL) = JNP
                  END IF
                  LINK3(J+NUMLNK,JEL) = JNP + NR-1
                end do
              END IF
              JEL = JEL + INCEL(IEL)
            end do

         ELSE

C         --Handle center element, different for corner elements
C         --NOTE: assumes LINK(IX1,i) and LINK(IX4,i) in same column
C         --   and LINK(IX2,i) and LINK(IX3,i) in next column

            IF (NUMLNK .NE. 4) THEN
               CALL PRTERR('FATAL',
     $              'Option not implemented in NEWEL1')
               STOP 'Unimplemented Option'
            END IF
            NE4 = NREL(IEL) / N4
            NCORN = INT (NE4/2) + 1
            INP1 = LINK(IX1,IEL)
            JNP1 = IXNP(INP1)
            INP2 = LINK(IX2,IEL)
            JNP2 = IXNP(INP2)
            INP3 = LINK(IX3,IEL)
            JNP3 = IXNP(INP3)
            INP4 = LINK(IX4,IEL)
            JNP4 = IXNP(INP4)
            I1 = 0
            I2 = 0
            DO NR = 1, NREL(IEL)
               IF (NR .EQ. NCORN) THEN
                  I3 = I2 + 2
                  IF (ROT360 .AND. I3 .GE. NRNP(INP2)) I3 = 0
                  LINK3(IX1,JEL) = JNP2 + I3
                  LINK3(IX4,JEL) = JNP3 + I3
               ELSE
                  I3 = I1 + 1
                  IF (ROT360 .AND. I3 .GE. NRNP(INP1)) I3 = 0
                  LINK3(IX1,JEL) = JNP1 + I3
                  LINK3(IX4,JEL) = JNP4 + I3
               END IF
               I4 = I2 + 1
               IF (ROT360 .AND. I4 .GE. NRNP(INP2)) I4 = 0
               LINK3(IX2,JEL) = JNP2 + I4
               LINK3(IX3,JEL) = JNP3 + I4
               LINK3(IX1+4,JEL) = JNP1 + I1
               LINK3(IX2+4,JEL) = JNP2 + I2
               LINK3(IX3+4,JEL) = JNP3 + I2
               LINK3(IX4+4,JEL) = JNP4 + I1
               IF (NR .EQ. NCORN) THEN
                  NCORN = NCORN + NE4
                  I1 = I1
                  I2 = I2 + 2
               ELSE
                  I1 = I1 + 1
                  I2 = I2 + 1
               END IF
               JEL = JEL + 1
             end do
         END IF
       end do

C   --Attributes - repeat attributes for the next plates/slices

      DO IEL = 1, NUMELB
        JEL = IXEL(IEL) - JELOFF
        DO NR = NRSTRT, NREND
          CALL CPYREA (NUMATR, ATRIB(1,IEL), ATRIB3(1,JEL))
          JEL = JEL + INCEL(IEL)
        end do
      end do

      RETURN
      END
