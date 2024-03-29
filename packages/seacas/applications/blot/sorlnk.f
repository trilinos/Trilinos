C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE SORLNK (A, NSETS, NEWELB,
     &   LENF, NLNKF, LINKF, IF2EL, IF2EL2, IE2ELB)
C=======================================================================

C   --*** SORLNK *** (MESH) Sort faces
C   --   Written by Amy Gilkey - revised 03/10/88
C   --
C   --SORLNK resorts all the faces into new element block sets.  This
C   --is done by coping the faces, then renaming the arrays.
C   --
C   --This routine uses MDFIND to find the following dynamic memory arrays:
C   --   LINKF - for length only
C   --   IF2EL - for length only
C   --
C   --Parameters:
C   --   A - IN - the dynamic memory base array
C   --   NSETS - IN - the number of element block sets in LINKF
C   --   NEWELB - IN - the new element block set for each face
C   --   LENF - IN/OUT - the cumulative face counts by element block for LINKF
C   --   NLNKF - IN - the number of nodes for each element block
C   --   LINKF - IN/OUT - the unsorted connectivity for all faces
C   --   IF2EL - IN/OUT - the element number of each face in LINKF
C   --   IF2EL2 - IN/OUT - the secondary element number of each face in LINKF
C   --   IE2ELB - IN - the element for each element block
C   --
C   --Common Variables:
C   --   Uses NELBLK of /DBNUMS/

      common /debugc/ cdebug
      common /debugn/ idebug
      character*8 cdebug

      include 'dbnums.blk'
      include 'd3nums.blk'

      DIMENSION A(*)
      INTEGER NEWELB(*)
      INTEGER LENF(0:NSETS)
      INTEGER NLNKF(NELBLK)
      INTEGER LINKF(*)
      INTEGER IF2EL(*), IF2EL2(*)
      INTEGER IE2ELB(*)

      LOGICAL SOMMOV

C   --Check that some faces will move

      SOMMOV = .FALSE.
      DO 110 IELB = 1, NSETS
         DO 100 IFAC = LENF(IELB-1)+1, LENF(IELB)
            IF (IELB .NE. IABS (NEWELB(IFAC))) THEN
               SOMMOV = .TRUE.
               GOTO 120
            END IF
  100    CONTINUE
  110 CONTINUE
  120 CONTINUE

      IF (SOMMOV) THEN
         CALL MDRSRV ('XLENF', KLENF, 1+NSETS)
         CALL MDFIND ('LINKF', IDUM, LENLNK)
         CALL MDRSRV ('XLINKF', KLINKF, LENLNK)
         CALL MDFIND ('IF2EL', IDUM, LENF2E)
         CALL MDRSRV ('XIF2EL', KIF2EL, LENF2E)
         IF (IS3DIM) CALL MDRSRV ('XIF2E2', KIF2E2, LENF2E)
         CALL MDSTAT (NERR, MEM)
         IF (NERR .GT. 0) GOTO 130

         CALL CPYINT (1+NSETS, LENF, A(KLENF))
         CALL CPYINT (LENLNK, LINKF, A(KLINKF))
         CALL CPYINT (LENF2E, IF2EL, A(KIF2EL))
         IF (IS3DIM) CALL CPYINT (LENF2E, IF2EL2, A(KIF2E2))

         CALL SORLNX (NSETS, NEWELB, IE2ELB,
     &      NLNKF, A(KLENF), A(KLINKF), A(KIF2EL), A(KIF2E2),
     &      LENF, LINKF, IF2EL, IF2EL2)

         CALL MDDEL ('XLENF')
         CALL MDDEL ('XLINKF')
         CALL MDDEL ('XIF2EL')
         IF (IS3DIM) CALL MDDEL ('XIF2E2')
         CALL MDSTAT (NERR, MEM)
         IF (NERR .GT. 0) GOTO 130
      END IF
  130 CONTINUE
      RETURN
      END
