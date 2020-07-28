C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE FNDE2F (IEL, LENF, IF2EL, NQARY, IFACES, IELB)
C=======================================================================

C   --*** FNDE2F *** (MESH) Find faces that make up element
C   --   Written by Amy Gilkey - revised 07/06/87
C   --
C   --FNDE2F finds the surface faces that make up the given element.
C   --
C   --Parameters:
C   --   IEL - IN - the element number
C   --   LENF - IN - the cumulative face counts by element block
C   --   IF2EL - IN - the element number of each face
C   --   NQARY - IN/OUT - input as the maximum length of the IFACES array;
C   --      output as the length of the IFACES array
C   --   IFACES - OUT - the surface faces that make up the element
C   --   IELB - OUT - the element block of the element quarilaterals
C   --
C   --Common Variables:
C   --   Uses NDIM, NELBLK of /DBNUMS/

      include 'dbnums.blk'

      INTEGER LENF(0:NELBLK)
      INTEGER IF2EL(*)
      INTEGER IFACES(*)

      NQARY = 0
      DO 110 IELB = 1, NELBLK
         DO 100 IFAC = LENF(IELB-1)+1, LENF(IELB)
            IF (IEL .EQ. IF2EL(IFAC)) THEN
               NQARY = NQARY + 1
               IFACES(NQARY) = IFAC
            END IF
  100    CONTINUE
         IF (NQARY .GT. 0) GOTO 120
  110 CONTINUE
      IELB = 1

  120 CONTINUE
      RETURN
      END
