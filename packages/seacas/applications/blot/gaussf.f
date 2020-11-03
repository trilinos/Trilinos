C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE GAUSSF (IS3DIM, NLNKF, LINKF1, XN, YN, ZN,
     &   XGAUSS, YGAUSS, ZGAUSS)
C=======================================================================

C   --*** GAUSSF *** (DETOUR) Calculate gauss point for a face
C   --   Written by Amy Gilkey - revised 10/20/87
C   --
C   --GAUSSF calculates the 4 element gauss points given the coordinates
C   --of the 4 nodes.
C   --
C   --Parameters:
C   --   IS3DIM - IN - true iff 3D versus 2D
C   --   NLNKF - IN - the number of nodes per face; must be 4
C   --   LINKF1 - IN - the connectivity for the face
C   --   XN, YN, ZN - IN - the nodal coordinates
C   --   XGAUSS, YGAUSS, ZGAUSS - OUT - the gauss point coordinates

      PARAMETER (S = .577350269189626)

      LOGICAL IS3DIM
      INTEGER LINKF1(NLNKF)
      REAL XN(*), YN(*), ZN(*)
      REAL XGAUSS(4), YGAUSS(4), ZGAUSS(4)

      REAL WT(4,4)
      LOGICAL FIRST
      SAVE FIRST, WT
C      --FIRST - true iff first time in routine
C      --WT - the gauss point constants

      DATA FIRST /.TRUE./

      IF (FIRST) THEN
         FIRST = .FALSE.

C      --Determine gauss point constants

         DO 100 N = 1, 4
            IF (N .EQ. 1) THEN
               SI = -S
               TI = -S
            ELSE IF (N .EQ. 2) THEN
               SI = +S
               TI = -S
            ELSE IF (N .EQ. 3) THEN
               SI = +S
               TI = +S
            ELSE IF (N .EQ. 4) THEN
               SI = -S
               TI = +S
            END IF
            WT(1,N) = .25 * (1-SI) * (1-TI)
            WT(2,N) = .25 * (1+SI) * (1-TI)
            WT(3,N) = .25 * (1+SI) * (1+TI)
            WT(4,N) = .25 * (1-SI) * (1+TI)
  100    CONTINUE
      END IF

C   --Determine gauss point coordinates

      DO 110 N = 1, 4
         XGAUSS(N) = WT(1,N)*XN(LINKF1(1)) + WT(2,N)*XN(LINKF1(2))
     &      + WT(3,N)*XN(LINKF1(3)) + WT(4,N)*XN(LINKF1(4))
         YGAUSS(N) = WT(1,N)*YN(LINKF1(1)) + WT(2,N)*YN(LINKF1(2))
     &      + WT(3,N)*YN(LINKF1(3)) + WT(4,N)*YN(LINKF1(4))
         IF (IS3DIM) THEN
            ZGAUSS(N) = WT(1,N)*ZN(LINKF1(1)) + WT(2,N)*ZN(LINKF1(2))
     &         + WT(3,N)*ZN(LINKF1(3)) + WT(4,N)*ZN(LINKF1(4))
         END IF
  110 CONTINUE

      RETURN
      END
